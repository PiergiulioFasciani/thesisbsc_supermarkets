#!/usr/bin/env Rscript

# ==============================================================================
# STACKED PPML EVENT-STUDY ANALYSIS
# ==============================================================================
# This script implements the public, methodology-backed thesis pipeline:
# 1. Load configuration and panel data
# 2. Apply time-window filtering and aggregation
# 3. Build stacked cohort samples using not-yet-treated controls
# 4. Estimate PPML event studies with site-clustered standard errors
# 5. Run Wald pretrend tests, support-weighted post-treatment ATT, and HonestDiD
# 6. Export tables, plots, and a compact results file
# ======================================================================

suppressPackageStartupMessages({
  library(fixest)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(yaml)
  library(zoo)
  library(car)
  library(HonestDiD)
})

cm_roman_family <- "CM Roman"

register_cm_roman_font <- function() {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  font_dir <- "/usr/local/texlive/2025basic/texmf-dist/fonts/opentype/public/lm"
  font_files <- c(
    plain = file.path(font_dir, "lmroman10-regular.otf"),
    bold = file.path(font_dir, "lmroman10-bold.otf"),
    italic = file.path(font_dir, "lmroman10-italic.otf"),
    bolditalic = file.path(font_dir, "lmroman10-bolditalic.otf")
  )

  if (!all(file.exists(font_files))) {
    return(invisible(FALSE))
  }

  tryCatch(
    {
      systemfonts::register_font(
        name = cm_roman_family,
        plain = unname(font_files["plain"]),
        bold = unname(font_files["bold"]),
        italic = unname(font_files["italic"]),
        bolditalic = unname(font_files["bolditalic"])
      )
      invisible(TRUE)
    },
    error = function(e) invisible(FALSE)
  )
}

register_cm_roman_font()

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

sanitize_filename <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

theme_aer <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = cm_roman_family, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0, family = cm_roman_family),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0, family = cm_roman_family),
      axis.title = element_text(face = "bold", size = 11, family = cm_roman_family),
      plot.caption = element_text(size = 8.5, color = "grey35", hjust = 0, family = cm_roman_family),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

event_study_title <- function(outcome, shift_q = 0) {
  base_title <- switch(
    outcome,
    count = "Effect of Supermarket Entry on Retail Density",
    entries = "Effect of Supermarket Entry on Retail Entry",
    exits = "Effect of Supermarket Entry on Retail Exit",
    paste("Effect of Supermarket Entry on", outcome)
  )

  if (shift_q > 0) {
    paste0(base_title, " (Six-Quarter Anticipation Shift)")
  } else {
    base_title
  }
}

run_label_for_shift <- function(shift_q) {
  if (shift_q > 0) paste0("shift", shift_q, "q") else "baseline"
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("  STACKED PPML EVENT-STUDY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

config_path <- "config.yml"
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

cat("Loading configuration from:", config_path, "\n")
config <- yaml::read_yaml(config_path)

analysis_config <- config$analysis %||% config

# workflow-level configs for metadata and downstream reporting
panel_config <- config$panel %||% list()
sampler_config <- config$sampler %||% list()

# Anticipation settings: whether to run the shifted (anticipation) robustness
# analysis and by how many aggregated periods. Defaults: enabled = TRUE, periods = 6
anticipation_enabled <- if (is.null(analysis_config$anticipation) || is.null(analysis_config$anticipation$enabled)) TRUE else as.logical(analysis_config$anticipation$enabled)
anticipation_periods <- as.integer(analysis_config$anticipation$periods %||% 6)

time_aggregation <- tolower(analysis_config$time_aggregation %||% "quarterly")
estimator <- tolower(analysis_config$estimator %||% "manual")
pre_periods <- analysis_config$pre_periods %||% 12
post_periods <- analysis_config$post_periods %||% 12
cluster_var <- analysis_config$cluster_var %||% "site_id"

min_year <- if (!is.null(analysis_config$start_date) && analysis_config$start_date != "") {
  as.numeric(format(as.Date(analysis_config$start_date), "%Y"))
} else {
  NULL
}

max_year <- if (!is.null(analysis_config$end_date) && analysis_config$end_date != "") {
  as.numeric(format(as.Date(analysis_config$end_date), "%Y"))
} else {
  NULL
}

pretrend_periods <- NULL
honestdid_pre_periods <- NULL
configured_anticipation_shift_q <- 6   # quarters (set to 0 to skip robustness run)

cat("Configuration:\n")
cat("  Estimator:", toupper(estimator), "(stacked cohort design)\n")
cat("  Time aggregation:", toupper(time_aggregation), "\n")
cat("  Time window:", sprintf("-%d to +%d %s\n", pre_periods, post_periods, time_aggregation))
if (!is.null(min_year)) cat("  Min year:", min_year, "\n")
if (!is.null(max_year)) cat("  Max year:", max_year, "\n")
cat("  Cluster variable:", cluster_var, "\n")
cat(sprintf("  Anticipation enabled: %s\n", ifelse(anticipation_enabled, "TRUE", "FALSE")))
cat(sprintf("  Anticipation periods (in %s): %d\n", time_aggregation, anticipation_periods))

# Print sampler / panel metadata for reproducibility
cat("\n  CONFIG METADATA (from config.yml):\n")
cat("    Sampler input_file:", sampler_config$input_file %||% "", "\n")
cat("    Sampler AOI center (lat,lon):", sampler_config$center_lat %||% "", sampler_config$center_lon %||% "", "\n")
cat("    Sampler radius_meters:", sampler_config$radius_meters %||% "", "\n")
cat("    Panel treatment_type:", panel_config$treatment_type %||% "", "\n")
cat("    Panel aot_radius_m:", panel_config$aot_radius_m %||% "", "\n")
cat("    Panel outcome_types:", if (!is.null(panel_config$outcome_types)) paste(panel_config$outcome_types, collapse = ", ") else "", "\n\n")
cat("  Pretrend test periods: full available pre-periods\n\n")

cat("Loading panel data...\n")
if (is.null(analysis_config$input_file) || analysis_config$input_file == "") {
  panel_files <- Sys.glob("data/output/panels/panel_*/panel.csv")
  if (length(panel_files) == 0) {
    stop("No panel directories found in data/output/panels/")
  }
  panel_info <- file.info(panel_files)
  latest_panel_file <- rownames(panel_info)[order(panel_info$mtime, decreasing = TRUE)][1]
  input_file <- latest_panel_file
} else {
  input_file <- analysis_config$input_file
}

if (!file.exists(input_file)) {
  stop("Panel file not found: ", input_file)
}

cat("Input:", input_file, "\n")
panel <- fread(input_file)

expected_cols <- c(
  "site_id", "month", "poi_poc", "count", "entries", "exits",
  "date_of_treatment", "cohort_id", "time_to_treatment"
)

missing_cols <- setdiff(expected_cols, names(panel))
extra_cols <- setdiff(names(panel), expected_cols)

if (length(missing_cols) > 0) {
  stop("Panel is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (length(extra_cols) > 0) {
  warning("Panel contains extra columns not used by the methodology: ",
          paste(extra_cols, collapse = ", "))
}

panel$month <- as.Date(paste0(panel$month, "-01"))

cat("Panel dimensions (before filtering):", nrow(panel), "rows x", ncol(panel), "columns\n")
cat("Sites:", length(unique(panel$site_id)), "\n")
cat("Months:", length(unique(panel$month)), "\n")
cat("Date range:", as.character(min(panel$month)), "to", as.character(max(panel$month)), "\n")
cat("Cohorts:", length(unique(panel$cohort_id)), "\n\n")

if (!is.null(min_year) || !is.null(max_year)) {
  cat("Applying time-window filter...\n")
  original_nrow <- nrow(panel)
  panel$year <- as.numeric(format(panel$month, "%Y"))

  if (!is.null(min_year)) {
    panel <- panel[year >= min_year]
    cat("  Filtered to year >=", min_year, "\n")
  }
  if (!is.null(max_year)) {
    panel <- panel[year <= max_year]
    cat("  Filtered to year <=", max_year, "\n")
  }

  cat("  Observations:", original_nrow, "->", nrow(panel),
      sprintf("(%.1f%% retained)\n", 100 * nrow(panel) / original_nrow))
  cat("  New date range:", as.character(min(panel$month)), "to", as.character(max(panel$month)), "\n\n")
}

aggregate_panel <- function(panel, time_aggregation) {
  cat("Converting to", toupper(time_aggregation), "aggregation...\n")

  if (time_aggregation == "monthly") {
    panel <- panel %>%
      mutate(
        period_index = as.integer(format(month, "%Y")) * 12 + as.integer(format(month, "%m")),
        period_to_treatment = time_to_treatment,
        period_date = month
      )
    time_label <- "Months"
  } else if (time_aggregation == "quarterly") {
    panel <- panel %>%
      mutate(
        month_num = as.integer(format(month, "%m")),
        period_index = as.integer(format(month, "%Y")) * 4 + ((month_num - 1) %/% 3) + 1,
        quarter = as.yearqtr(month),
        quarter_start = as.Date(quarter),
        period_to_treatment = floor(time_to_treatment / 3),
        period_date = quarter_start
      )
    time_label <- "Quarters"
  } else if (time_aggregation == "yearly") {
    panel <- panel %>%
      mutate(
        period_index = as.integer(format(month, "%Y")),
        year = as.numeric(format(month, "%Y")),
        year_start = as.Date(paste0(year, "-01-01")),
        treat_year = as.numeric(format(as.Date(paste0(date_of_treatment, "-01")), "%Y")),
        period_to_treatment = year - treat_year,
        period_date = year_start
      )
    time_label <- "Years"
  } else {
    stop("Unsupported time aggregation: ", time_aggregation)
  }

  cat("Aggregating per site-period...\n")
  if (time_aggregation == "quarterly") {
    panel_agg <- panel %>%
      group_by(site_id, date_of_treatment, cohort_id, period_to_treatment, period_date, period_index, poi_poc) %>%
      summarise(
        count = mean(count, na.rm = TRUE),
        entries = sum(entries, na.rm = TRUE),
        exits = sum(exits, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (time_aggregation == "yearly") {
    panel_agg <- panel %>%
      group_by(site_id, date_of_treatment, cohort_id, period_to_treatment, period_date, period_index, poi_poc) %>%
      summarise(
        count = mean(count, na.rm = TRUE),
        entries = sum(entries, na.rm = TRUE),
        exits = sum(exits, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    panel_agg <- panel %>%
      group_by(site_id, date_of_treatment, cohort_id, period_to_treatment, period_date, period_index, poi_poc) %>%
      summarise(
        count = sum(count, na.rm = TRUE),
        entries = sum(entries, na.rm = TRUE),
        exits = sum(exits, na.rm = TRUE),
        .groups = "drop"
      )
  }

  panel_agg <- panel_agg %>%
    mutate(post = ifelse(period_to_treatment >= 1, 1, 0))

  panel_wide <- panel_agg %>%
    pivot_wider(
      names_from = poi_poc,
      values_from = c(count, entries, exits),
      names_glue = "{poi_poc}_{.value}"
    ) %>%
    rename(
      POI_count = POI_count,
      POC_count = POC_count,
      POI_entries = POI_entries,
      POC_entries = POC_entries,
      POI_exits = POI_exits,
      POC_exits = POC_exits,
      month = period_date
    ) %>%
    rename(time_to_treatment_period = period_to_treatment)

  panel_wide$time_period_label <- format(panel_wide$month, "%Y-%m")

  list(panel_agg = panel_wide, time_label = time_label)
}

panel_info <- aggregate_panel(panel, time_aggregation)
panel_agg <- panel_info$panel_agg
time_label <- panel_info$time_label

cat("Aggregated to", nrow(panel_agg), "site-period observations\n\n")
cat("Reshaping to long format...\n")

panel_long <- panel_agg %>%
  pivot_longer(
    cols = c(POI_count, POI_entries, POI_exits, POC_count, POC_entries, POC_exits),
    names_to = "outcome_var",
    values_to = "outcome_value"
  ) %>%
  mutate(
    is_poi = ifelse(grepl("^POI_", outcome_var), 1, 0),
    outcome_type = gsub("^(POI|POC)_", "", outcome_var),
    is_treated_site = ifelse(!is.na(date_of_treatment) & !is.na(cohort_id) & cohort_id != "" & cohort_id != "never_treated" & cohort_id != "0", 1, 0),
    treated_post = is_poi * post,
    cohort_year = ifelse(is_treated_site == 1 & !is.na(cohort_id) & cohort_id != "", as.numeric(substr(cohort_id, 1, 4)), NA_real_),
    cohort_quarter = ifelse(is_treated_site == 1 & !is.na(cohort_id) & cohort_id != "", as.numeric(gsub(".*Q", "", cohort_id)), NA_real_),
    cohort_period = ifelse(!is.na(cohort_year) & !is.na(cohort_quarter), cohort_year * 4 + cohort_quarter, NA_real_)
  ) %>%
  group_by(site_id) %>%
  mutate(
    time_index = as.numeric(month - min(month)) / 365.25,
    treat_cohort = ifelse(any(is_treated_site == 1) & any(!is.na(date_of_treatment)),
                          as.Date(paste0(first(na.omit(date_of_treatment)), "-01")),
                          as.Date(NA))
  ) %>%
  ungroup() %>%
  mutate(
    cohort_id_numeric = as.numeric(factor(cohort_id)),
    cohort_trend = cohort_id_numeric * time_index,
    year_month = format(month, "%Y-%m"),
    treat_cohort = as.Date(treat_cohort, origin = "1970-01-01")
  )

cat("Reshaped to long format:", nrow(panel_long), "observations\n")
cat("  Treated sites (unique):", length(unique(panel_long$site_id[panel_long$is_treated_site == 1])), "\n")
cat("  Control sites (unique):", length(unique(panel_long$site_id[panel_long$is_treated_site == 0])), "\n")

n_control_sites <- length(unique(panel_long$site_id[panel_long$is_treated_site == 0]))
n_treated_sites <- length(unique(panel_long$site_id[panel_long$is_treated_site == 1]))

if (n_control_sites == 0) {
  cat("\n⚠️  WARNING: No never-treated sites detected. The stacked estimator will rely on not-yet-treated controls.\n\n")
  cohort_periods <- panel_long %>%
    filter(is_treated_site == 1, !is.na(cohort_period)) %>%
    pull(cohort_period) %>%
    unique() %>%
    sort()
  if (length(cohort_periods) < 3) {
    stop("Insufficient variation for stacked estimator - need never-treated sites or more cohort variation")
  }
  cat("  Found", length(cohort_periods), "treatment cohorts - proceeding with stacked design\n\n")
}

run_stacked_ppml_analysis <- function(shift_q, base_output_dir = NULL, base_plots_dir = NULL, base_csv_dir = NULL, consolidated_sink = FALSE) {
  cat("Anticipation shift (quarters):", shift_q, "\n")
  cat("Reference period (shifted): -1\n")
  cat("Reference period (original):", -(1 + shift_q), "\n")

  shifted_k_min <- -pre_periods
  shifted_k_max <- post_periods + shift_q
  window_label <- sprintf("%d to +%d", shifted_k_min, shifted_k_max)
  specification_label <- if (shift_q > 0) paste("Shift", shift_q, "quarters") else "Standard stacked PPML"
  if (shift_q > 0) {
    cat("Observation window (original):", sprintf("-%d to +%d", pre_periods, post_periods), "\n")
    cat("Observation window (shifted):", window_label, "\n")
  }

  if (is.null(base_output_dir)) {
    suffix <- "_baseline"
    if (!is.null(min_year) || !is.null(max_year)) {
      if (!is.null(min_year) && !is.null(max_year)) {
        suffix <- paste0(suffix, sprintf("_%d_%d", min_year, max_year))
      } else if (!is.null(min_year)) {
        suffix <- paste0(suffix, sprintf("_%d_onwards", min_year))
      } else {
        suffix <- paste0(suffix, sprintf("_up_to_%d", max_year))
      }
    }
    if (time_aggregation != "quarterly") {
      suffix <- paste0(suffix, "_", time_aggregation)
    }
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("data/output/analysis", paste0("stacked_ppml_es_", timestamp, suffix))
  } else {
    output_dir <- base_output_dir
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(base_plots_dir)) {
    plots_dir <- file.path(output_dir, "plots")
  } else {
    plots_dir <- base_plots_dir
  }
  
  if (is.null(base_csv_dir)) {
    csv_dir <- file.path(output_dir, "csv")
  } else {
    csv_dir <- base_csv_dir
  }
  
  plots_subfolder <- if (shift_q > 0) "anticipation" else "standard"
  csv_suffix <- if (shift_q > 0) paste0("_shift", shift_q, "q") else "_baseline"
  run_label <- run_label_for_shift(shift_q)
  dir.create(file.path(plots_dir, plots_subfolder), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(csv_dir, "standard"), recursive = TRUE, showWarnings = FALSE)

  cat("Output directory:", output_dir, "\n")
  cat("  Plots:", plots_dir, "\n")
  cat("  CSV exports:", csv_dir, "\n\n")

  results_file <- file.path(output_dir, "stacked_ppml_event_study_results.txt")
  
  if (!consolidated_sink) {
    sink(results_file)
    on.exit({ while (sink.number() > 0) sink() }, add = TRUE)
  }

  if (!consolidated_sink) {
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  STACKED PPML EVENT-STUDY ANALYSIS\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Input:", input_file, "\n")
    cat("Output:", output_dir, "\n")
    cat("Data aggregation:", toupper(time_aggregation), "\n")
    cat("Cluster variable:", cluster_var, "\n")
    cat("Pretrend test periods: full available pre-periods\n\n")
    cat("Estimation event-time window:", window_label, "\n\n")

    cat("Specification:\n")
    cat("  - is_poi = 1 for POI rows, 0 for POC rows\n")
    cat("  - is_treated_site = 1 for treated sites, 0 for never-treated or not-yet-treated controls\n")
    cat("  - Event-time coefficients interact relative time with treated_poi = treated_in_stack × is_poi\n")
    cat("  - Stacked sample uses not-yet-treated sites as controls\n\n")
  } else {
    cat("\n───────────────────────────────────────────────────────────────\n")
    cat("ROBUSTNESS: ANTICIPATION SHIFT ANALYSIS (shift_q = ", shift_q, " quarters)\n", sep = "")
    if (shift_q > 0) {
      cat("  Original window baseline: -", pre_periods, " to +", post_periods, "\n", sep = "")
      cat("  Added pre-treatment periods for robustness: +", shift_q, "\n", sep = "")
      cat("  Shifted window used for estimation: ", window_label, "\n", sep = "")
    }
    cat("───────────────────────────────────────────────────────────────\n\n")
  }

build_stacked_dataset <- function(panel_long, outcome_type, k_min, k_max, anticipation_shift_q = 0) {
  cat("Building stacked dataset for outcome:", outcome_type, "\n")
  df <- panel_long %>% filter(outcome_type == !!outcome_type)
  treated_cohorts <- df %>% filter(is_treated_site == 1, !is.na(cohort_period)) %>% pull(cohort_period) %>% unique() %>% sort()

  if (length(treated_cohorts) == 0) {
    stop("No treatment cohorts found for outcome: ", outcome_type)
  }

  never_treated_sites <- df %>% filter(is_treated_site == 0) %>% pull(site_id) %>% unique()
  cat("  Treated cohorts found:", length(treated_cohorts), "\n")
  cat("  Never-treated sites:", length(never_treated_sites), "\n")

  stack_list <- list()
  skipped_cohorts <- character()

  for (g in treated_cohorts) {
    treated_sites <- df %>% filter(is_treated_site == 1, cohort_period == g) %>% pull(site_id) %>% unique()
    control_sites <- df %>% filter(is_treated_site == 0 | cohort_period > g) %>% pull(site_id) %>% unique()
    stack_sites <- unique(c(treated_sites, control_sites))

    if (length(treated_sites) == 0 || length(stack_sites) == 0) {
      skipped_cohorts <- c(skipped_cohorts, as.character(g))
      next
    }

    stack_df <- df %>% filter(site_id %in% stack_sites) %>% mutate(stack_g = g)
    stack_df <- stack_df %>% mutate(
      treated_in_stack = ifelse(site_id %in% treated_sites, 1, 0),
      treated_poi = treated_in_stack * is_poi,
      rel_time = period_index - g,
      rel_time_shifted = rel_time + anticipation_shift_q
    ) %>% filter(rel_time_shifted >= k_min, rel_time_shifted <= k_max)

    if (nrow(stack_df) == 0) {
      skipped_cohorts <- c(skipped_cohorts, as.character(g))
      next
    }

    if (!any(stack_df$treated_in_stack == 1 & stack_df$rel_time_shifted == -1, na.rm = TRUE)) {
      skipped_cohorts <- c(skipped_cohorts, as.character(g))
      next
    }

    stack_df <- stack_df %>%
      mutate(
        rel_time_factor = factor(rel_time_shifted, levels = c(setdiff(as.character(k_min:k_max), "-1"), "-1")),
        rel_time_factor = relevel(rel_time_factor, ref = "-1"),
        site_stack = interaction(stack_g, site_id, drop = TRUE),
        time_stack = interaction(stack_g, month, drop = TRUE)
      )

    stack_list[[length(stack_list) + 1]] <- stack_df
  }

  if (length(stack_list) == 0) {
    stop("No valid stacks created for outcome: ", outcome_type)
  }

  stacked <- bind_rows(stack_list)
  stacked$rel_time_factor <- factor(stacked$rel_time_shifted, levels = c(setdiff(as.character(k_min:k_max), "-1"), "-1"))
  stacked$rel_time_factor <- relevel(stacked$rel_time_factor, ref = "-1")

  if (!"treated_poi" %in% names(stacked)) stop("treated_poi was not created")
  if (sum(stacked$treated_poi, na.rm = TRUE) == 0) stop("treated_poi has no treated POI observations")

  cat("  Stacked dataset created:", nrow(stacked), "observations\n")
  cat("  Stacks used:", length(unique(stacked$stack_g)), "\n")
  cat("  Treated sites:", n_distinct(stacked$site_id[stacked$treated_in_stack == 1]), "\n")
  cat("  Control sites:", n_distinct(stacked$site_id[stacked$treated_in_stack == 0]), "\n")
  cat("  treated_poi rows:", sum(stacked$treated_poi == 1, na.rm = TRUE), "\n")
  cat("  treated_poi sites:", n_distinct(stacked$site_id[stacked$treated_poi == 1]), "\n")
  cat("  Original rel-time range:", paste(range(stacked$rel_time, na.rm = TRUE), collapse = " to "), "\n")
  cat("  Shifted rel-time range:", paste(range(stacked$rel_time_shifted, na.rm = TRUE), collapse = " to "), time_label, "\n\n")

  stacked
}

estimate_stacked_es <- function(stacked_df, ref = -1, method = "ppml", cluster_var = "site_id") {
  cat("Estimating stacked event study using", toupper(method), "...\n")
  vcov_formula <- as.formula(paste0("~", cluster_var))
  ref_char <- as.character(ref)

  model <- if (method == "ppml") {
    fepois(
      outcome_value ~ is_poi + treated_poi + i(rel_time_factor, treated_poi, ref = ref_char) |
        site_stack + time_stack,
      data = stacked_df,
      vcov = vcov_formula
    )
  } else if (method == "ols") {
    feols(
      outcome_value ~ is_poi + treated_poi + i(rel_time_factor, treated_poi, ref = ref_char) |
        site_stack + time_stack,
      data = stacked_df,
      vcov = vcov_formula
    )
  } else {
    stop("Method must be 'ppml' or 'ols'")
  }

  cat("  Model estimated successfully\n")
  coef_names <- names(coef(model))
  interaction_coefs <- coef_names[grepl("rel_time_factor::", coef_names) & grepl(":treated_poi", coef_names)]
  if (length(interaction_coefs) > 0) {
    times_estimated <- sort(as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", interaction_coefs)))
    cat("  Event-time coefficients identified:", length(times_estimated), "\n")
    cat("  Range:", min(times_estimated), "to", max(times_estimated), "\n\n")
  } else {
    cat("  No event-time coefficients identified\n\n")
  }
  model
}

extract_stacked_es_coefs <- function(model, k_min = -12, k_max = 12) {
  coef_data <- coef(model)
  se_data <- se(model)
  time_vars <- names(coef_data)[grepl("rel_time_factor::", names(coef_data)) & grepl(":treated_poi", names(coef_data))]

  full_grid <- data.frame(time = k_min:k_max)
  if (length(time_vars) == 0) {
    warning("No interaction coefficients found")
    return(full_grid %>% mutate(coef = NA_real_, se = NA_real_, pct_effect = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, is_significant = FALSE))
  }

  times_estimated <- as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", time_vars))
  coef_est <- data.frame(
    time = times_estimated,
    coef = as.numeric(coef_data[time_vars]),
    se = as.numeric(se_data[time_vars]),
    row.names = NULL
  )

  coef_df <- full_grid %>% left_join(coef_est, by = "time") %>%
    mutate(
      coef = ifelse(time == -1, 0, coef),
      se = ifelse(time == -1, 0, se),
      pct_effect = ifelse(is.na(coef), NA_real_, (exp(coef) - 1) * 100),
      ci_lower = ifelse(is.na(coef) | is.na(se), NA_real_, (exp(coef - 1.96 * se) - 1) * 100),
      ci_upper = ifelse(is.na(coef) | is.na(se), NA_real_, (exp(coef + 1.96 * se) - 1) * 100),
      is_significant = ifelse(is.na(ci_lower) | is.na(ci_upper), FALSE, !(ci_lower <= 0 & ci_upper >= 0))
    )

  expected_times <- setdiff(k_min:k_max, -1)
  missing_times <- setdiff(expected_times, times_estimated)
  if (length(missing_times) > 0) {
    cat("  WARNING: Missing coefficients for event-times:", paste(head(missing_times, 10), collapse = ", "))
    if (length(missing_times) > 10) cat(", ... (", length(missing_times), "total)", sep = "")
    cat("\n")
  }

  coef_df
}

get_treated_support_by_event_time <- function(stacked_df, k_min = NULL, k_max = NULL) {
  df <- stacked_df %>% filter(treated_in_stack == 1)
  time_var <- if ("rel_time_shifted" %in% names(df)) "rel_time_shifted" else "rel_time"
  df <- df %>% mutate(time = .data[[time_var]])
  if (!is.null(k_min)) df <- df %>% filter(.data[[time_var]] >= k_min)
  if (!is.null(k_max)) df <- df %>% filter(.data[[time_var]] <= k_max)

  df %>%
    group_by(time) %>%
    summarise(
      n_treated_obs = n(),
      n_stacks = n_distinct(stack_g),
      .groups = "drop"
    )
}

aggregate_event_time_att <- function(coef_df, stacked_df, k_min = NULL, k_max = NULL) {
  support_df <- get_treated_support_by_event_time(stacked_df, k_min = k_min, k_max = k_max)

  agg_df <- coef_df %>%
    filter(!is.na(pct_effect), time != -1) %>%
    left_join(support_df, by = "time") %>%
    mutate(
      n_treated_obs = ifelse(is.na(n_treated_obs), 0, n_treated_obs),
      n_stacks = ifelse(is.na(n_stacks), 0, n_stacks)
    ) %>%
    filter(n_treated_obs > 0) %>%
    mutate(weight_event_time = n_treated_obs / sum(n_treated_obs, na.rm = TRUE))

  agg_df
}

compute_weighted_post_att <- function(model, stacked_df, start_period = 1, end_period = NULL) {
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  coef_names <- names(coefs)
  interaction_terms <- grep("rel_time_factor::[-]?[0-9]+:treated_poi", coef_names, value = TRUE)

  if (length(interaction_terms) == 0) {
    return(list(
      att_pct = NA_real_, theta_log = NA_real_, n_periods = 0, periods = integer(), total_treated_weight = 0,
      weights = data.frame(), l_vector = NULL, coef_names = character(), att_se_pct = NA_real_,
      att_ci_lower = NA_real_, att_ci_upper = NA_real_, estimand_label = "Support-weighted post-treatment PPML effect"
    ))
  }

  times <- as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", interaction_terms))
  keep_idx <- times >= start_period
  if (!is.null(end_period)) keep_idx <- keep_idx & times <= end_period

  times_kept <- times[keep_idx]
  coef_names_kept <- interaction_terms[keep_idx]
  if (length(times_kept) == 0) {
    return(list(
      att_pct = NA_real_, theta_log = NA_real_, n_periods = 0, periods = integer(), total_treated_weight = 0,
      weights = data.frame(), l_vector = NULL, coef_names = character(), att_se_pct = NA_real_,
      att_ci_lower = NA_real_, att_ci_upper = NA_real_, estimand_label = "Support-weighted post-treatment PPML effect"
    ))
  }

  support_df <- get_treated_support_by_event_time(stacked_df, k_min = start_period, k_max = end_period)
  weight_df <- data.frame(
    time = times_kept,
    coef_name = coef_names_kept,
    beta = as.numeric(coefs[coef_names_kept]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(support_df, by = "time") %>%
    mutate(
      n_treated_obs = ifelse(is.na(n_treated_obs), 0, n_treated_obs),
      n_stacks = ifelse(is.na(n_stacks), 0, n_stacks)
    ) %>%
    filter(n_treated_obs > 0)

  if (nrow(weight_df) == 0) {
    return(list(
      att_pct = NA_real_, theta_log = NA_real_, n_periods = 0, periods = integer(), total_treated_weight = 0,
      weights = data.frame(), l_vector = NULL, coef_names = character(), att_se_pct = NA_real_,
      att_ci_lower = NA_real_, att_ci_upper = NA_real_, estimand_label = "Support-weighted post-treatment PPML effect"
    ))
  }

  weight_df <- weight_df %>% mutate(weight = n_treated_obs / sum(n_treated_obs, na.rm = TRUE))
  w_vec <- weight_df$weight
  beta_vec <- weight_df$beta
  coef_names_final <- weight_df$coef_name
  vcov_sub <- vcov_mat[coef_names_final, coef_names_final, drop = FALSE]
  theta_log <- as.numeric(sum(w_vec * beta_vec))
  beta_var <- as.numeric(t(w_vec) %*% vcov_sub %*% w_vec)
  beta_se <- sqrt(beta_var)

  att_pct <- (exp(theta_log) - 1) * 100
  att_se_pct <- exp(theta_log) * 100 * beta_se
  att_ci_lower <- (exp(theta_log - 1.96 * beta_se) - 1) * 100
  att_ci_upper <- (exp(theta_log + 1.96 * beta_se) - 1) * 100

  list(
    att_pct = att_pct,
    theta_log = theta_log,
    n_periods = nrow(weight_df),
    periods = weight_df$time,
    total_treated_weight = sum(weight_df$n_treated_obs),
    weights = weight_df %>% select(time, coef_name, beta, n_treated_obs, weight),
    l_vector = w_vec,
    coef_names = coef_names_final,
    att_se_pct = att_se_pct,
    att_ci_lower = att_ci_lower,
    att_ci_upper = att_ci_upper,
    estimand_label = "Support-weighted post-treatment PPML effect"
  )
}

plot_stacked_es <- function(coef_df, title, output_path, time_label = "Quarters", k_min = -12, k_max = 12, shift_q = 0) {
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    warning("No coefficients to plot")
    return(NULL)
  }

  plot_df <- coef_df %>% filter(!is.na(coef))
  if (nrow(plot_df) == 0) {
    warning("No valid coefficients to plot")
    return(NULL)
  }

  p <- ggplot(plot_df, aes(x = time, y = pct_effect)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.45) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_line(color = "black", linewidth = 0.45) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), color = "black", width = 0.16, linewidth = 0.45) +
    geom_point(aes(fill = is_significant), shape = 21, color = "black", size = 1.9, stroke = 0.4) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    scale_x_continuous(limits = c(k_min, k_max), breaks = seq(k_min, k_max, by = 2)) +
    labs(
      x = paste0("Event Time (", time_label, ")"),
      y = "Percentage Effect"
    ) +
    theme_aer()

  if (shift_q > 0) {
    p <- p + geom_vline(xintercept = shift_q - 0.5, linetype = "dotted", color = "black", linewidth = 0.5)
  }

  ggsave(output_path, p, width = 7, height = 5, dpi = 300, bg = "white", device = ragg::agg_png)
  p
}

test_pretrends <- function(model, outcome_name, coef_df, pretrend_periods = NULL) {
  cat(sprintf("%s:\n", outcome_name))
  coef_names <- names(coef(model))
  available_pretrend_coefs <- coef_names[grepl("rel_time_factor::", coef_names) & grepl(":treated_poi", coef_names)]
  available_pretrend_periods <- sort(unique(as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", available_pretrend_coefs))))
  available_pretrend_periods <- available_pretrend_periods[available_pretrend_periods < 0]

  if (is.null(pretrend_periods)) {
    pretrend_periods <- available_pretrend_periods
  } else {
    pretrend_periods <- sort(unique(pretrend_periods[pretrend_periods < 0]))
    pretrend_periods <- intersect(available_pretrend_periods, pretrend_periods)
  }

  pre_coefs <- paste0("rel_time_factor::", pretrend_periods, ":treated_poi")
  pre_coefs <- pre_coefs[pre_coefs %in% coef_names]

  if (length(pre_coefs) == 0) {
    cat("  No pretrend coefficients found for periods:", paste(pretrend_periods, collapse = ", "), "\n\n")
    return(list(p_value = NA_real_, f_stat = NA_real_, periods_tested = 0, periods = integer(), periods_label = "none"))
  }

  cat(sprintf("  Pretrend test (periods %s):\n", paste(pretrend_periods, collapse = ", ")))
  wald_result <- tryCatch(wald(model, keep = pre_coefs), error = function(e) {
    cat("  Error in Wald test:", e$message, "\n")
    NULL
  })

  if (!is.null(wald_result)) {
    f_stat <- wald_result$stat
    p_val <- wald_result$p
    cat(sprintf("    F-statistic: %.2f\n", f_stat))
    cat(sprintf("    p-value: %.4f\n", p_val))
    if (!is.null(coef_df)) {
      lead_df <- coef_df %>% filter(time %in% pretrend_periods)
      if (nrow(lead_df) > 0) {
        cat(sprintf("    Max absolute lead estimate: %.3f%%\n", max(abs(lead_df$pct_effect), na.rm = TRUE)))
      }
    }
    if (p_val > 0.05) {
      cat("  No statistically significant pre-trend deviations detected (p > 0.05)\n")
    } else if (p_val > 0.01) {
      cat("  ⚠ Pre-trend deviations are marginally significant (0.01 < p < 0.05)\n")
    } else {
      cat("  ✗ Pre-trend deviations are statistically significant (p < 0.01)\n")
    }
    result <- list(
      p_value = p_val,
      f_stat = f_stat,
      periods_tested = length(pre_coefs),
      periods = pretrend_periods,
      periods_label = paste(pretrend_periods, collapse = ", ")
    )
  } else {
    result <- list(
      p_value = NA_real_,
      f_stat = NA_real_,
      periods_tested = length(pre_coefs),
      periods = pretrend_periods,
      periods_label = paste(pretrend_periods, collapse = ", ")
    )
  }

  cat("\n")
  result
}

run_honestdid_sensitivity <- function(model, outcome_name, output_dir,
                                      pre_periods_vec = NULL,
                                      post_periods_vec = 1:10,
                                      Mbar_grid = seq(0, 1, by = 0.1),
                                      subtitle = "",
                                      att_result = NULL,
                                      shift_q = 0) {
  cat("HonestDiD Sensitivity Analysis:", outcome_name, "\n")
  cat("  Subtitle:", subtitle, "\n")

  post_names <- paste0("rel_time_factor::", post_periods_vec, ":treated_poi")
  available_coefs <- names(coef(model))
  available_pre_names <- available_coefs[grepl("rel_time_factor::", available_coefs) & grepl(":treated_poi", available_coefs)]
  available_pre_periods <- sort(unique(as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", available_pre_names))))
  available_pre_periods <- available_pre_periods[available_pre_periods < 0]

  if (is.null(pre_periods_vec)) {
    pre_periods_vec <- available_pre_periods
  } else {
    pre_periods_vec <- sort(unique(pre_periods_vec[pre_periods_vec < 0]))
    pre_periods_vec <- intersect(available_pre_periods, pre_periods_vec)
  }

  pre_names <- paste0("rel_time_factor::", pre_periods_vec, ":treated_poi")
  pre_present <- pre_names[pre_names %in% available_coefs]
  post_present <- post_names[post_names %in% available_coefs]
  pre_missing <- pre_periods_vec[!pre_names %in% available_coefs]
  post_missing <- post_periods_vec[!post_names %in% available_coefs]

  if (length(pre_missing) > 0) cat("  NOTE: Missing pre-period coefficients for k =", paste(pre_missing, collapse = ", "), "\n")
  if (length(post_missing) > 0) cat("  NOTE: Missing post-period coefficients for k =", paste(post_missing, collapse = ", "), "\n")

  if (length(pre_present) < 2 || length(post_present) < 1) {
    cat("  ERROR: Insufficient coefficients for HonestDiD\n\n")
    return(invisible(list(results = NULL, breakdown_Mbar = NA_real_, theta_pct = NA_real_, theta_log = NA_real_, l_vec = NULL, pre_present = pre_present, post_present = post_present, restriction = "relative_magnitudes")))
  }

  cat("  Using", length(pre_present), "pre-periods and", length(post_present), "post-periods\n")

  if (!is.null(att_result)) {
    cat("\n  === ESTIMAND ALIGNMENT CHECK ===\n")
    cat("  Printed ATT uses periods:", paste(att_result$periods, collapse = ", "), "\n")
    cat("  HonestDiD post-periods:", paste(post_periods_vec, collapse = ", "), "\n")
    if (!is.na(att_result$att_ci_lower) && !is.na(att_result$att_ci_upper)) {
      includes_zero <- (att_result$att_ci_lower <= 0 && att_result$att_ci_upper >= 0)
      cat("  Conventional 95% CI: [", sprintf("%.3f", att_result$att_ci_lower), ",", sprintf("%.3f", att_result$att_ci_upper), "]\n", sep = "")
      cat("  Conventional CI includes 0:", includes_zero, "\n")
    }
    cat("  ================================\n\n")
  }

  all_names <- c(pre_present, post_present)
  betahat <- coef(model)[all_names]
  sigma <- vcov(model)[all_names, all_names, drop = FALSE]
  numPrePeriods <- length(pre_present)
  numPostPeriods <- length(post_present)
  post_times <- as.numeric(gsub(".*rel_time_factor::([-]?[0-9]+):.*", "\\1", post_present))

  if (!is.null(att_result) && !is.null(att_result$periods)) {
    if (!identical(sort(unique(att_result$periods)), sort(unique(post_times)))) {
      warning("HonestDiD post periods differ from att_result$periods")
    }
  }

  if (!is.null(att_result) && !is.null(att_result$weights)) {
    weight_map <- att_result$weights %>% select(time, weight)
    l_df <- data.frame(time = post_times, coef_name = post_present, stringsAsFactors = FALSE) %>%
      left_join(weight_map, by = "time") %>%
      mutate(weight = ifelse(is.na(weight), 0, weight))
    l_vec <- l_df$weight
    if (all(is.na(l_vec)) || sum(l_vec, na.rm = TRUE) == 0) {
      l_vec <- rep(1 / length(post_present), length(post_present))
    } else {
      l_vec <- l_vec / sum(l_vec, na.rm = TRUE)
    }
  } else {
    l_vec <- rep(1 / length(post_present), length(post_present))
  }

  if (length(l_vec) != numPostPeriods) {
    warning("length(l_vec) != numPostPeriods; reverting to equal weights")
    l_vec <- rep(1 / numPostPeriods, numPostPeriods)
  }
  if (abs(sum(l_vec) - 1) > 1e-8) {
    warning("sum(l_vec) does not equal 1; renormalizing")
    l_vec <- l_vec / sum(l_vec)
  }

  cat("  Pre-periods used:", paste(pre_present, collapse = ", "), "\n")
  cat("  Post-periods used:", paste(post_present, collapse = ", "), "\n")
  cat("  l_vec:", paste(sprintf("%.4f", l_vec), collapse = ", "), "\n")
  cat("  sum(l_vec):", sprintf("%.8f", sum(l_vec)), "\n")
  cat("  Mbar_grid:", paste(sprintf("%.3f", Mbar_grid), collapse = ", "), "\n")

  theta_log <- as.numeric(sum(l_vec * betahat[post_present]))
  theta_pct <- (exp(theta_log) - 1) * 100
  cat(sprintf("  Support-weighted PPML point estimate: %.2f%%\n", theta_pct))

  honest_ci <- tryCatch({
    createSensitivityResults_relativeMagnitudes(
      betahat = betahat,
      sigma = sigma,
      numPrePeriods = numPrePeriods,
      numPostPeriods = numPostPeriods,
      l_vec = l_vec,
      Mbarvec = Mbar_grid,
      method = "C-LF",
      bound = "deviation from parallel trends"
    )
  }, error = function(e) {
    cat("  Error in relative-magnitudes HonestDiD:", e$message, "\n\n")
    NULL
  })

  if (is.null(honest_ci)) {
    return(invisible(list(results = NULL, breakdown_Mbar = NA_real_, theta_pct = theta_pct, theta_log = theta_log, l_vec = l_vec, pre_present = pre_present, post_present = post_present, restriction = "relative_magnitudes")))
  }

  if (is.data.frame(honest_ci)) {
    honest_results <- honest_ci
  } else if (!is.null(honest_ci$results) && is.data.frame(honest_ci$results)) {
    honest_results <- honest_ci$results
  } else if (!is.null(honest_ci$lb) && !is.null(honest_ci$ub)) {
    honest_results <- data.frame(Mbar = Mbar_grid, lb = honest_ci$lb, ub = honest_ci$ub)
  } else {
    warning("createSensitivityResults_relativeMagnitudes() returned bounds without columns named lb and ub")
    return(invisible(list(results = NULL, breakdown_Mbar = NA_real_, theta_pct = theta_pct, theta_log = theta_log, l_vec = l_vec, pre_present = pre_present, post_present = post_present, restriction = "relative_magnitudes")))
  }

  if (!all(c("lb", "ub") %in% names(honest_results))) {
    warning("createSensitivityResults_relativeMagnitudes() returned bounds without columns named lb and ub")
    return(invisible(list(results = NULL, breakdown_Mbar = NA_real_, theta_pct = theta_pct, theta_log = theta_log, l_vec = l_vec, pre_present = pre_present, post_present = post_present, restriction = "relative_magnitudes")))
  }

  if (!"Mbar" %in% names(honest_results)) {
    if ("M" %in% names(honest_results)) {
      honest_results$Mbar <- honest_results$M
    } else {
      honest_results$Mbar <- Mbar_grid[seq_len(nrow(honest_results))]
    }
  }

  honest_results <- honest_results %>%
    mutate(
      lb_pct = (exp(lb) - 1) * 100,
      ub_pct = (exp(ub) - 1) * 100,
      includes_zero = lb_pct <= 0 & ub_pct >= 0
    )

  breakdown_candidates <- honest_results %>% filter(includes_zero)
  if (nrow(breakdown_candidates) > 0) {
    breakdown_Mbar <- min(breakdown_candidates$Mbar)
    cat(sprintf("  Breakdown Mbar: %.3f (smallest relaxation where CI includes 0)\n", breakdown_Mbar))
  } else {
    breakdown_Mbar <- NA_real_
    cat("  Breakdown Mbar: None (CI never includes 0 for tested Mbar values)\n")
  }
  min_M_zero <- breakdown_Mbar

  run_label <- run_label_for_shift(shift_q)
  plot_path <- file.path(output_dir, paste0(sanitize_filename(outcome_name), "_honestdid_", run_label, ".png"))
  p <- ggplot(honest_results, aes(x = Mbar)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.45) +
    geom_ribbon(aes(ymin = lb_pct, ymax = ub_pct), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = lb_pct), color = "black", linewidth = 0.65) +
    geom_line(aes(y = ub_pct), color = "black", linewidth = 0.65) +
    geom_hline(yintercept = theta_pct, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(
      x = expression(M),
      y = "Percentage Effect"
    ) +
    theme_aer()

  if (!is.na(min_M_zero)) {
    p <- p +
      geom_vline(xintercept = min_M_zero, linetype = "dotted", color = "black", linewidth = 0.5)
  }

  ggsave(plot_path, p, width = 7, height = 5, dpi = 300, bg = "white", device = ragg::agg_png)
  cat("  Plot saved:", plot_path, "\n\n")

  invisible(list(
    results = honest_results,
    breakdown_Mbar = breakdown_Mbar,
    min_M_zero = min_M_zero,
    theta_pct = theta_pct,
    theta_log = theta_log,
    l_vec = l_vec,
    pre_present = pre_present,
    post_present = post_present,
    restriction = "relative_magnitudes"
  ))
}

outcomes <- c("count", "entries", "exits")
standard_results <- list()
summary_rows <- list()

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PART 1: STANDARD STACKED EVENT STUDIES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

  for (out in outcomes) {
    if (!consolidated_sink) {
      cat("─────────────────────────────────────────────────────────────\n")
      cat("OUTCOME:", toupper(out), "\n")
      cat("─────────────────────────────────────────────────────────────\n\n")
    } else {
      cat("  OUTCOME:", toupper(out), "\n\n")
    }

    stacked_df <- build_stacked_dataset(
      panel_long,
      out,
      shifted_k_min,
      shifted_k_max,
      anticipation_shift_q = shift_q
    )
    model <- estimate_stacked_es(stacked_df, ref = -1, method = "ppml", cluster_var = cluster_var)

    if (!consolidated_sink) {
      cat("Model Summary:\n")
      print(summary(model))
      cat("\n")
    }

    coef_df <- extract_stacked_es_coefs(model, k_min = shifted_k_min, k_max = shifted_k_max)
    csv_path <- file.path(csv_dir, "standard", paste0(out, "_coefficients", csv_suffix, ".csv"))
    write.csv(coef_df, csv_path, row.names = FALSE)
    if (!consolidated_sink) cat("Coefficients saved to:", csv_path, "\n\n")

    if (!consolidated_sink) cat("PRETREND TEST:\n")
    pretrend_result <- test_pretrends(model, toupper(out), coef_df)

    plot_run_dir <- file.path(plots_dir, plots_subfolder)
    run_label <- run_label_for_shift(shift_q)
    plot_path <- file.path(plot_run_dir, paste0(sanitize_filename(out), "_event_study_", run_label, ".png"))
    plot_title <- event_study_title(out, shift_q)
    if (consolidated_sink) sink()
    plot_stacked_es(coef_df, plot_title, plot_path, time_label = time_label, k_min = shifted_k_min, k_max = shifted_k_max, shift_q = shift_q)
    if (consolidated_sink) sink(results_file, append = TRUE)
    if (!consolidated_sink) cat("Plot saved to:", plot_path, "\n\n")

    honest_post_max <- shifted_k_max
  att_result <- compute_weighted_post_att(model, stacked_df, start_period = 1, end_period = honest_post_max)
  cat(sprintf("%s (periods 1 to %d): ", att_result$estimand_label, honest_post_max))
  if (!is.na(att_result$att_pct)) {
    cat(sprintf("%.2f%% (SE: %.2f%%, 95%% CI: [%.2f%%, %.2f%%])\n", att_result$att_pct, att_result$att_se_pct, att_result$att_ci_lower, att_result$att_ci_upper))
  } else {
    cat("not available\n")
  }
  cat(sprintf("  (identified over %d post-treatment event-times; total treated obs weight = %d)\n\n", att_result$n_periods, att_result$total_treated_weight))

  event_time_att_df <- aggregate_event_time_att(coef_df, stacked_df, k_min = shifted_k_min, k_max = shifted_k_max)
  write.csv(event_time_att_df, file.path(csv_dir, "standard", paste0(out, "_event_time_support_weighted", csv_suffix, ".csv")), row.names = FALSE)

  cat("HONESTDID SENSITIVITY ANALYSIS:\n")
  if (consolidated_sink) sink()
  honestdid_result <- run_honestdid_sensitivity(
    model, toupper(out), file.path(plots_dir, plots_subfolder),
    post_periods_vec = 1:honest_post_max,
    subtitle = if (shift_q > 0) paste("Anticipation shift:", shift_q, "quarters | Window:", window_label) else paste("Standard stacked PPML event study | Window:", window_label),
    att_result = att_result,
    shift_q = shift_q
  )
  if (consolidated_sink) sink(results_file, append = TRUE)

  standard_results[[out]] <- list(
    model = model,
    coef_df = coef_df,
    pretrend = pretrend_result,
    att = att_result,
    honestdid = honestdid_result
  )

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    specification = specification_label,
    event_window = window_label,
    outcome = out,
    att_pct = att_result$att_pct,
    att_n_periods = att_result$n_periods,
    att_se_pct = att_result$att_se_pct,
    att_ci_lower = att_result$att_ci_lower,
    att_ci_upper = att_result$att_ci_upper,
    pretrend_f = pretrend_result$f_stat,
    pretrend_p = pretrend_result$p_value,
    pretrend_periods = pretrend_result$periods_label,
    honestdid_M0_lb_pct = if (!is.null(honestdid_result$results) && nrow(honestdid_result$results) > 0) honestdid_result$results$lb_pct[1] else NA_real_,
    honestdid_M0_ub_pct = if (!is.null(honestdid_result$results) && nrow(honestdid_result$results) > 0) honestdid_result$results$ub_pct[1] else NA_real_,
    breakdown_Mbar = honestdid_result$breakdown_Mbar,
    honestdid_min_M_zero = honestdid_result$min_M_zero,
    stringsAsFactors = FALSE
  )
}

  if (!consolidated_sink) {
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("PART 1 COMPLETE: Standard stacked PPML event studies\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")

    summary_table <- bind_rows(summary_rows)
    cat("SUMMARY TABLE:\n")
    print(summary_table, row.names = FALSE)
    cat("\n")

    summary_csv <- file.path(csv_dir, "complete_summary.csv")
    write.csv(summary_table, summary_csv, row.names = FALSE)
    cat("Summary table saved to:", summary_csv, "\n\n")

    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  KEY FINDINGS AND INTERPRETATION\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")

    summary_table <- bind_rows(summary_rows)
    for (i in seq_len(nrow(summary_table))) {
      row <- summary_table[i, ]
      cat(sprintf("─────────────────────────────────────────────────────────────\n"))
      cat(sprintf("%s - %s:\n", row$specification, toupper(row$outcome)))
      cat(sprintf("  Event window: %s\n", row$event_window))
      cat(sprintf("  ATT: %.2f%% (over %d periods)\n", row$att_pct, row$att_n_periods))
      if (!is.na(row$pretrend_p)) {
        pretrend_symbol <- ifelse(row$pretrend_p > 0.05, "✓", ifelse(row$pretrend_p > 0.01, "⚠", "✗"))
        cat(sprintf("  Pretrend test %s: F=%.2f, p=%.4f (periods %s)\n", pretrend_symbol, row$pretrend_f, row$pretrend_p, row$pretrend_periods))
      } else {
        cat("  Pretrend test: Not available\n")
      }

      if (!is.na(row$honestdid_M0_lb_pct) && !is.na(row$honestdid_M0_ub_pct)) {
        m0_includes_zero <- (row$honestdid_M0_lb_pct <= 0 & row$honestdid_M0_ub_pct >= 0)
        m0_symbol <- ifelse(!m0_includes_zero, "✓", "✗")
        cat(sprintf("  HonestDiD M=0 %s: [%.2f%%, %.2f%%]\n", m0_symbol, row$honestdid_M0_lb_pct, row$honestdid_M0_ub_pct))
        if (!is.na(row$honestdid_min_M_zero)) {
          cat(sprintf("  Minimum M for 0: %.3f\n", row$honestdid_min_M_zero))
        } else if (!is.na(row$breakdown_Mbar)) {
          cat(sprintf("  Minimum M for 0: %.3f\n", row$breakdown_Mbar))
        } else {
          cat("  Minimum M for 0: None within tested grid\n")
        }
      } else {
        cat("  HonestDiD: Not available\n")
      }

      pretrend_pass <- !is.na(row$pretrend_p) && row$pretrend_p > 0.05
      m0_pass <- !is.na(row$honestdid_M0_lb_pct) && !is.na(row$honestdid_M0_ub_pct) && !(row$honestdid_M0_lb_pct <= 0 & row$honestdid_M0_ub_pct >= 0)

      if (pretrend_pass && m0_pass) {
        cat("  Overall: ✓ ROBUST - passes pretrend and M=0 checks\n")
      } else if (!pretrend_pass && !m0_pass) {
        cat("  Overall: ✗ CAUTION - fails pretrend and M=0 checks\n")
      } else {
        cat("  Overall: ⚠ MIXED - inspect the validity checks individually\n")
      }
      cat("\n")
    }

    cat("═══════════════════════════════════════════════════════════════\n")
    cat("INTERPRETATION NOTES:\n")
    cat("  - Pretrend test: tests all available pre-treatment periods in the stacked model\n")
    cat("  - HonestDiD Mbar=0: robust CI under exact post-treatment parallel trends in the relative-magnitudes framework\n")
    cat("  - Minimum M for 0: smallest relative-magnitudes relaxation under which the robust CI includes zero\n")
    cat("  - Support-weighted post-treatment effect: weighted average of post-treatment PPML event-time effects\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")

    sink()

    cat("Stacked PPML event-study analysis complete.\n")
    cat("Output directory:", output_dir, "\n")
    cat("  - Results text:", results_file, "\n")
    cat("  - Summary table:", summary_csv, "\n")
    cat("  - Event-study plots:", file.path(plots_dir, "standard"), "\n")
    cat("  - HonestDiD plots:", file.path(plots_dir, "standard"), "\n")
    cat("  - CSV exports:", csv_dir, "\n")
  } else {
    return(invisible(list(summary_rows = summary_rows)))
  }

  invisible(list(
    output_dir = output_dir,
    results_file = results_file,
    plots_dir = plots_dir,
    csv_dir = csv_dir,
    summary_rows = summary_rows,
    shift_q = shift_q
  ))
}

baseline_run <- run_stacked_ppml_analysis(0)

if (anticipation_enabled && !is.na(anticipation_periods) && anticipation_periods > 0) {
  sink(baseline_run$results_file, append = TRUE)
  on.exit({ while (sink.number() > 0) sink() }, add = TRUE)
  
  shift_run <- run_stacked_ppml_analysis(
    anticipation_periods,
    base_output_dir = baseline_run$output_dir,
    base_plots_dir = baseline_run$plots_dir,
    base_csv_dir = baseline_run$csv_dir,
    consolidated_sink = TRUE
  )
  
  combined_summary <- bind_rows(baseline_run$summary_rows, shift_run$summary_rows)
  summary_csv <- file.path(baseline_run$csv_dir, "complete_summary.csv")
  write.csv(combined_summary, summary_csv, row.names = FALSE)
  cat("\nCombined summary table saved to:", summary_csv, "\n")

  comparison_df <- combined_summary %>%
    select(specification, outcome, att_pct, pretrend_p, breakdown_Mbar, honestdid_min_M_zero) %>%
    pivot_wider(
      names_from = specification,
      values_from = c(att_pct, pretrend_p, breakdown_Mbar, honestdid_min_M_zero)
    )
  comparison_csv <- file.path(baseline_run$csv_dir, "anticipation_shift_comparison.csv")
  write.csv(comparison_df, comparison_csv, row.names = FALSE)
  cat("Comparison table saved to:", comparison_csv, "\n\n")

  cat("Simple baseline vs shift interpretation:\n")
  baseline_label <- "Standard stacked PPML"
  shift_label <- paste("Shift", anticipation_periods, "quarters")
  for (outcome_name in unique(combined_summary$outcome)) {
    outcome_rows <- combined_summary %>% filter(outcome == outcome_name)
    baseline_row <- outcome_rows %>% filter(specification == baseline_label) %>% slice(1)
    shift_row <- outcome_rows %>% filter(specification == shift_label) %>% slice(1)

    if (nrow(baseline_row) == 1 && nrow(shift_row) == 1) {
      if (!is.na(shift_row$pretrend_p) && !is.na(baseline_row$pretrend_p) && shift_row$pretrend_p > baseline_row$pretrend_p) {
        cat("  - Shift improves pretrend balance for", toupper(outcome_name), "\n")
      }
      if (!is.na(shift_row$att_pct) && !is.na(baseline_row$att_pct) && abs(shift_row$att_pct) > abs(baseline_row$att_pct)) {
        cat("  - Shift increases estimated magnitude for", toupper(outcome_name), "\n")
      }
    }
  }
  cat("\n")
  
  cat("\n\n─────────────────────────────────────────────────────────────\n")
  cat("ANALYSIS COMPLETE - Baseline + Anticipation Shift\n")
  cat("─────────────────────────────────────────────────────────────\n\n")
  cat("Output location:", baseline_run$output_dir, "\n")
  cat("  - Consolidated results:", baseline_run$results_file, "\n")
  cat("  - Summary table:", summary_csv, "\n")
  cat("  - Standard (baseline) graphs:", file.path(baseline_run$plots_dir, "standard"), "\n")
  cat("  - Anticipation shift graphs:", file.path(baseline_run$plots_dir, "anticipation"), "\n")
  cat("  - CSV exports:", baseline_run$csv_dir, "\n")
} else {
  shift_run <- NULL
}
