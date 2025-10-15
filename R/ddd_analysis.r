#!/usr/bin/env Rscript
# =============================================================================
# Triple-Difference (DDD) Analysis for Supermarket Entry Effects
# =============================================================================
#
# Purpose: Estimate baseline DDD, by-ring DDD, and event-study DDD
#          using prepared area × category × month panel data
#
# Inputs:  Panel CSV from data/processed/ddd_panel_*.csv
# Outputs: Tables and figures to output/ddd_results_TIMESTAMP/
#
# Author:  Piergiulio Fasciani (with the aid of GitHub Copilot)
# Date:    October 2025
# =============================================================================

# ========================== SETUP & PACKAGES ==========================

set.seed(42)  # Reproducibility

# Suppress startup messages
options(warn = -1)
suppressPackageStartupMessages({
  library(data.table)    # Fast data manipulation
  library(fixest)        # High-dimensional fixed effects
  library(ggplot2)       # Plotting
  library(kableExtra)    # Table formatting
  library(broom)         # Model tidying
})
options(warn = 0)

# String concatenation helper
`%+%` <- function(a, b) paste0(a, b)

# Console output helpers
banner <- function(msg) {
  cat("\n" %+% strrep("=", 70) %+% "\n")
  cat(msg %+% "\n")
  cat(strrep("=", 70) %+% "\n")
}

info <- function(msg) cat("[INFO] " %+% msg %+% "\n")
warn <- function(msg) cat("[WARN] " %+% msg %+% "\n")
error <- function(msg) {
  cat("[ERROR] " %+% msg %+% "\n")
  stop(msg, call. = FALSE)
}
success <- function(msg) cat("[OK] " %+% msg %+% "\n")

# ========================== INPUT HANDLING ==========================

banner("DDD ANALYSIS: TRIPLE-DIFFERENCE ESTIMATION")

# Get input file from command line or use most recent
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  panel_file <- args[1]
} else {
  # Find most recent panel file
  panel_files <- list.files("data/processed", pattern = "^ddd_panel_.*\\.csv$", 
                           full.names = TRUE)
  if (length(panel_files) == 0) {
    error("No panel files found in data/processed/")
  }
  panel_file <- panel_files[length(panel_files)]
  info("No input specified; using most recent: " %+% basename(panel_file))
}

if (!file.exists(panel_file)) {
  error("Panel file not found: " %+% panel_file)
}

# Create output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path("output", "ddd_results_" %+% timestamp)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

info("Output directory: " %+% output_dir)

# ========================== DATA LOADING ==========================

banner("STEP 1: LOADING AND VALIDATING PANEL DATA")

info("Reading: " %+% panel_file)
panel <- fread(panel_file)
info("Loaded " %+% nrow(panel) %+% " rows, " %+% ncol(panel) %+% " columns")

# ========================== PRE-FLIGHT VALIDATION ==========================

info("Running pre-flight checks...")

# Required columns
required_cols <- c("area_id", "site_id", "month", "category", "A", "ring", 
                   "ring_min_m", "ring_max_m", "y", "t0", "E", "event_time", "ev_cap")

missing_cols <- setdiff(required_cols, names(panel))
if (length(missing_cols) > 0) {
  error("Missing required columns: " %+% paste(missing_cols, collapse = ", "))
}
success("All required columns present")

# Convert date columns
panel[, month := as.Date(month)]
panel[, t0 := as.Date(t0)]

# Ensure A and E are 0/1 integers
panel[, A := as.integer(A)]
panel[, E := as.integer(E)]

if (!all(panel$A %in% c(0, 1))) {
  error("Column A must be 0 or 1")
}
if (!all(panel$E %in% c(0, 1))) {
  error("Column E must be 0 or 1")
}
success("Column types validated")

# Check uniqueness of (area_id, category, month)
dups <- panel[, .N, by = .(area_id, category, month)][N > 1]
if (nrow(dups) > 0) {
  error("Duplicate rows found in (area_id, category, month):\n" %+% 
        paste(head(dups, 5), collapse = "\n"))
}
success("Panel uniqueness: OK (no duplicates)")

# Check for missing outcomes
na_y <- sum(is.na(panel$y))
if (na_y > 0) {
  error("Found " %+% na_y %+% " missing values in outcome y")
}
success("No missing outcomes")

# Exposure coherence check
panel[, E_expected := as.integer(month >= t0)]
incoherent <- panel[E != E_expected]
if (nrow(incoherent) > 0) {
  warn("Found " %+% nrow(incoherent) %+% " rows where E != 1[month >= t0]")
  print(head(incoherent[, .(area_id, month, t0, E, E_expected)]))
  error("Exposure coherence failed - please fix upstream")
}
panel[, E_expected := NULL]
success("Exposure coherence: OK (E = 1[month >= t0])")

# Ring consistency
valid_rings <- c("in", "mid", "out")
invalid_rings <- unique(panel$ring)[!unique(panel$ring) %in% valid_rings]
if (length(invalid_rings) > 0) {
  error("Invalid ring values: " %+% paste(invalid_rings, collapse = ", "))
}

ring_checks <- list(
  "in" = c(0, 500),
  "mid" = c(500, 1000),
  "out" = c(1000, 1500)
)

for (r in names(ring_checks)) {
  expected_min <- ring_checks[[r]][1]
  expected_max <- ring_checks[[r]][2]
  
  wrong <- panel[ring == r & (ring_min_m != expected_min | ring_max_m != expected_max)]
  if (nrow(wrong) > 0) {
    error("Ring '" %+% r %+% "' has inconsistent numeric cutoffs")
  }
}
success("Ring definitions: OK")

# Composition check
a_dist <- table(panel$A)
if (length(a_dist) != 2) {
  error("Both A=0 and A=1 must be present in the data")
}

info("Treatment composition:")
info("  A=1 (affected):  " %+% a_dist["1"] %+% " obs (" %+% 
     round(100 * a_dist["1"] / sum(a_dist), 1) %+% "%)")
info("  A=0 (placebo):   " %+% a_dist["0"] %+% " obs (" %+% 
     round(100 * a_dist["0"] / sum(a_dist), 1) %+% "%)")

# ========================== SAMPLE SUMMARY ==========================

banner("STEP 2: SAMPLE SUMMARY")

n_sites <- uniqueN(panel$site_id)
n_areas <- uniqueN(panel$area_id)
n_months <- uniqueN(panel$month)
n_categories <- uniqueN(panel$category)
n_obs <- nrow(panel)

date_range <- range(panel$month)

info("Sample dimensions:")
info("  Sites (clusters):     " %+% n_sites)
info("  Areas:                " %+% n_areas)
info("  Months:               " %+% n_months)
info("  Categories:           " %+% n_categories)
info("  Total observations:   " %+% format(n_obs, big.mark = ","))
info("  Date range:           " %+% date_range[1] %+% " to " %+% date_range[2])

# Ring distribution
ring_dist <- panel[, .(areas = uniqueN(area_id)), by = ring]
info("Ring distribution:")
for (i in 1:nrow(ring_dist)) {
  info("  " %+% ring_dist$ring[i] %+% ": " %+% ring_dist$areas[i] %+% " areas")
}

# Category distribution
cat_dist <- panel[, .(obs = .N, affected = unique(A)), by = category]
setorder(cat_dist, -affected, category)
info("Category distribution:")
for (i in 1:nrow(cat_dist)) {
  label <- ifelse(cat_dist$affected[i] == 1, "affected", "placebo")
  info("  " %+% cat_dist$category[i] %+% " (A=" %+% cat_dist$affected[i] %+% ", " %+% 
       label %+% "): " %+% format(cat_dist$obs[i], big.mark = ","))
}

# Eligibility check (if available)
if ("eligible" %in% names(panel)) {
  elig_sites <- panel[, .(eligible = unique(eligible)), by = site_id][eligible == 1, .N]
  info("Eligible sites (>=12 months pre/post): " %+% elig_sites %+% "/" %+% n_sites)
}

# ========================== MODEL 1: BASELINE DDD ==========================

banner("STEP 3: BASELINE DDD (ALL RINGS POOLED)")

info("Specification: y ~ E:A | area_id^category + category^month + area_id^month")
info("Clustering: site_id")

# Create interaction variable
panel[, EA := E * A]

# Estimate baseline model
tryCatch({
  model_baseline <- feols(
    y ~ EA | area_id^category + category^month + area_id^month,
    data = panel,
    cluster = ~site_id
  )
  
  success("Baseline DDD model estimated")
  
  # Extract coefficient
  coef_baseline <- coef(model_baseline)["EA"]
  se_baseline <- se(model_baseline)["EA"]
  t_baseline <- coef_baseline / se_baseline
  p_baseline <- 2 * pt(-abs(t_baseline), df = model_baseline$nobs - 1)
  
  info("Results:")
  info("  E×A coefficient:  " %+% sprintf("%.4f", coef_baseline))
  info("  Std. error:       " %+% sprintf("%.4f", se_baseline))
  info("  t-statistic:      " %+% sprintf("%.2f", t_baseline))
  info("  p-value:          " %+% sprintf("%.4f", p_baseline))
  info("  Significance:     " %+% 
       ifelse(p_baseline < 0.01, "***", ifelse(p_baseline < 0.05, "**", 
              ifelse(p_baseline < 0.1, "*", ""))))
  
  # Model diagnostics
  info("Model diagnostics:")
  info("  N observations:   " %+% format(model_baseline$nobs, big.mark = ","))
  info("  N clusters:       " %+% model_baseline$ncluster)
  info("  R² within:        " %+% sprintf("%.4f", r2(model_baseline, type = "wr2")))
  
}, error = function(e) {
  error("Baseline model failed: " %+% conditionMessage(e))
})

# ========================== MODEL 2: DDD BY RING ==========================

banner("STEP 4: DDD BY RING (DISTANCE DECAY)")

info("Specification: y ~ E:A:ring | area_id^category + category^month + area_id^month")
info("Reference ring: out")

# Create ring interactions with E and A
panel[, EA_in := E * A * (ring == "in")]
panel[, EA_mid := E * A * (ring == "mid")]
panel[, EA_out := E * A * (ring == "out")]

tryCatch({
  # Estimate by-ring model (with "out" as reference, so include in and mid)
  model_byring <- feols(
    y ~ EA_in + EA_mid | area_id^category + category^month + area_id^month,
    data = panel,
    cluster = ~site_id
  )
  
  success("By-ring DDD model estimated")
  
  # Extract coefficients
  coef_in <- coef(model_byring)["EA_in"]
  coef_mid <- coef(model_byring)["EA_mid"]
  coef_out <- 0  # Reference category
  
  se_in <- se(model_byring)["EA_in"]
  se_mid <- se(model_byring)["EA_mid"]
  
  info("Results (distance decay):")
  info("  In ring (0-500m):      " %+% sprintf("%.4f", coef_in) %+% 
       " (SE: " %+% sprintf("%.4f", se_in) %+% ")")
  info("  Mid ring (500-1000m):  " %+% sprintf("%.4f", coef_mid) %+% 
       " (SE: " %+% sprintf("%.4f", se_mid) %+% ")")
  info("  Out ring (1000-1500m): " %+% sprintf("%.4f", coef_out) %+% " (reference)")
  
  # Check distance decay pattern
  if (coef_in > coef_mid && coef_mid > coef_out) {
    success("Distance decay pattern CONFIRMED: in > mid > out")
  } else if (coef_in > coef_mid || coef_mid > coef_out) {
    info("Distance decay pattern PARTIAL: not monotonic")
  } else {
    warn("Distance decay pattern NOT observed: consider alternative specification")
  }
  
}, error = function(e) {
  error("By-ring model failed: " %+% conditionMessage(e))
})

# ========================== MODEL 3: EVENT-STUDY DDD ==========================

banner("STEP 5: EVENT-STUDY DDD")

info("Specification: y ~ i(ev_cap, A, ref=-1) | area_id^category + category^month + area_id^month")
info("Event time window: -24 to +24 months")
info("Reference period: -1 month before entry")

tryCatch({
  # Estimate event-study model
  model_eventstudy <- feols(
    y ~ i(ev_cap, A, ref = -1) | area_id^category + category^month + area_id^month,
    data = panel,
    cluster = ~site_id
  )
  
  success("Event-study DDD model estimated")
  
  # Extract coefficients for plotting
  es_coefs <- coef(model_eventstudy)
  es_ses <- se(model_eventstudy)
  
  # Get event times from coefficient names
  coef_names <- names(es_coefs)
  event_times <- as.numeric(gsub("ev_cap::(-?\\d+):A", "\\1", coef_names))
  
  # Create data frame for plotting
  es_data <- data.table(
    event_time = event_times,
    coef = es_coefs,
    se = es_ses,
    ci_lower = es_coefs - 1.96 * es_ses,
    ci_upper = es_coefs + 1.96 * es_ses
  )
  
  # Add reference period
  es_data <- rbind(
    data.table(event_time = -1, coef = 0, se = 0, ci_lower = 0, ci_upper = 0),
    es_data
  )
  setorder(es_data, event_time)
  
  info("Event-study coefficients extracted: " %+% nrow(es_data) %+% " periods")
  
  # Pre-trend test (joint test that all pre-event coefficients = 0)
  pre_coefs <- es_data[event_time < 0 & event_time != -1, coef]
  
  if (length(pre_coefs) > 0) {
    # Use Wald test for joint significance - get coefficient names for pre-period
    pre_coef_names <- coef_names[event_times < 0]
    
    if (length(pre_coef_names) > 0) {
      wald_test <- wald(model_eventstudy, keep = pre_coef_names)
      
      info("Pre-trend test (joint test: all pre-event = 0):")
      info("  F-statistic: " %+% sprintf("%.2f", wald_test$stat))
      info("  p-value:     " %+% sprintf("%.4f", wald_test$p))
      
      if (wald_test$p > 0.05) {
        success("Pre-trends: NOT rejected (parallel trends plausible)")
      } else {
        warn("Pre-trends: REJECTED (parallel trends assumption may be violated)")
      }
    }
  }
  
  # ========================== EVENT-STUDY PLOT ==========================
  
  info("Creating event-study plot...")
  
  p_eventstudy <- ggplot(es_data, aes(x = event_time, y = coef)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray40", size = 0.5) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", size = 0.7) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "steelblue") +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2) +
    scale_x_continuous(breaks = seq(-24, 24, by = 6)) +
    labs(
      title = "Event-Study: Effect of Supermarket Entry on Small Retail (DDD)",
      subtitle = paste0(
        "FE: area×category + category×month + area×month | ",
        "Cluster: site_id (N=", model_eventstudy$ncluster, ") | ",
        "Reference: t=-1"
      ),
      x = "Months Since Supermarket Entry",
      y = "Effect on Affected vs. Placebo Categories (DDD)",
      caption = paste0(
        "Note: 95% confidence intervals shown. ",
        "Vertical line marks entry (t=0). ",
        "Affected categories: bakery, butcher, greengrocer, convenience. ",
        "Placebo: funeral_home, hairdresser, restaurant."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray40"),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  plot_file <- file.path(figures_dir, "event_study_ddd.png")
  ggsave(plot_file, p_eventstudy, width = 10, height = 6, dpi = 300)
  success("Event-study plot saved: " %+% plot_file)
  
}, error = function(e) {
  error("Event-study model failed: " %+% conditionMessage(e))
})

# ========================== SAVE TABLES ==========================

banner("STEP 6: SAVING RESULTS TABLES")

# Function to format model output
format_model_table <- function(model, title, file_name) {
  
  # Extract key statistics
  coefs <- coef(model)
  ses <- se(model)
  pvals <- summary(model)$coeftable[, 4]
  
  # Create significance stars
  stars <- sapply(pvals, function(p) {
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    if (p < 0.1) return("*")
    return("")
  })
  
  # Build table
  results <- data.frame(
    Variable = names(coefs),
    Coefficient = sprintf("%.4f", coefs),
    Std_Error = sprintf("%.4f", ses),
    t_statistic = sprintf("%.2f", coefs / ses),
    p_value = sprintf("%.4f", pvals),
    Sig = stars,
    stringsAsFactors = FALSE
  )
  
  # Add model statistics
  footer <- paste0(
    "Observations: ", format(model$nobs, big.mark = ","), "\n",
    "Clusters (sites): ", model$ncluster, "\n",
    "R² within: ", sprintf("%.4f", r2(model, type = "wr2")), "\n",
    "Fixed Effects: area×category, category×month, area×month\n",
    "Standard errors clustered by site_id\n",
    "Significance: * p<0.1, ** p<0.05, *** p<0.01"
  )
  
  # Save as CSV
  write.csv(results, file.path(tables_dir, file_name), row.names = FALSE)
  
  # Also save formatted LaTeX
  kbl(results, format = "latex", booktabs = TRUE, caption = title) %>%
    kable_styling(latex_options = c("striped", "hold_position")) %>%
    add_footnote(footer, notation = "none") %>%
    save_kable(file.path(tables_dir, gsub(".csv", ".tex", file_name)))
  
  success("Saved: " %+% file_name)
}

# Save baseline model
format_model_table(model_baseline, 
                  "Baseline DDD: E×A Effect on Small Retail Stock",
                  "baseline_ddd.csv")

# Save by-ring model
format_model_table(model_byring,
                  "DDD by Ring: Distance Decay Pattern",
                  "byring_ddd.csv")

# Save event-study coefficients
es_table <- es_data[, .(
  event_time = event_time,
  coefficient = sprintf("%.4f", coef),
  std_error = sprintf("%.4f", se),
  ci_lower = sprintf("%.4f", ci_lower),
  ci_upper = sprintf("%.4f", ci_upper)
)]

write.csv(es_table, file.path(tables_dir, "eventstudy_ddd.csv"), row.names = FALSE)
success("Saved: eventstudy_ddd.csv")

# ========================== SUMMARY REPORT ==========================

banner("STEP 7: GENERATING SUMMARY REPORT")

report_lines <- c(
  "# DDD ANALYSIS SUMMARY REPORT",
  "# Generated: " %+% Sys.time(),
  "# Panel file: " %+% basename(panel_file),
  "",
  "## SAMPLE OVERVIEW",
  "Sites (clusters):     " %+% n_sites,
  "Areas:                " %+% n_areas,
  "Months:               " %+% n_months,
  "Categories:           " %+% n_categories,
  "Total observations:   " %+% format(n_obs, big.mark = ","),
  "Date range:           " %+% date_range[1] %+% " to " %+% date_range[2],
  "",
  "Treatment shares:",
  "  A=1 (affected):     " %+% round(100 * a_dist["1"] / sum(a_dist), 1) %+% "%",
  "  A=0 (placebo):      " %+% round(100 * a_dist["0"] / sum(a_dist), 1) %+% "%",
  "",
  "## MODEL SPECIFICATIONS",
  "Fixed Effects: area×category + category×month + area×month",
  "Clustering: site_id",
  "Ring definitions: in=[0,500]m, mid=(500,1000]m, out=(1000,1500]m",
  "Panel window: 2018-01 to 2025-09",
  "Event time: -24 to +24 months (capped)",
  "Reference period: t=-1",
  "",
  "## BASELINE DDD RESULTS",
  "E×A coefficient:  " %+% sprintf("%.4f", coef_baseline) %+% 
    " (SE: " %+% sprintf("%.4f", se_baseline) %+% ")",
  "p-value:          " %+% sprintf("%.4f", p_baseline),
  "Interpretation:   " %+% ifelse(p_baseline < 0.05, "SIGNIFICANT", "NOT significant"),
  "",
  "## BY-RING DDD RESULTS (Distance Decay)",
  "In ring (0-500m):      " %+% sprintf("%.4f", coef_in) %+% 
    " (SE: " %+% sprintf("%.4f", se_in) %+% ")",
  "Mid ring (500-1000m):  " %+% sprintf("%.4f", coef_mid) %+% 
    " (SE: " %+% sprintf("%.4f", se_mid) %+% ")",
  "Out ring (1000-1500m): 0.0000 (reference)",
  "Pattern:               " %+% ifelse(coef_in > coef_mid && coef_mid > 0, 
                                       "Distance decay CONFIRMED", 
                                       "Distance decay pattern unclear"),
  "",
  "## EVENT-STUDY RESULTS",
  "Pre-trend test p-value: " %+% sprintf("%.4f", wald_test$p),
  "Pre-trends status:      " %+% ifelse(wald_test$p > 0.05, 
                                        "NOT rejected (parallel trends OK)", 
                                        "REJECTED (caution advised)"),
  "",
  "## OUTPUT FILES",
  "Tables:",
  "  - baseline_ddd.csv / .tex",
  "  - byring_ddd.csv / .tex",
  "  - eventstudy_ddd.csv",
  "",
  "Figures:",
  "  - event_study_ddd.png",
  "",
  "## INTERPRETATION NOTES",
  "- Baseline E×A: Overall effect of supermarket entry on affected vs placebo",
  "- By-ring: Distance decay suggests localized spillover effects",
  "- Event-study: Dynamic effects over time; pre-trends validate identification",
  "- Clustering at site level accounts for spatial correlation",
  "",
  "================================ END OF REPORT ================================"
)

report_file <- file.path(output_dir, "SUMMARY_REPORT.txt")
writeLines(report_lines, report_file)
success("Summary report saved: " %+% report_file)

# ========================== FINAL MESSAGE ==========================

banner("ANALYSIS COMPLETE")

info("All results saved to: " %+% output_dir)
info("Tables directory:     " %+% tables_dir)
info("Figures directory:    " %+% figures_dir)

cat("\nKey findings:\n")
cat("  Baseline DDD (E×A):    " %+% sprintf("%.4f", coef_baseline) %+% 
    " (p=" %+% sprintf("%.3f", p_baseline) %+% ")\n")
cat("  Distance decay:        " %+% 
    ifelse(coef_in > coef_mid && coef_mid > 0, "Confirmed ✓", "Not confirmed ✗") %+% "\n")
cat("  Pre-trends:            " %+% 
    ifelse(wald_test$p > 0.05, "Valid ✓", "Violated ✗") %+% "\n")

cat("\n" %+% strrep("=", 70) %+% "\n")
success("DDD analysis pipeline completed successfully!")
cat(strrep("=", 70) %+% "\n\n")
