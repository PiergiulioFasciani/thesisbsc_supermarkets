# ================= run_all_proximity_models_single_txt.R =================
options(stringsAsFactors = FALSE, scipen = 999, digits = 5)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# --- 0) Packages ---
need <- c("lmtest","sandwich")
to_install <- need[!sapply(need, requireNamespace, quietly=TRUE)]
if (length(to_install)) install.packages(to_install)
library(lmtest); library(sandwich)

# --- 1) Load CSV (from CLI arg or picker) ---
args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) > 0 && nzchar(args[1])) args[1] else ""

if (!nzchar(csv_path)) {
  message("No CSV path passed via command line, opening interactive picker...")
  csv_path <- tryCatch({
    if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable())
      rstudioapi::selectFile("Select quadtree CSV", label="Select")
    else if (.Platform$OS.type=="windows")
      utils::choose.files("Select quadtree CSV", multi=FALSE)
    else file.choose()
  }, error=function(e) "")
}

if (!nzchar(csv_path) || !file.exists(csv_path)) {
  stop("Fatal: CSV file not found or specified path is empty. Path: '", csv_path, "'")
}
message("Loading CSV: ", csv_path)
df <- read.csv(csv_path, check.names = TRUE)

# --- 2) Output directory + one summary file ---
out_dir <- file.path(dirname(csv_path),
                     paste0("model_summaries_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "proximity_models_all.txt")
con <- file(out_file, open = "wt", encoding = "UTF-8")
on.exit(close(con), add = TRUE)

# --- Helpers to write to CLI + single file ---
say <- function(...) { msg <- paste(..., collapse=""); cat(msg, "\n"); writeLines(msg, con) }
dump_obj <- function(obj, header=NULL) {
  if (!is.null(header)) { say(""); say(header) }
  cap <- capture.output(print(obj))
  cat(paste(cap, collapse="\n"), "\n")
  writeLines(cap, con)
}
qfun <- function(x) as.numeric(quantile(x[is.finite(x)], c(.25,.50,.75), na.rm=TRUE))

# --- 3) Small utilities & distance sources ---
has <- function(...) any(c(...) %in% names(df))
v <- function(x) if (x %in% names(df)) df[[x]] else rep(NA_real_, nrow(df))
coalesce_vec <- function(...) { xs <- list(...); out <- xs[[1]]
if (length(xs) > 1) for (i in 2:length(xs)) out <- ifelse(is.finite(out), out, xs[[i]]); out }

# Polygon distances with “outside” fallback WHEN original==0
poly_super <- coalesce_vec(
  ifelse(is.finite(v("dist_poly_to_supermarket_m")) & v("dist_poly_to_supermarket_m")==0 &
           is.finite(v("dist_poly_to_supermarket_outside_m")),
         v("dist_poly_to_supermarket_outside_m"),
         v("dist_poly_to_supermarket_m"))
)
poly_pt    <- coalesce_vec(
  ifelse(is.finite(v("dist_poly_to_pt_m")) & v("dist_poly_to_pt_m")==0 &
           is.finite(v("dist_poly_to_pt_outside_m")),
         v("dist_poly_to_pt_outside_m"),
         v("dist_poly_to_pt_m"))
)
poly_duomo <- v("dist_poly_to_duomo_m")

# Row-wise: prefer POI→facility, else polygon (above)
super_raw <- if (has("poi2super_min_m")) coalesce_vec(v("poi2super_min_m"), poly_super) else poly_super
pt_raw    <- if (has("poi2pt_min_m"))    coalesce_vec(v("poi2pt_min_m"),    poly_pt)    else poly_pt
duomo_raw <- if (has("poi2duomo_min_m")) coalesce_vec(v("poi2duomo_min_m"), poly_duomo) else poly_duomo

# Controls if present
cnt_300  <- if (has("super_cnt_300m")) v("super_cnt_300m") else NULL
cnt_600  <- if (has("super_cnt_600m")) v("super_cnt_600m") else NULL
log_area <- if (has("tile_area_m2"))   log(pmax(v("tile_area_m2"), 1))    else NULL

# Must-haves
stopifnot(all(c("has_poi","center_lon","center_lat") %in% names(df)))
df$has_poi <- as.integer(df$has_poi)

# --- 4) Features for LOG specs (uniform logs) ---
df$lsuper <- log1p(pmax(super_raw, 0))
df$lpt    <- log1p(pmax(pt_raw,    0))
df$lduomo <- log1p(pmax(duomo_raw, 0))
if (!is.null(cnt_300)) df$lcnt0_300   <- log1p(pmax(cnt_300,  0))
if (!is.null(cnt_600)) df$lcnt300_600 <- log1p(pmax(cnt_600,  0))
if (!is.null(log_area)) df$log_area   <- log_area

# --- 5) Build datasets for three specs ---
# FULL_LOG needs lpt + lsuper + lduomo (+ optional controls)
keep_full <- Reduce(`&`, lapply(
  c("has_poi","lsuper","lpt","lduomo",
    intersect(c("lcnt0_300","lcnt300_600","log_area"), names(df))),
  function(nm) is.finite(df[[nm]])
))
df_full <- df[keep_full, , drop=FALSE]

# MIN_LOG needs only lsuper + lduomo
keep_min <- Reduce(`&`, lapply(c("has_poi","lsuper","lduomo"),
                               function(nm) is.finite(df[[nm]])))
df_min <- df[keep_min, , drop=FALSE]

# SIMPLE_POLY_LINEAR uses polygon-only linear distances (per 100m / km)
df_simple <- data.frame(
  has_poi = df$has_poi,
  center_lon = df$center_lon,
  center_lat = df$center_lat,
  dist_pt_100m    = poly_pt    / 100,
  dist_super_100m = poly_super / 100,
  dist_duomo_km   = poly_duomo / 1000
)
keep_simple <- Reduce(`&`, lapply(df_simple, is.finite))
df_simple <- df_simple[keep_simple, , drop=FALSE]

# --- 6) Clusters (~1 km) ---
clust_full <- with(df_full, interaction(
  floor((center_lon - min(center_lon, na.rm=TRUE))/0.012),
  floor((center_lat - min(center_lat, na.rm=TRUE))/0.009), drop=TRUE))
clust_min <- with(df_min, interaction(
  floor((center_lon - min(center_lon, na.rm=TRUE))/0.012),
  floor((center_lat - min(center_lat, na.rm=TRUE))/0.009), drop=TRUE))
clust_simple <- with(df_simple, interaction(
  floor((center_lon - min(center_lon, na.rm=TRUE))/0.012),
  floor((center_lat - min(center_lat, na.rm=TRUE))/0.009), drop=TRUE))

# --- 7) Formulas ---
rhs_full <- "lpt + lsuper*lduomo + I(lsuper^2) + I(lduomo^2)"
if ("lcnt0_300"   %in% names(df_full)) rhs_full <- paste(rhs_full, "+ lcnt0_300")
if ("lcnt300_600" %in% names(df_full)) rhs_full <- paste(rhs_full, "+ lcnt300_600")
if ("log_area"    %in% names(df_full)) rhs_full <- paste(rhs_full, "+ log_area")

rhs_min <- "lsuper*lduomo"  # minimal spec
form_lpm_full   <- as.formula(paste0("has_poi ~ ", rhs_full))
form_logit_full <- as.formula(paste0("has_poi ~ ", rhs_full))
form_lpm_min    <- as.formula(paste0("has_poi ~ ", rhs_min))
form_logit_min  <- as.formula(paste0("has_poi ~ ", rhs_min))

form_simple <- has_poi ~ dist_pt_100m + dist_super_100m*dist_duomo_km +
  I(dist_super_100m^2) + I(dist_duomo_km^2)

# --- 8) Fit models ---
lpm_full     <- lm(form_lpm_full, data=df_full)
logit_full   <- glm(form_logit_full, data=df_full, family=binomial("logit"))
lpm_min      <- lm(form_lpm_min, data=df_min)
logit_min    <- glm(form_logit_min, data=df_min, family=binomial("logit"))
lpm_simple   <- lm(form_simple, data=df_simple)
logit_simple <- glm(form_simple, data=df_simple, family=binomial("logit"))

# --- 9) Robust and cluster-robust tables ---
hc1   <- function(m) vcovHC(m, type="HC1")
cl_vc <- function(m, cl) vcovCL(m, cluster=cl, type="HC1")

# FULL_LOG
lpm_full_hc1   <- coeftest(lpm_full,   vcov = hc1(lpm_full))
lpm_full_cl    <- coeftest(lpm_full,   vcov = cl_vc(lpm_full,   clust_full))
logit_full_hc1 <- coeftest(logit_full, vcov = hc1(logit_full))
logit_full_cl  <- coeftest(logit_full, vcov = cl_vc(logit_full, clust_full))

# MIN_LOG
lpm_min_hc1    <- coeftest(lpm_min,    vcov = hc1(lpm_min))
lpm_min_cl     <- coeftest(lpm_min,    vcov = cl_vc(lpm_min,    clust_min))
logit_min_hc1  <- coeftest(logit_min,  vcov = hc1(logit_min))
logit_min_cl   <- coeftest(logit_min,  vcov = cl_vc(logit_min,  clust_min))

# SIMPLE_POLY_LINEAR
lpm_simple_hc1   <- coeftest(lpm_simple,   vcov = hc1(lpm_simple))
lpm_simple_cl    <- coeftest(lpm_simple,   vcov = cl_vc(lpm_simple,   clust_simple))
logit_simple_hc1 <- coeftest(logit_simple, vcov = hc1(logit_simple))
logit_simple_cl  <- coeftest(logit_simple, vcov = cl_vc(logit_simple, clust_simple))

# --- 10) Predicted probabilities (grids) ---
# FULL_LOG: quartiles of lpt, lsuper, lduomo
grid_full <- expand.grid(
  lpt    = qfun(df_full$lpt),
  lsuper = qfun(df_full$lsuper),
  lduomo = qfun(df_full$lduomo)
)
if ("lcnt0_300"   %in% names(df_full)) grid_full$lcnt0_300   <- mean(df_full$lcnt0_300,   na.rm=TRUE)
if ("lcnt300_600" %in% names(df_full)) grid_full$lcnt300_600 <- mean(df_full$lcnt300_600, na.rm=TRUE)
if ("log_area"    %in% names(df_full)) grid_full$log_area    <- mean(df_full$log_area,    na.rm=TRUE)
grid_full$lpm_hat   <- pmin(1, pmax(0, as.numeric(predict(lpm_full,   newdata=grid_full, type="response"))))
grid_full$logit_hat <- as.numeric(plogis(predict(logit_full, newdata=grid_full, type="link")))
grid_full$pt_m    <- round(exp(grid_full$lpt)    - 1)
grid_full$super_m <- round(exp(grid_full$lsuper) - 1)
grid_full$duomo_m <- round(exp(grid_full$lduomo) - 1)

# MIN_LOG: quartiles of lsuper, lduomo
grid_min <- expand.grid(
  lsuper = qfun(df_min$lsuper),
  lduomo = qfun(df_min$lduomo)
)
grid_min$lpm_hat   <- pmin(1, pmax(0, as.numeric(predict(lpm_min,   newdata=grid_min, type="response"))))
grid_min$logit_hat <- as.numeric(plogis(predict(logit_min, newdata=grid_min, type="link")))
grid_min$super_m <- round(exp(grid_min$lsuper) - 1)
grid_min$duomo_m <- round(exp(grid_min$lduomo) - 1)

# SIMPLE_POLY_LINEAR: quartiles of linear distances
grid_simple <- expand.grid(
  dist_pt_100m    = qfun(df_simple$dist_pt_100m),
  dist_super_100m = qfun(df_simple$dist_super_100m),
  dist_duomo_km   = qfun(df_simple$dist_duomo_km)
)
grid_simple$lpm_hat   <- pmin(1, pmax(0, as.numeric(predict(lpm_simple,   newdata=grid_simple, type="response"))))
grid_simple$logit_hat <- as.numeric(predict(logit_simple, newdata=grid_simple, type="response"))

# --- 11) Write ONE consolidated report ---
say("CSV: ", csv_path)
say("Output: ", out_file)
say(sprintf("Rows kept — FULL_LOG: %d | MIN_LOG: %d | SIMPLE: %d (of %d)",
            nrow(df_full), nrow(df_min), nrow(df_simple), nrow(df)))
say(sprintf("Y mean — FULL_LOG: %.4f | MIN_LOG: %.4f | SIMPLE: %.4f",
            mean(df_full$has_poi), mean(df_min$has_poi), mean(df_simple$has_poi)))
say("Distance sources for LOG specs: row-wise POI→facility preferred; polygon with outside fallback otherwise.")
say("RHS — FULL_LOG: ", rhs_full)
say("RHS — MIN_LOG : ", rhs_min)
say("RHS — SIMPLE_POLY_LINEAR: dist_pt_100m + dist_super_100m + dist_duomo_km + squares + interaction")

dump_obj(lpm_full_hc1,   "\n=== FULL_LOG — LPM (HC1) ===")
dump_obj(lpm_full_cl,    "=== FULL_LOG — LPM (cluster ~1km) ===")
dump_obj(logit_full_hc1, "=== FULL_LOG — Logit (HC1) ===")
dump_obj(logit_full_cl,  "=== FULL_LOG — Logit (cluster ~1km) ===")
dump_obj(grid_full[, c("pt_m","super_m","duomo_m","lpm_hat","logit_hat")],
         "=== FULL_LOG — Predicted probabilities (quartiles) ===")

dump_obj(lpm_min_hc1,    "\n=== MIN_LOG — LPM (HC1) ===")
dump_obj(lpm_min_cl,     "=== MIN_LOG — LPM (cluster ~1km) ===")
dump_obj(logit_min_hc1,  "=== MIN_LOG — Logit (HC1) ===")
dump_obj(logit_min_cl,   "=== MIN_LOG — Logit (cluster ~1km) ===")
dump_obj(grid_min[, c("super_m","duomo_m","lpm_hat","logit_hat")],
         "=== MIN_LOG — Predicted probabilities (quartiles) ===")

dump_obj(lpm_simple_hc1,   "\n=== SIMPLE_POLY_LINEAR — LPM (HC1) ===")
dump_obj(lpm_simple_cl,    "=== SIMPLE_POLY_LINEAR — LPM (cluster ~1km) ===")
dump_obj(logit_simple_hc1, "=== SIMPLE_POLY_LINEAR — Logit (HC1) ===")
dump_obj(logit_simple_cl,  "=== SIMPLE_POLY_LINEAR — Logit (cluster ~1km) ===")
dump_obj(grid_simple[, c("dist_pt_100m","dist_super_100m","dist_duomo_km","lpm_hat","logit_hat")],
         "=== SIMPLE_POLY_LINEAR — Predicted probabilities (quartiles) ===")

say("\nAll done. Single consolidated report at:")
say(out_file)
# ================= END run_all_proximity_models_single_txt.R =================
