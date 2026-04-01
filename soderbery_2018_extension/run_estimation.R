#' ===========================================================================
#' run_estimation.R
#'
#' Unified estimation script for Soderbery (2018) heterogeneous elasticities
#' applied to BACI data. Replaces the separate run_regional.R and
#' run_country.R scripts.
#'
#' Set scope = "regional" or scope = "country" below to control behavior:
#'   - regional: Aggregates to 20 Soderbery Table 1 entities, runs serially.
#'   - country:  Keeps ~200 individual countries, runs in parallel.
#'
#' TWO-STEP PROCEDURE:
#'   Step 1 (Discovery): Estimate with paper defaults as structural
#'     parameters. Discovers BACI-specific distributional properties.
#'   Step 2 (Final): Fix structural defaults (sigma_V, gamma_V) at
#'     Step 1 medians. Re-estimate. These are the final results.
#'
#'   Structural defaults are updated exactly once, not iterated.
#'   Starting values are updated freely for optimizer speed.
#'
#' Both steps detect previously completed results and skip re-estimation.
#' The serial path (ncores = 1) supports checkpoint/resume.
#'
#' Sources feen94_het_baci.R for all core functions.
#'
#' CITATION:
#'   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
#'   Tariffs," JIE, 114, 2018, pp. 44-62.
#'
#' Last updated: 2026-03-31
#' ===========================================================================

setwd("C:/Users/User/Project")
source("feen94_het_baci.R")


# ===========================================================================
#  CONFIG
#
#  Set scope to "regional" or "country" to control aggregation and
#  parallelization. All other settings apply to both scopes.
# ===========================================================================

scope <- "country"    # "regional" or "country"

config <- list(
  # --- BACI column names ---
  value    = "v",  quan = "q",  good = "k",
  importer = "j",  exporter = "i",  time = "t",

  # --- Data ---
  filepath   = "BACI_HS92_V202601/",
  minyear    = 1995,
  maxyear    = NULL,       # NULL = use all available years
  agg_level  = "hs4",

  # --- Aggregation ---
  #   Determined by scope above. Do not set manually.
  use_regions       = NULL,   # Set automatically from scope
  custom_region_map = NULL,

  # --- Filtering ---
  min_exporters        = 2,
  min_destinations     = 2,
  min_periods          = 3,
  uv_outlier_threshold = 2.0,

  # --- Starting values and structural defaults ---
  # Step 1 uses paper Table 2 medians (Comtrade, 1991-2007) as both
  # starting values and structural defaults. Step 2 replaces these
  # with BACI-specific medians from Step 1.
  sigma_start     = 2.88,
  gamma_start     = 0.69,
  sigma_V_default = 2.88,
  gamma_V_default = 0.69,

  # --- Post-estimation trimming ---
  tail_trim_pct = 0.005
)


# ===========================================================================
#  SCOPE-DEPENDENT SETTINGS
# ===========================================================================

if (!scope %in% c("regional", "country")) {
  stop("scope must be 'regional' or 'country', got: ", scope)
}

config$use_regions <- (scope == "regional")

# Regional runs serially (only ~20 entities per product — parallelization
# overhead exceeds computation time). Country-level runs in parallel.
if (scope == "regional") {
  ncores <- 1L
} else {
  #ncores <- max(1L, parallel::detectCores() - 2L)
  ncores <- 4L
}

cat(sprintf("Scope: %s | Cores: %d\n", scope, ncores))


# ===========================================================================
#  OUTPUT PREFIX
# ===========================================================================

out_prefix <- build_output_prefix(config, scope = scope)
cat(sprintf("Output prefix: %s\n\n", out_prefix))


# ===========================================================================
#  LOAD REGIONAL STARTING VALUES (country-level only)
#
#  If regional results exist, they provide per-cell optimizer starting
#  values via country-to-region lookup. This is a speed optimization —
#  it affects where the optimizer begins searching, not the structural
#  parameters in the moment equations.
# ===========================================================================

if (scope == "country") {
  regional_prefix <- build_output_prefix(config, scope = "regional")
  regional_file <- NULL
  for (step in c("step2", "step1")) {
    candidate <- paste0(regional_prefix, "_", step, ".rds")
    if (file.exists(candidate)) { regional_file <- candidate; break }
  }
  # Fallback: legacy naming
  if (is.null(regional_file)) {
    for (p in c("step2", "step1", sprintf("p%d", 8:1))) {
      candidate <- sprintf("elast_regional_%s.rds", p)
      if (file.exists(candidate)) { regional_file <- candidate; break }
    }
  }

  if (!is.null(regional_file)) {
    cat(sprintf("Loading regional results from: %s\n", regional_file))
    regional_results <- readRDS(regional_file)
    config <- init_from_regional(config, regional_results)
  } else {
    cat("No regional results found. Using global defaults for starting values.\n")
    cat("  (Run with scope='regional' first for faster country-level convergence.)\n\n")
  }
}


# ===========================================================================
#  STEP 1: DISCOVERY
# ===========================================================================

step1_file <- paste0(out_prefix, "_step1.rds")
if (file.exists(step1_file)) {
  cat("\n========== STEP 1: LOADING PREVIOUS RESULTS ==========\n")
  cat(sprintf("  Found existing Step 1 results: %s\n", step1_file))
  cat("  Skipping Step 1 estimation. Delete this file to re-run.\n\n")
  results_step1 <- readRDS(step1_file)
} else {
  cat("\n========== STEP 1: DISCOVERY ==========\n")
  cat("  Structural defaults: sigma_V=2.88, gamma_V=0.69 (paper Table 2)\n")
  cat(sprintf("  Purpose: discover BACI-specific distributional properties (%s).\n\n", scope))

  results_step1 <- estimate_all_parallel(config, ncores = ncores)
  saveRDS(results_step1, step1_file)
}

med_sigma_1 <- median(results_step1$sigma, na.rm = TRUE)
med_gamma_1 <- median(results_step1$gamma, na.rm = TRUE)


# ===========================================================================
#  STEP 2: FINAL ESTIMATION
# ===========================================================================

step2_file <- paste0(out_prefix, "_step2.rds")
if (file.exists(step2_file)) {
  cat("\n========== STEP 2: LOADING PREVIOUS RESULTS ==========\n")
  cat(sprintf("  Found existing Step 2 results: %s\n", step2_file))
  cat("  Skipping Step 2 estimation. Delete this file to re-run.\n\n")
  results_step2 <- readRDS(step2_file)
} else {
  cat("\n========== STEP 2: FINAL ESTIMATION ==========\n")
  cat(sprintf("  Structural defaults updated from Step 1: sigma_V=%.3f, gamma_V=%.3f\n",
              med_sigma_1, med_gamma_1))
  cat("  These are FIXED for this pass (no further iteration).\n")
  cat("  Starting values also updated for faster optimizer convergence.\n\n")

  config$sigma_V_default <- med_sigma_1
  config$gamma_V_default <- med_gamma_1
  config$sigma_start     <- med_sigma_1
  config$gamma_start     <- med_gamma_1

  results_step2 <- estimate_all_parallel(config, ncores = ncores)
  saveRDS(results_step2, step2_file)
}

med_sigma_2 <- median(results_step2$sigma, na.rm = TRUE)
med_gamma_2 <- median(results_step2$gamma, na.rm = TRUE)


# ===========================================================================
#  SENSITIVITY COMPARISON
# ===========================================================================

cat("\n========== SENSITIVITY: STEP 1 vs STEP 2 ==========\n\n")
cat(sprintf("  Step 1 medians (Comtrade defaults): sigma=%.3f, gamma=%.3f\n",
            med_sigma_1, med_gamma_1))
cat(sprintf("  Step 2 medians (BACI defaults):     sigma=%.3f, gamma=%.3f\n",
            med_sigma_2, med_gamma_2))
cat(sprintf("  Shift: sigma %+.3f (%.1f%%), gamma %+.3f (%.1f%%)\n",
            med_sigma_2 - med_sigma_1,
            100 * (med_sigma_2 - med_sigma_1) / med_sigma_1,
            med_gamma_2 - med_gamma_1,
            100 * (med_gamma_2 - med_gamma_1) / med_gamma_1))

merged <- merge(
  results_step1[, .(importer, exporter, good, sigma_1 = sigma, gamma_1 = gamma)],
  results_step2[, .(importer, exporter, good, sigma_2 = sigma, gamma_2 = gamma)],
  by = c("importer", "exporter", "good"))

if (nrow(merged) > 10L) {
  cor_sigma <- cor(merged$sigma_1, merged$sigma_2, use = "complete.obs")
  cor_gamma <- cor(merged$gamma_1, merged$gamma_2, use = "complete.obs")
  cat(sprintf("  Cross-step correlation: sigma=%.4f, gamma=%.4f\n", cor_sigma, cor_gamma))
  cat(sprintf("  Matched cells: %s\n", format(nrow(merged), big.mark = ",")))
  cat("\n  Note: Low correlations indicate sensitivity to structural defaults.\n")
  cat("  The Step 2 results use BACI-appropriate defaults and are preferred.\n")
}


# ===========================================================================
#  FINAL OUTPUT
# ===========================================================================

cat(sprintf("\n  FINAL medians: sigma=%.3f, gamma=%.3f\n", med_sigma_2, med_gamma_2))
cat(sprintf("\nFinal results: %s_step2.rds\n", out_prefix))
fwrite(results_step2, paste0(out_prefix, "_step2.csv"))
cat(sprintf("CSV export: %s_step2.csv\n", out_prefix))

# --- Summary report ---
write_estimation_summary(results_step2, config, out_prefix,
                         step1_results = results_step1, scope = scope)
