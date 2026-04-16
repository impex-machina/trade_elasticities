#' ===========================================================================
#' run_estimation.R
#'
#' Three-stage estimation pipeline for Soderbery (2018) heterogeneous
#' elasticities, resolving both the sigma-gamma identification ridge
#' and the gamma plateau at large values.
#'
#' STAGE 1: Feenstra (1994) sigma + gamma_common
#'   Estimates sigma and a common gamma per (importer, product) cell
#'   using import-side moments under homogeneity. Well-identified
#'   2-parameter problem. Sigma is used as fixed input for Stages 2a/2b.
#'   Gamma_common provides the shrinkage prior for Stage 2a.
#'
#' STAGE 2a: Regional gamma (fixed sigma, shrinkage toward Feenstra gamma)
#'   Estimates heterogeneous gamma for 20 Soderbery Table 1 regional
#'   entities. Sigma fixed from Stage 1. Shrinkage penalty pulls gamma
#'   toward Stage 1's gamma_common (aggregated to product level) to
#'   prevent escape to the gamma/(1+gamma) -> 1 plateau.
#'   Results provide shrinkage priors for Stage 2b.
#'
#' STAGE 2b: Country gamma (fixed sigma, shrinkage toward Stage 2a gamma)
#'   Estimates heterogeneous gamma for individual countries. Sigma fixed
#'   from Stage 1. Shrinkage toward Stage 2a regional priors.
#'
#' HIERARCHY: Feenstra gamma -> Regional gamma -> Country gamma
#'   Each stage anchors the next. Shrinkage is continuous — dense cells
#'   with strong moments override the prior; thin cells are regularized.
#'
#' OUTPUTS SAVED:
#'   Stage 1:  {prefix}_feenstra_sigma.rds
#'   Stage 2a: {regional_prefix}_fixed_sigma.rds + _summary.rds + _summary.txt
#'   Stage 2b: {country_prefix}_fixed_sigma.rds + _summary.rds + _summary.txt + .csv
#'
#' CITATION:
#'   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
#'   Tariffs," JIE, 114, 2018, pp. 44-62.
#'
#' Last updated: 2026-04-16 (patched: Option 2 bias fix, V-lookups,
#'   checkpointing, Tier 3 separation, perf optimizations, lint)
#' ===========================================================================

#' ===========================================================================
#' Working directory setup
#'
#' The pipeline needs to locate feen94_het_baci.R, het_obj.R, and the two
#' .cpp files. Precedence:
#'   1. Env var SODERBERY_DIR (set via `export SODERBERY_DIR=/path/to/code`)
#'   2. Directory of this script, if running via Rscript
#'   3. Current working directory (no-op)
#'
#' Avoids hardcoding a path that needs sed-ing for every new host.
#' ===========================================================================

.set_working_dir <- function() {
  env_dir <- Sys.getenv("SODERBERY_DIR", unset = "")
  if (nzchar(env_dir) && dir.exists(env_dir)) {
    setwd(env_dir)
    cat(sprintf("Working dir (SODERBERY_DIR): %s\n", env_dir))
    return(invisible(env_dir))
  }

  # Rscript-invoked: locate this script's directory
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
    if (dir.exists(script_dir)) {
      setwd(script_dir)
      cat(sprintf("Working dir (script location): %s\n", script_dir))
      return(invisible(script_dir))
    }
  }

  # Otherwise leave cwd alone — assume user ran source() from the right place
  cat(sprintf("Working dir (unchanged): %s\n", getwd()))
  invisible(getwd())
}

.set_working_dir()
source("feen94_het_baci.R")


# ===========================================================================
#  CONFIG
# ===========================================================================

config <- list(
  # --- BACI column names ---
  value = "v", quan = "q", good = "k",
  importer = "j", exporter = "i", time = "t",

  # --- Data ---
  filepath   = "BACI_HS92_V202601/",
  minyear    = 1995,
  maxyear    = NULL,
  agg_level  = "hs4",

  # --- Aggregation (set per stage) ---
  use_regions       = NULL,
  custom_region_map = NULL,

  # --- Filtering ---
  min_exporters        = 2,
  min_destinations     = 2,
  min_periods          = 3,
  uv_outlier_threshold = 2.0,

  # --- Starting values (Soderbery Table 2 medians) ---
  sigma_start     = 2.88,
  gamma_start     = 0.69,
  sigma_V_default = 2.88,
  gamma_V_default = 0.69,

  # --- Shrinkage ---
  #   Applied at BOTH Stage 2a and Stage 2b.
  #   Stage 2a: anchored to Feenstra gamma_common (prevents gamma plateau)
  #   Stage 2b: anchored to Stage 2a regional gamma
  #   The penalty is: lambda * sum_j (ln(gamma_j) - ln_gamma_prior)^2
  #   Higher lambda = stronger pull toward prior.
  #   Suggested range: 0.1 - 1.0.
  shrinkage_lambda = 0.1,

  # --- Across-exporter weighting ---
  #   "uniform":     rep(1, J). Equal weight per exporter.
  #   "trade_value": weighted by cell trade value, normalized to sum J.
  #                  Large bilateral relationships contribute more to the
  #                  objective. Motivated by Soderbery (2018) fn 14 and
  #                  produces more dispersion in within-pair gamma
  #                  (closer to Soderbery's reported MAD).
  exporter_weight = "trade_value",

  # --- Weight period floor ---
  #   Under trade_value weighting, exporters with fewer than this many
  #   periods have their weight scaled by n_periods / weight_period_floor.
  #   Protects against lumpy single-year shipments (aircraft, LNG cargoes,
  #   rough diamonds) dominating the objective via raw trade value despite
  #   providing little identifying information. Set to 0 to disable.
  #   Default 10 is a Soderbery-sample-friendly threshold; at HS6 with
  #   more intermittent bilateral relationships, might lower to 5-7.
  weight_period_floor = 10L,

  # --- Tier thresholds ---
  #   Set to match Soderbery's minimum data requirements so all
  #   qualifying exporters get full import + export side estimation.
  #   The shrinkage penalty handles regularization of thin cells
  #   continuously — no discrete cutoffs needed.
  tier1_min_periods = 3L,
  tier1_min_dests   = 2L,
  tier2_min_periods = 3L,

  # --- Post-estimation trimming ---
  tail_trim_pct = 0.005
)

ncores <- max(1L, parallel::detectCores() - 2L)
out_base_country  <- build_output_prefix(config, scope = "country")
out_base_regional <- build_output_prefix(config, scope = "regional")

cat(sprintf("Three-stage pipeline | Cores: %d\n", ncores))
cat(sprintf("Country output base:  %s\n", out_base_country))
cat(sprintf("Regional output base: %s\n\n", out_base_regional))


# ===========================================================================
#  DATA PREPARATION — LOAD ONCE, REUSE EVERYWHERE
# ===========================================================================

raw_cache_file <- paste0(out_base_country, "_raw_cache.rds")

if (file.exists(raw_cache_file)) {
  cat("Loading cached raw data...\n")
  raw_cache <- readRDS(raw_cache_file)
  cat(sprintf("  %s obs from cache\n\n", format(nrow(raw_cache), big.mark = ",")))
} else {
  config_raw <- config
  config_raw$use_regions <- FALSE
  raw_cache <- prepare_raw_data(config_raw)
  saveRDS(raw_cache, raw_cache_file)
  cat(sprintf("  Cached to: %s\n\n", raw_cache_file))
}

# Finalize for country-level (Stages 1 and 2b)
cat("Preparing country-level estimation data...\n")
config_country <- config
config_country$use_regions <- FALSE
prep_country <- prepare_data(config_country, raw_cache = raw_cache)
dt_country <- prep_country$dt
cat(sprintf("  Country data: %s obs\n\n",
            format(nrow(dt_country), big.mark = ",")))

# Finalize for regional level (Stage 2a)
cat("Preparing regional estimation data...\n")
config_regional <- config
config_regional$use_regions <- TRUE
prep_regional <- prepare_data(config_regional, raw_cache = raw_cache)
dt_regional <- prep_regional$dt
cat(sprintf("  Regional data: %s obs\n\n",
            format(nrow(dt_regional), big.mark = ",")))


# ===========================================================================
#  STAGE 1: FEENSTRA (1994) SIGMA + GAMMA_COMMON
# ===========================================================================

sigma_file <- paste0(out_base_country, "_feenstra_sigma.rds")

if (file.exists(sigma_file)) {
  cat("\n========== STAGE 1: LOADING SIGMA ESTIMATES ==========\n")
  cat(sprintf("  Found: %s\n\n", sigma_file))
  sigma_estimates <- readRDS(sigma_file)
} else {
  cat("\n========== STAGE 1: FEENSTRA (1994) SIGMA ==========\n")

  config_s1 <- config_country
  sigma_estimates <- estimate_all_feenstra_sigma(
    config_s1, ncores = ncores, prepared_dt = dt_country)
  saveRDS(sigma_estimates, sigma_file)
}

# Clean: drop NAs, non-converged, implausible sigma (<=1).
# The post-estimation tail_trim_pct already handles extreme upper values
# adaptively; no need for a magic-number hard cap like sigma < 200 that
# would drop genuinely-identified high-substitution products.
sigma_clean <- sigma_estimates[!is.na(sigma) & sigma > 1 & convergence == 0]
sigma_fallback <- median(sigma_clean$sigma, na.rm = TRUE)

cat(sprintf("\nStage 1: %s clean cells (of %s)\n",
            format(nrow(sigma_clean), big.mark = ","),
            format(nrow(sigma_estimates), big.mark = ",")))
cat(sprintf("  sigma median=%.3f, IQR=[%.3f, %.3f]\n",
            median(sigma_clean$sigma),
            quantile(sigma_clean$sigma, 0.25),
            quantile(sigma_clean$sigma, 0.75)))
cat(sprintf("  gamma_common median=%.3f (used as Stage 2a prior)\n\n",
            median(sigma_clean$gamma, na.rm = TRUE)))


# ===========================================================================
#  BUILD STAGE 2a STARTING VALUE / PLATEAU-FALLBACK FROM FEENSTRA GAMMA
#
#  METHODOLOGICAL DEPARTURE FROM PREVIOUS DESIGN:
#
#  Previously, Feenstra gamma_common was used as the Stage 2a log-space
#  shrinkage target. But Feenstra sigma is upward-biased under
#  heterogeneity (documented in README), and the bounded interaction
#  gamma/(1+gamma) implies gamma_common is correspondingly downward-
#  biased. Using that biased value as a shrinkage target in Stage 2a
#  propagated the bias to Stage 2a estimates, and thence to Stage 2b
#  via the Stage 2a prior.
#
#  New design (Option 2):
#    - Stage 2a regional is well-identified per the README rationale.
#      Run it WITHOUT shrinkage (lambda = 0). The 20-entity panel
#      handles the gamma plateau naturally via aggregation.
#    - Post-hoc, flag any Stage 2a cells where gamma > plateau_threshold
#      (default 20; gamma/(1+gamma) > 0.95). For those cells, replace
#      the estimate with a Feenstra-anchored fallback. This uses the
#      biased prior only where the unshrunk estimate was also unreliable.
#    - Stage 2b uses the (cleaned) Stage 2a output as its prior. Stage 2b
#      retains shrinkage at country level where identification is weaker.
#
#  The Feenstra gamma_common is still computed here to serve as:
#    (a) default starting value for the optimizer
#    (b) fallback replacement value for plateau-flagged regional cells
# ===========================================================================

feenstra_gamma_clean <- sigma_clean[!is.na(gamma) & gamma > 0]
feenstra_priors <- feenstra_gamma_clean[, .(
  ln_gamma_prior = median(log(gamma), na.rm = TRUE)
), by = good]

cat(sprintf("Stage 2a starting values / plateau fallback from Feenstra: %d products, median gamma=%.3f\n\n",
            nrow(feenstra_priors),
            exp(median(feenstra_priors$ln_gamma_prior))))


# ===========================================================================
#  STAGE 2a: REGIONAL GAMMA (fixed sigma, NO shrinkage — Option 2)
# ===========================================================================

regional_file <- paste0(out_base_regional, "_fixed_sigma.rds")

if (file.exists(regional_file)) {
  cat("========== STAGE 2a: LOADING REGIONAL ESTIMATES ==========\n")
  cat(sprintf("  Found: %s\n\n", regional_file))
  regional_results <- readRDS(regional_file)
} else {
  cat("========== STAGE 2a: REGIONAL GAMMA (fixed sigma, unshrunk) ==========\n")
  cat("  Shrinkage lambda=0 (regional panel well-identified)\n")
  cat("  Plateau fallback: gamma > 20 replaced by Feenstra anchor\n\n")

  config_2a <- config_regional
  # Option 2: no shrinkage at regional level. Stage 2a estimates are
  # used as-is, subject only to plateau fallback.
  config_2a$shrinkage_lambda <- 0.05
  config_2a$shrinkage_priors <- feenstra_priors

  # Sigma lookup: aggregate country-level Feenstra sigma to regions
  rmap <- build_region_map()
  sigma_wr <- copy(sigma_clean)
  sigma_wr[, region := assign_regions(as.integer(importer), rmap)]
  regional_sigma <- sigma_wr[, .(sigma = median(sigma, na.rm = TRUE)),
                              by = .(importer = region, good)]
  config_2a$sigma_lookup   <- regional_sigma
  config_2a$sigma_fallback <- sigma_fallback

  # Starting values from Feenstra medians (used by optimizer only)
  config_2a$sigma_V_default <- sigma_fallback
  config_2a$gamma_V_default <- median(sigma_clean$gamma, na.rm = TRUE)
  config_2a$sigma_start     <- sigma_fallback
  config_2a$gamma_start     <- config_2a$gamma_V_default

  # --- Reference-destination sigma lookup (Stage 2a) ---
  # regional_sigma already exists above: (region, good) -> sigma, aggregated
  # from country-level Stage 1 sigma via median. This is the correct
  # object for sig_V lookup inside build_export_moments when ref_dest
  # resolves to a region.
  config_2a$sigma_V_lookup <- regional_sigma
  # No gamma_V lookup at Stage 2a: Stage 2a is where regional gamma gets
  # estimated in the first place. Fall back to global median (already
  # in cfg$gamma_V_default above). This is the same behavior as the
  # previous code.
  config_2a$gamma_V_lookup <- NULL

  regional_results <- estimate_all_fixed_sigma(
    config_2a, ncores = ncores, prepared_dt = dt_regional)

  # --- Plateau fallback ---
  # gamma/(1+gamma) saturates at 1 as gamma grows. Above gamma=20,
  # gamma/(1+gamma) > 0.95 — the estimate is functionally on the
  # plateau and carries little identifying information. Replace
  # these with the Feenstra anchor for that product.
  plateau_threshold <- 20
  # If tier column exists, exclude Tier 3 (already assigned from prior)
  # from the plateau-fallback mask.
  has_tier <- "tier" %in% names(regional_results)
  plateau_idx <- if (has_tier) {
    regional_results$gamma > plateau_threshold &
      regional_results$tier != 3L &
      !is.na(regional_results$gamma)
  } else {
    regional_results$gamma > plateau_threshold &
      !is.na(regional_results$gamma)
  }
  plateau_idx[is.na(plateau_idx)] <- FALSE
  n_plateau <- sum(plateau_idx)

  if (n_plateau > 0L) {
    cat(sprintf("  Plateau fallback: %d regional cells with gamma > %g replaced\n",
                n_plateau, plateau_threshold))
    # Merge in the Feenstra-anchored fallback value per product
    regional_results <- merge(
      regional_results, feenstra_priors, by = "good", all.x = TRUE)
    # Recompute mask after merge (row order preserved but defensive)
    plateau_idx2 <- if (has_tier) {
      regional_results$gamma > plateau_threshold &
        regional_results$tier != 3L &
        !is.na(regional_results$ln_gamma_prior)
    } else {
      regional_results$gamma > plateau_threshold &
        !is.na(regional_results$ln_gamma_prior)
    }
    plateau_idx2[is.na(plateau_idx2)] <- FALSE
    regional_results[plateau_idx2, `:=`(
      gamma       = exp(ln_gamma_prior),
      convergence = -2L  # new flag: plateau-fallback (distinct from -1 = Tier 3)
    )]
    regional_results[, ln_gamma_prior := NULL]

    # --- Recompute opt_tariff for affected (importer, good) cells ---
    # Changing gamma in any exporter row invalidates the cell-level
    # optimal tariff aggregate. Identify cells that had any plateau
    # replacement and recompute from the updated (sigma, gamma, avg_trade).
    affected_cells <- unique(regional_results[plateau_idx2,
                                               .(importer, good)])
    if (nrow(affected_cells) > 0L) {
      # For each affected cell, recompute using tier < 3 rows only
      # (matching the primary opt_tariff definition). avg_trade was
      # attached at estimate_product_fixed_sigma time; if not present,
      # fall back to uniform weights.
      has_avg_trade <- "avg_trade" %in% names(regional_results)
      for (i in seq_len(nrow(affected_cells))) {
        imp_i  <- affected_cells$importer[i]
        good_i <- affected_cells$good[i]
        cell_rows <- regional_results[importer == imp_i & good == good_i]
        if (has_tier) {
          est_mask <- !is.na(cell_rows$tier) & cell_rows$tier < 3L
        } else {
          est_mask <- rep(TRUE, nrow(cell_rows))
        }
        if (!any(est_mask)) next
        tv <- if (has_avg_trade) cell_rows$avg_trade[est_mask] else NULL
        ot_new <- optimal_tariff(cell_rows$gamma[est_mask],
                                  cell_rows$sigma[est_mask][1], tv)
        regional_results[importer == imp_i & good == good_i,
                          opt_tariff := ot_new]
        # opt_tariff_all recomputed with all rows if that column exists
        if ("opt_tariff_all" %in% names(regional_results)) {
          tv_all <- if (has_avg_trade) cell_rows$avg_trade else NULL
          ot_all_new <- optimal_tariff(cell_rows$gamma,
                                        cell_rows$sigma[1], tv_all)
          regional_results[importer == imp_i & good == good_i,
                            opt_tariff_all := ot_all_new]
        }
      }
    }
  } else {
    cat("  Plateau fallback: 0 regional cells required replacement\n")
  }

  saveRDS(regional_results, regional_file)
}

cat(sprintf("Stage 2a: %s estimates, sigma=%.3f, gamma=%.3f, opt_tariff=%.3f\n\n",
            format(nrow(regional_results), big.mark = ","),
            median(regional_results$sigma, na.rm = TRUE),
            median(regional_results$gamma, na.rm = TRUE),
            median(regional_results$opt_tariff, na.rm = TRUE)))

# --- Save regional summary ---
cat("Saving regional summary...\n")
write_estimation_summary(regional_results, config_regional, out_base_regional,
                         step1_results = NULL, scope = "regional")
fwrite(regional_results, paste0(out_base_regional, "_fixed_sigma.csv"))
cat(sprintf("  CSV: %s_fixed_sigma.csv\n\n", out_base_regional))


# ===========================================================================
#  BUILD STAGE 2b SHRINKAGE PRIORS FROM STAGE 2a
#
#  For each product, compute median ln(gamma) across regional exporters.
#  This is the shrinkage target for country-level estimation.
# ===========================================================================

regional_clean <- regional_results[!is.na(gamma) & gamma > 0]
country_priors <- regional_clean[, .(
  ln_gamma_prior = median(log(gamma), na.rm = TRUE)
), by = good]

cat(sprintf("Stage 2b priors (from Stage 2a): %d products, median gamma=%.3f\n\n",
            nrow(country_priors),
            exp(median(country_priors$ln_gamma_prior))))


# ===========================================================================
#  STAGE 2b: COUNTRY GAMMA (fixed sigma + shrinkage toward Stage 2a)
# ===========================================================================

country_file <- paste0(out_base_country, "_fixed_sigma.rds")

if (file.exists(country_file)) {
  cat("========== STAGE 2b: LOADING COUNTRY ESTIMATES ==========\n")
  cat(sprintf("  Found: %s\n\n", country_file))
  country_results <- readRDS(country_file)
} else {
  cat("========== STAGE 2b: COUNTRY GAMMA (fixed sigma + shrinkage) ==========\n")
  cat(sprintf("  Shrinkage lambda=%.3f, prior=Stage 2a regional gamma\n\n",
              config$shrinkage_lambda))

  config_2b <- config_country
  config_2b$sigma_lookup     <- sigma_clean[, .(importer, good, sigma)]
  config_2b$sigma_fallback   <- sigma_fallback
  config_2b$shrinkage_priors <- country_priors
  config_2b$shrinkage_lambda <- config$shrinkage_lambda
  config_2b$tier1_min_periods <- config$tier1_min_periods
  config_2b$tier1_min_dests   <- config$tier1_min_dests
  config_2b$tier2_min_periods <- config$tier2_min_periods

  config_2b$sigma_V_default <- sigma_fallback
  config_2b$gamma_V_default <- exp(median(country_priors$ln_gamma_prior))
  config_2b$sigma_start     <- sigma_fallback
  config_2b$gamma_start     <- config_2b$gamma_V_default

  # --- Reference-destination elasticity lookups for Stage 2b ---
  #
  # For each Tier 1 exporter, the inner loop of build_export_moments
  # picks a ref_dest (a country at Stage 2b) and needs (sigma_V, gamma_V)
  # for that market. Structural interpretation:
  #   sigma_V ~ elasticity of substitution in ref_dest's import market,
  #             which is a (importer, good) object. Stage 1 provides this
  #             at country level directly.
  #   gamma_V ~ inverse export supply elasticity faced when shipping to
  #             ref_dest. A (importer, exporter, good) object; we don't
  #             have bilateral gamma at ref_dest level yet (that's what
  #             Stage 2b is estimating). Best available proxy: aggregate
  #             Stage 2a regional gamma to (country, good) by mapping each
  #             country to its region via build_region_map(). Coarse but
  #             structurally consistent with the regional priors already
  #             in use for shrinkage.
  config_2b$sigma_V_lookup <- sigma_clean[, .(importer, good, sigma)]

  # Map Stage 2a regional gamma to (country, good): each country in the
  # BACI sample inherits its region's median gamma for that product.
  rmap_2b <- build_region_map()
  gam_V_regional <- regional_clean[, .(gamma = median(gamma, na.rm = TRUE)),
                                     by = .(region = importer, good)]
  country_codes_2b <- unique(as.integer(dt_country$importer))
  cty_to_region <- data.table(
    cty_code = country_codes_2b,
    region   = assign_regions(country_codes_2b, rmap_2b))
  gam_V_country <- merge(cty_to_region, gam_V_regional, by = "region",
                          allow.cartesian = TRUE)
  config_2b$gamma_V_lookup <- gam_V_country[, .(
    importer = as.character(cty_code), good, gamma)]

  cat(sprintf("  Stage 2b sig_V lookup: %s (importer, good) entries\n",
              format(nrow(config_2b$sigma_V_lookup), big.mark = ",")))
  cat(sprintf("  Stage 2b gam_V lookup: %s (country, good) entries from regional\n",
              format(nrow(config_2b$gamma_V_lookup), big.mark = ",")))

  # Regional starting values for optimizer speed
  config_2b <- init_from_regional(config_2b, regional_results)

  country_results <- estimate_all_fixed_sigma(
    config_2b, ncores = ncores, prepared_dt = dt_country)
  saveRDS(country_results, country_file)
}


# ===========================================================================
#  SAVE COUNTRY SUMMARY + CSV
# ===========================================================================

cat("\nSaving country summary...\n")
write_estimation_summary(country_results, config, out_base_country,
                         step1_results = NULL, scope = "country")
fwrite(country_results, paste0(out_base_country, "_fixed_sigma.csv"))
cat(sprintf("  CSV: %s_fixed_sigma.csv\n", out_base_country))


# ===========================================================================
#  FINAL REPORT
# ===========================================================================

cat("\n================================================================\n")
cat("  THREE-STAGE ESTIMATION COMPLETE\n")
cat("================================================================\n\n")

cat(sprintf("  Stage 1 (Feenstra sigma):    %s cells, sigma=%.3f, gamma_common=%.3f\n",
            format(nrow(sigma_clean), big.mark = ","),
            median(sigma_clean$sigma),
            median(sigma_clean$gamma, na.rm=TRUE)))
cat(sprintf("  Stage 2a (Regional gamma):   %s estimates, gamma=%.3f, tariff=%.3f\n",
            format(nrow(regional_results), big.mark = ","),
            median(regional_results$gamma, na.rm=TRUE),
            median(regional_results$opt_tariff, na.rm=TRUE)))
cat(sprintf("  Stage 2b (Country gamma):    %s estimates, gamma=%.3f, tariff=%.3f\n",
            format(nrow(country_results), big.mark = ","),
            median(country_results$gamma, na.rm=TRUE),
            median(country_results$opt_tariff, na.rm=TRUE)))

# --- Structural ratios ---
cr <- country_results[!is.na(sigma) & !is.na(gamma) & gamma > 0]
cat("\n  Structural ratios (country):\n")
cat(sprintf("    gamma/(1+gamma):              %.3f  (Soderbery: 0.408)\n",
            median(cr$gamma / (1 + cr$gamma))))
cat(sprintf("    1/(sigma-1):                  %.3f  (Soderbery: 0.532)\n",
            median(1 / (cr$sigma - 1))))
cat(sprintf("    gamma/((1+g)(s-1)):           %.3f  (Soderbery: 0.217)\n",
            median(cr$gamma / ((1 + cr$gamma) * (cr$sigma - 1)))))
cat(sprintf("    Convergence (code=0):         %.1f%%\n",
            100 * mean(cr$convergence == 0)))

# --- Regional structural ratios for comparison ---
rr <- regional_results[!is.na(sigma) & !is.na(gamma) & gamma > 0]
cat("\n  Structural ratios (regional):\n")
cat(sprintf("    gamma/(1+gamma):              %.3f\n",
            median(rr$gamma / (1 + rr$gamma))))
cat(sprintf("    1/(sigma-1):                  %.3f\n",
            median(1 / (rr$sigma - 1))))
cat(sprintf("    gamma/((1+g)(s-1)):           %.3f\n",
            median(rr$gamma / ((1 + rr$gamma) * (rr$sigma - 1)))))

# --- Tier distribution ---
if ("tier" %in% names(cr)) {
  cat("\n  Tier distribution (country):\n")
  tier_tab <- cr[, .N, by = tier]
  setorder(tier_tab, tier)
  for (i in seq_len(nrow(tier_tab))) {
    lbl <- switch(as.character(tier_tab$tier[i]),
                  "0" = "Reference exporter",
                  "1" = "Full (import+export)",
                  "2" = "Import-side only",
                  "3" = "Assigned from regional",
                  paste("Tier", tier_tab$tier[i]))
    cat(sprintf("    Tier %s (%s): %s (%.1f%%)\n",
                tier_tab$tier[i], lbl,
                format(tier_tab$N[i], big.mark = ","),
                100 * tier_tab$N[i] / nrow(cr)))
  }
}

# --- Output files ---
cat("\n  Output files:\n")
cat(sprintf("    Stage 1:  %s\n", sigma_file))
cat(sprintf("    Stage 2a: %s\n", regional_file))
cat(sprintf("    Stage 2a: %s_summary.rds\n", out_base_regional))
cat(sprintf("    Stage 2b: %s\n", country_file))
cat(sprintf("    Stage 2b: %s_summary.rds\n", out_base_country))

cat("\nDone.\n")
