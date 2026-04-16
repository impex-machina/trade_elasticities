#' ===========================================================================
#' feen94_het_baci.R
#'
#' Core function library for Soderbery (2018) heterogeneous elasticity
#' estimation, adapted for CEPII BACI data.
#'
#' This file contains ALL shared functions. It does NOT contain config
#' or estimation calls. Source this, then use:
#'   - run_estimation.R   (three-stage pipeline: Feenstra sigma ->
#'                         regional gamma -> country gamma)
#'
#' CITATION:
#'   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
#'   Tariffs," JIE, 114, 2018, pp. 44-62.
#'
#' Last updated: 2026-04-16 (patched: exp_lookup perf, V-lookups,
#'   period-floor weights, Tier 3 separation, plateau fallback,
#'   checkpointing, lambda calibration diagnostics, lint)
#' ===========================================================================

library(data.table)
library(parallel)

# ===========================================================================
#  OBJECTIVE FUNCTION LOADING
#
#  Try Rcpp version first (5-10x faster). Fall back to pure R if Rcpp
#  is unavailable or compilation fails (e.g., no Rtools on Windows).
#  Both produce identical results — Rcpp is a speed optimization only.
# ===========================================================================

.het_obj_rcpp_loaded <- FALSE

tryCatch({
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    cpp_file <- file.path(getwd(), "het_obj_rcpp.cpp")
    if (file.exists(cpp_file)) {
      Rcpp::sourceCpp(cpp_file)
      # het_obj_rcpp is now in the global environment
      het_obj <- het_obj_rcpp
      .het_obj_rcpp_loaded <- TRUE
      cat("  Objective function: Rcpp (compiled C++)\n")
    }
  }
}, error = function(e) {
  cat(sprintf("  Rcpp compilation failed: %s\n", conditionMessage(e)))
})

if (!.het_obj_rcpp_loaded) {
  source("het_obj.R")
  cat("  Objective function: pure R\n")
}

cat("HETEROGENEOUS IMPORT DEMAND AND EXPORT SUPPLY ELASTICITIES\n")
cat("  Soderbery (2018), JIE 114, pp. 44-62.\n")
cat("  Applied to CEPII BACI data.\n")
cat("  PLEASE CITE ACCORDINGLY.\n\n")

# ===========================================================================
#  REGION MAPPING
#
#  Soderbery (2018) Table 1 keeps 13 large countries as individual units
#  and groups remaining countries into 7 regions following UN M49
#  (https://unstats.un.org/unsd/methods/m49/m49regin.htm).
#
#  BACI uses Comtrade numeric country codes. This mapping covers the
#  major trading nations. Unmapped codes are assigned to "OTHER".
#  Users should verify against their BACI country_codes metadata file
#  and adjust as needed.
# ===========================================================================

build_region_map <- function() {

  # --- 13 individual countries (BACI/Comtrade numeric codes) ---
  # Verified against BACI country_codes_V202601.csv
  individual <- data.table(
    cty_code = c(36L, 76L, 124L, 156L, 344L, 446L, 276L, 251L, 826L,
                 699L, 380L, 392L, 484L, 643L, 842L),
    region   = c("AUS","BRA","CAN","CHN","CHN","CHN","DEU","FRA","GBR",
                 "IND","ITA","JPN","MEX","RUS","USA")
  )
  # Notes: 380 = Italy (not 381); 344 = Hong Kong, 446 = Macao → CHN
  #        Historical codes 278/280 (German states) → DEU below

  # --- African (AFR) ---
  afr_codes <- c(
    12L,24L,204L,72L,854L,108L,120L,132L,140L,148L,174L,178L,180L,
    262L,818L,226L,232L,231L,266L,270L,288L,324L,624L,384L,404L,
    426L,430L,434L,450L,454L,466L,478L,480L,504L,508L,516L,562L,
    566L,646L,678L,686L,690L,694L,706L,710L,728L,729L,736L,748L,
    834L,768L,788L,800L,894L,716L,
    # Additional BACI codes
    175L,  # Mayotte
    654L,  # Saint Helena
    711L   # Southern African Customs Union (...1999)
  )

  # --- Asian (ASA, excl CHN/HK/Macao, IND, JPN) ---
  asa_codes <- c(
    4L,48L,50L,64L,96L,104L,116L,360L,364L,368L,376L,400L,398L,
    414L,417L,418L,422L,458L,462L,496L,524L,512L,586L,608L,634L,
    682L,702L,410L,144L,760L,762L,764L,626L,795L,784L,860L,
    704L,887L,
    # Additional BACI codes
    275L,  # State of Palestine
    408L,  # Dem. People's Rep. of Korea
    490L   # Other Asia, nes (includes Taiwan in Comtrade)
  )

  # --- Caribbean (CAR) ---
  car_codes <- c(
    28L,44L,52L,84L,192L,212L,214L,308L,332L,388L,659L,662L,670L,
    740L,780L,535L,
    # Additional BACI codes
    60L,   # Bermuda
    92L,   # British Virgin Islands
    136L,  # Cayman Islands
    500L,  # Montserrat
    530L,  # Netherlands Antilles (...2010)
    531L,  # Curacao
    533L,  # Aruba
    534L,  # Saint Maarten
    652L,  # Saint Barthelemy
    660L,  # Anguilla
    796L   # Turks and Caicos
  )

  # --- Northern/Western Europe (NWU, excl DEU, FRA, GBR) ---
  nwu_codes <- c(
    40L,56L,208L,233L,246L,352L,372L,428L,440L,442L,528L,579L,
    752L,757L,724L,620L,300L,470L,
    # Additional BACI codes
    20L,   # Andorra
    58L,   # Belgium-Luxembourg (...1998)
    292L,  # Gibraltar
    304L,  # Greenland
    666L,  # Saint Pierre and Miquelon
    697L   # Europe EFTA, nes
  )
  # Notes: 579 = Norway (not 578); 757 = Switzerland (not 756)

  # --- Oceania (OCE, excl AUS) ---
  oce_codes <- c(
    554L,598L,242L,90L,882L,776L,548L,583L,584L,585L,520L,
    # Additional BACI codes
    16L,   # American Samoa
    162L,  # Christmas Islands
    166L,  # Cocos Islands
    184L,  # Cook Islands
    258L,  # French Polynesia
    296L,  # Kiribati
    540L,  # New Caledonia
    570L,  # Niue
    574L,  # Norfolk Islands
    580L,  # N. Mariana Islands
    772L,  # Tokelau
    798L,  # Tuvalu
    849L,  # US Misc. Pacific Islands
    876L   # Wallis and Futuna
  )

  # --- South American (SAM, excl BRA, MEX) ---
  sam_codes <- c(
    32L,68L,152L,170L,188L,218L,222L,320L,328L,340L,558L,591L,
    600L,604L,858L,862L,
    # Additional BACI codes
    238L   # Falkland Islands
  )

  # --- Southern/Eastern Europe (SEU, excl ITA, RUS) ---
  seu_codes <- c(
    8L,70L,100L,191L,196L,203L,268L,348L,498L,499L,616L,642L,
    688L,703L,705L,792L,804L,807L,112L,51L,31L,
    # Additional BACI codes
    200L,  # Czechoslovakia (...1992)
    674L,  # San Marino
    891L   # Serbia and Montenegro (...2005)
  )

  # --- Historical German codes → DEU ---
  hist_deu <- data.table(
    cty_code = c(278L, 280L),  # Dem. Rep. / Fed. Rep. of Germany (...1990)
    region   = c("DEU", "DEU")
  )

  regions <- rbindlist(list(
    data.table(cty_code = afr_codes, region = "AFR"),
    data.table(cty_code = asa_codes, region = "ASA"),
    data.table(cty_code = car_codes, region = "CAR"),
    data.table(cty_code = nwu_codes, region = "NWU"),
    data.table(cty_code = oce_codes, region = "OCE"),
    data.table(cty_code = sam_codes, region = "SAM"),
    data.table(cty_code = seu_codes, region = "SEU"),
    hist_deu
  ))

  # Remove any duplicates (a code appearing in both individual + region)
  regions <- regions[!cty_code %in% individual$cty_code]
  regions <- unique(regions, by = "cty_code")

  rbindlist(list(individual, regions))
}


#' Assign regions to a vector of BACI country codes.
#' Unmapped codes are assigned to "OTHER".
#'
#' @param codes Integer or character vector of BACI country codes.
#' @param custom_map Optional data.table with columns (cty_code, region)
#'   to override the built-in mapping.
#' @return Character vector of region labels, same length as codes.
assign_regions <- function(codes, custom_map = NULL) {
  rmap <- if (!is.null(custom_map)) custom_map else build_region_map()
  rmap[, cty_code := as.integer(cty_code)]
  lookup <- data.table(cty_code = as.integer(codes))
  merged <- rmap[lookup, on = "cty_code"]
  merged[is.na(region), region := "OTHER"]
  merged$region
}


# ===========================================================================
#  HS CODE UTILITIES
# ===========================================================================

#' Pad HS6 codes to ensure 6-digit strings with leading zeroes.
#'
#' BACI's product column (k) may be read as numeric by fread/read.csv,
#' which strips leading zeroes (e.g., 010110 -> 10110). This function
#' detects and fixes the issue.
#'
#' @param k Vector of HS6 codes (character or numeric).
#' @return Character vector of 6-digit zero-padded HS6 codes.
pad_hs6 <- function(k) {
  k <- as.character(k)
  # Detect if padding is needed: any code shorter than 6 digits
  needs_pad <- nchar(k) < 6L
  if (any(needs_pad)) {
    n_padded <- sum(needs_pad)
    cat(sprintf("  [HS6 padding] %s codes shorter than 6 digits; ",
                format(n_padded, big.mark = ",")))
    cat("padding with leading zeroes.\n")
    k <- formatC(as.integer(k), width = 6, format = "d", flag = "0")
  }
  k
}


#' Extract HS4 heading from HS6 codes.
#' @param hs6 Character vector of 6-digit HS6 codes.
#' @return Character vector of 4-digit HS4 codes.
hs6_to_hs4 <- function(hs6) {
  substr(hs6, 1, 4)
}


# ===========================================================================
#  CONFIG VALIDATION
# ===========================================================================

#' Validate config for obvious misconfigurations before expensive operations.
#'
#' Checks that required fields exist, types are correct, and values are
#' logically consistent. Stops with an informative error on failure.
#' Warns on likely-problematic but non-fatal settings.
#'
#' @param cfg Config list.
validate_config <- function(cfg) {

  # --- Required fields ---
  required <- c("filepath", "value", "quan", "good", "importer", "exporter",
                "time", "minyear", "agg_level", "use_regions",
                "min_exporters", "min_destinations", "min_periods",
                "sigma_start", "gamma_start", "sigma_V_default",
                "gamma_V_default", "tail_trim_pct")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0L) {
    stop("Config missing required fields: ", paste(missing, collapse = ", "))
  }

  # --- Filepath ---
  if (!file.exists(cfg$filepath)) {
    stop("Data path does not exist: ", cfg$filepath)
  }

  # --- Aggregation level ---
  if (!cfg$agg_level %in% c("hs4", "hs6")) {
    stop("agg_level must be 'hs4' or 'hs6', got: ", cfg$agg_level)
  }

  # --- Year range ---
  if (!is.numeric(cfg$minyear) || cfg$minyear < 1900 || cfg$minyear > 2100) {
    stop("minyear must be a reasonable year, got: ", cfg$minyear)
  }
  if (!is.null(cfg$maxyear) && !is.na(cfg$maxyear)) {
    if (!is.numeric(cfg$maxyear) || cfg$maxyear < cfg$minyear) {
      stop("maxyear must be >= minyear (", cfg$minyear, "), got: ", cfg$maxyear)
    }
    year_span <- cfg$maxyear - cfg$minyear + 1L
    if (year_span < cfg$min_periods + 1L) {
      warning(sprintf(paste("Year range (%d-%d = %d years) may be too short",
                            "for min_periods=%d (need %d+ years of data",
                            "after first-differencing)."),
                      cfg$minyear, cfg$maxyear, year_span,
                      cfg$min_periods, cfg$min_periods + 1L))
    }
  }

  # --- Structural defaults ---
  if (cfg$sigma_V_default <= 1) {
    stop("sigma_V_default must be > 1, got: ", cfg$sigma_V_default)
  }
  if (cfg$gamma_V_default <= 0) {
    stop("gamma_V_default must be > 0, got: ", cfg$gamma_V_default)
  }
  if (cfg$sigma_start <= 1) {
    stop("sigma_start must be > 1, got: ", cfg$sigma_start)
  }
  if (cfg$gamma_start <= 0) {
    stop("gamma_start must be > 0, got: ", cfg$gamma_start)
  }

  # --- Filtering ---
  if (cfg$min_exporters < 1L) {
    stop("min_exporters must be >= 1, got: ", cfg$min_exporters)
  }
  if (cfg$min_periods < 2L) {
    warning("min_periods < 2 means single-observation cells may be estimated. ",
            "This is unlikely to produce reliable estimates.")
  }

  # --- Trimming ---
  if (!is.na(cfg$tail_trim_pct) && (cfg$tail_trim_pct < 0 || cfg$tail_trim_pct >= 0.5)) {
    stop("tail_trim_pct must be in [0, 0.5), got: ", cfg$tail_trim_pct)
  }

  invisible(TRUE)
}


# ===========================================================================
#  DATA QUALITY TRACKER
# ===========================================================================

new_quality_log <- function() {
  env <- new.env(parent = emptyenv())
  env$steps <- list()
  env$add <- function(stage, n_obs, n_dropped = NA_integer_,
                      trade_value = NA_real_, detail = "") {
    entry <- list(stage = stage, n_obs = n_obs,
                  n_dropped = n_dropped, trade_value = trade_value,
                  detail = detail)
    env$steps[[length(env$steps) + 1L]] <- entry
    invisible(NULL)
  }
  env
}

print_quality_log <- function(qlog) {
  cat("\n")
  cat("=====================================================================\n")
  cat("  DATA QUALITY REPORT\n")
  cat("=====================================================================\n\n")

  has_tv <- any(sapply(qlog$steps, function(s) !is.na(s$trade_value)))

  if (has_tv) {
    cat(sprintf("  %-40s %12s %12s %16s\n", "Stage", "Obs", "Dropped", "Trade Val ($B)"))
    cat(sprintf("  %-40s %12s %12s %16s\n",
                paste(rep("-", 40), collapse = ""),
                paste(rep("-", 12), collapse = ""),
                paste(rep("-", 12), collapse = ""),
                paste(rep("-", 16), collapse = "")))
  } else {
    cat(sprintf("  %-45s %12s %12s\n", "Stage", "Obs", "Dropped"))
    cat(sprintf("  %-45s %12s %12s\n",
                paste(rep("-", 45), collapse = ""),
                paste(rep("-", 12), collapse = ""),
                paste(rep("-", 12), collapse = "")))
  }

  for (s in qlog$steps) {
    n_str <- format(s$n_obs, big.mark = ",")
    d_str <- if (is.na(s$n_dropped)) "" else format(s$n_dropped, big.mark = ",")
    if (has_tv) {
      tv_str <- if (is.na(s$trade_value)) "" else sprintf("%.1f", s$trade_value / 1e6)
      cat(sprintf("  %-40s %12s %12s %16s\n", s$stage, n_str, d_str, tv_str))
    } else {
      cat(sprintf("  %-45s %12s %12s\n", s$stage, n_str, d_str))
    }
    if (nchar(s$detail) > 0) {
      cat(sprintf("    %s\n", s$detail))
    }
  }

  # Overall retention
  n_start <- qlog$steps[[1]]$n_obs
  n_end   <- qlog$steps[[length(qlog$steps)]]$n_obs
  pct <- round(100 * n_end / n_start, 1)
  cat(sprintf("\n  Overall retention: %s / %s (%.1f%%)\n",
              format(n_end, big.mark = ","),
              format(n_start, big.mark = ","), pct))

  if (has_tv) {
    tv_start <- qlog$steps[[1]]$trade_value
    tv_end   <- qlog$steps[[length(qlog$steps)]]$trade_value
    if (!is.na(tv_start) && !is.na(tv_end) && tv_start > 0) {
      cat(sprintf("  Trade value retention: $%.1fB / $%.1fB (%.1f%%)\n",
                  tv_end / 1e6, tv_start / 1e6, 100 * tv_end / tv_start))
    }
  }
  cat("=====================================================================\n\n")
}


# ===========================================================================
#  HELPER FUNCTIONS
# ===========================================================================

#' Choose reference exporter within an import market.
#' Selects the largest, most persistent exporter.
choose_reference <- function(dt) {
  stats <- dt[, .(n_periods = uniqueN(t),
                   total_value = sum(cusval, na.rm = TRUE)),
              by = exporter]
  max_pd <- max(stats$n_periods)
  candidates <- stats[n_periods >= max_pd]
  candidates$exporter[which.max(candidates$total_value)]
}


#' Compute BW weights (paper p. 50, fn 14).
#' Weight = T^(3/2) * (1/x_t + 1/x_{t-1})^(-1/2)
bw_weight <- function(cusval_t, cusval_lag, T_count) {
  w <- T_count^1.5 * (1 / cusval_t + 1 / cusval_lag)^(-0.5)
  w[is.na(w) | !is.finite(w)] <- 1
  w
}


# ===========================================================================
#  DATA LOADING
# ===========================================================================

#' Load BACI data from either a single CSV/RDS or a directory of per-year CSVs.
#'
#' BACI is distributed as one CSV per year inside a zip archive.
#' This function handles both cases.
load_baci <- function(filepath) {
  if (dir.exists(filepath)) {
    # Match only BACI trade data files (BACI_HS*_Y####_V*.csv),
    # excluding metadata files like country_codes*.csv and product_codes*.csv
    csv_files <- list.files(filepath, pattern = "BACI_HS.*_Y\\d{4}_V.*\\.csv$",
                            full.names = TRUE, recursive = FALSE)
    if (length(csv_files) == 0L) {
      # Fallback: try all CSVs but warn
      csv_files <- list.files(filepath, pattern = "\\.csv$",
                              full.names = TRUE, recursive = FALSE)
      warning("No files matching BACI_HS*_Y*_V*.csv pattern found. ",
              "Loading all ", length(csv_files), " CSV files in directory. ",
              "Remove metadata CSVs (country_codes, product_codes) ",
              "from the directory to avoid errors.")
    }
    if (length(csv_files) == 0L) stop("No CSV files found in: ", filepath)
    cat(sprintf("  Loading %d BACI trade data files from: %s\n",
                length(csv_files), filepath))
    dt_list <- lapply(csv_files, fread, colClasses = list(character = "k"))
    raw <- rbindlist(dt_list)
  } else if (grepl("\\.csv$", filepath, ignore.case = TRUE)) {
    raw <- fread(filepath, colClasses = list(character = "k"))
  } else if (grepl("\\.rds$", filepath, ignore.case = TRUE)) {
    raw <- as.data.table(readRDS(filepath))
  } else if (grepl("\\.dta$", filepath, ignore.case = TRUE)) {
    if (!requireNamespace("haven", quietly = TRUE))
      stop("Package 'haven' required to read .dta files.")
    raw <- as.data.table(haven::read_dta(filepath))
  } else {
    stop("Unsupported file type: ", filepath)
  }
  raw
}


# ===========================================================================
#  ESTIMATE ONE (IMPORTER, PRODUCT) PAIR
# ===========================================================================

#' Lightweight failure indicator for cell-level diagnostics.
#' Returned instead of NULL so that estimate_product can log the reason.
cell_failure <- function(reason) {
  structure(list(reason = reason), class = "cell_failure")
}


#' Compute exporter weights for weighted least squares.
#'
#' Returns per-exporter weights for the objective function sum across
#' exporters. Two modes:
#'   "uniform"     — rep(1, J). Equal weight per exporter.
#'   "trade_value" — weighted by total cusval within the cell,
#'                   normalized to sum to J. Large bilateral
#'                   relationships contribute more to the objective.
#'
#' PERIOD FLOOR (trade_value mode only):
#'   An exporter's trade-value weight is further multiplied by
#'   min(n_periods, period_floor_ref) / period_floor_ref, so exporters
#'   with very short panels (e.g., a single-year large shipment — aircraft,
#'   rough diamonds, LNG cargo) get downweighted relative to persistent
#'   relationships with the same total value. Default reference is 10 years.
#'   Set cfg$weight_period_floor to disable (0) or use a different reference.
#'
#' @param dt_nonref Non-reference exporter data for the focal (importer, product).
#' @param exporter_order Character vector of exporters aligned with moments.
#' @param cfg Config list. Reads cfg$exporter_weight (default "uniform"),
#'   cfg$weight_period_floor (default 10).
#' @return Numeric vector of length(exporter_order).
compute_exporter_weights <- function(dt_nonref, exporter_order, cfg) {
  mode <- if (!is.null(cfg$exporter_weight)) cfg$exporter_weight else "uniform"
  J <- length(exporter_order)
  if (J == 0L) return(numeric(0))

  if (mode == "trade_value") {
    # Aggregate trade value AND period count per exporter in one pass
    stats <- dt_nonref[, .(wt = sum(cusval, na.rm = TRUE),
                            n_per = uniqueN(t)), by = exporter]
    setkey(stats, exporter)
    matched <- stats[J(exporter_order)]
    w     <- matched$wt
    n_per <- matched$n_per

    # Period-count floor: exporters with < period_floor_ref periods
    # get their weight scaled by n_periods / period_floor_ref.
    period_floor_ref <- if (!is.null(cfg$weight_period_floor))
      cfg$weight_period_floor else 10L
    if (period_floor_ref > 0L) {
      adj <- pmin(n_per, period_floor_ref) / period_floor_ref
      adj[is.na(adj)] <- 0
      w <- w * adj
    }

    tot <- sum(w, na.rm = TRUE)
    if (is.finite(tot) && tot > 0) {
      w <- w * J / tot
    }
    w[is.na(w) | !is.finite(w) | w <= 0] <- 1
    w
  } else {
    rep(1, J)
  }
}


#' Build export-side moment matrices (Eq. 11).
#'
#' For each non-reference exporter j, constructs the BW-weighted
#' time-averaged export-side moments by looking at j's export
#' shares across destinations. Returns the matrices, index map,
#' and reference-destination elasticities needed by the objective.
#'
#' Extracted from estimate_importer_product for clarity.
#' The computation is identical — this is a structural refactor only.
#'
#' @param exporter_order Character vector of non-reference exporters.
#' @param focal_importer The importer being estimated.
#' @param all_dt Full product-level data.table (all importers).
#'   Used only if exp_lookup is NULL.
#' @param cfg Config list. Reads cfg$sigma_V_default, cfg$gamma_V_default,
#'   and optionally cfg$sigma_V_lookup, cfg$gamma_V_lookup.
#' @param exp_lookup Optional named list of per-exporter data.tables
#'   (produced by compute_exporter_lookup). If provided, avoids repeated
#'   O(n) filtering of all_dt — key perf optimization for HS6.
#'   exp_lookup[[exp_j]] returns the subset of all_dt where exporter == exp_j.
#' @return Named list with exp_Y, exp_X, exp_jmap, sig_V, gam_V, wt_exp.
#'
#' REFERENCE-DESTINATION ELASTICITY LOOKUP:
#'   Soderbery's Eq. (11) parameterizes the export-side moments by
#'   (sigma_V, gamma_V) at the reference destination V. The baseline
#'   treatment uses a single global default for all exporters.
#'
#'   When cfg$sigma_V_lookup / cfg$gamma_V_lookup are provided, each
#'   Tier 1 exporter's reference destination gets its own (sigma, gamma)
#'   drawn from Stage 1 / Stage 2a estimates:
#'     - sigma_V_lookup: data.table(importer, good, sigma) — Stage 1 output
#'     - gamma_V_lookup: data.table(importer, good, gamma) — Stage 2a output
#'         (at Stage 2b; at Stage 2a, gamma_V_lookup is typically NULL
#'         and the global default is used)
#'   Lookup falls back to cfg$sigma_V_default / cfg$gamma_V_default when
#'   no match is found (thin bilateral relationships, non-overlapping
#'   product coverage, etc.).
build_export_moments <- function(exporter_order, focal_importer, all_dt, cfg,
                                  exp_lookup = NULL) {

  exp_Y_list   <- list()
  exp_X_list   <- list()
  exp_jmap_vec <- integer(0)
  sig_V_vec    <- numeric(0)
  gam_V_vec    <- numeric(0)

  # Helpers: lookup (importer, good) -> (sigma, gamma) with fallback.
  # Indexed by importer for O(log n) lookup via setkey; caller is
  # responsible for having keyed these lookups once upstream.
  sig_V_lkp <- cfg$sigma_V_lookup
  gam_V_lkp <- cfg$gamma_V_lookup

  for (j_idx in seq_along(exporter_order)) {
    exp_j <- exporter_order[j_idx]

    # Use pre-split lookup if available (HS6 perf optimization);
    # else fall back to O(n) filter of all_dt.
    if (!is.null(exp_lookup)) {
      exp_flows <- exp_lookup[[exp_j]]
      if (is.null(exp_flows) || nrow(exp_flows) == 0L) next
    } else {
      exp_flows <- all_dt[exporter == exp_j]
    }
    n_dest <- uniqueN(exp_flows$importer)
    if (n_dest < cfg$min_destinations) next

    dest_stats <- exp_flows[importer != focal_importer,
                            .(dest_val = sum(cusval, na.rm = TRUE)),
                            by = importer]
    if (nrow(dest_stats) == 0L) next

    ref_dest <- dest_stats$importer[which.max(dest_stats$dest_val)]

    ref_dest_vals <- exp_flows[importer == ref_dest,
                               .(t, ref_lp_exp = lp_dif,
                                 ref_ls_exp = ls_exp_dif)]

    focal_vals <- exp_flows[importer == focal_importer]
    focal_vals <- ref_dest_vals[focal_vals, on = "t"]
    focal_vals <- focal_vals[!is.na(ref_lp_exp) & !is.na(ref_ls_exp) &
                             !is.na(lp_dif) & !is.na(ls_exp_dif)]

    if (nrow(focal_vals) == 0L) next

    focal_vals[, `:=`(
      exp_y  = (lp_dif - ref_lp_exp)^2,
      exp_x1 = ls_exp_dif^2,
      exp_x2 = ls_exp_dif * lp_dif,
      exp_x3 = ref_ls_exp^2,
      exp_x4 = ref_ls_exp * lp_dif,
      exp_x5 = ref_ls_exp * ref_lp_exp,
      exp_x6 = ls_exp_dif * ref_lp_exp,
      exp_x7 = ls_exp_dif * ref_ls_exp,
      exp_x8 = ref_lp_exp^2,
      exp_x9 = ref_lp_exp * lp_dif
    )]
    focal_vals <- focal_vals[!is.na(exp_y)]
    if (nrow(focal_vals) == 0L) next

    setorder(focal_vals, t)
    focal_vals[, `:=`(cusval_lag = shift(cusval, 1L), pd_e = .N)]
    focal_vals[, bw_w_e := bw_weight(cusval, cusval_lag, pd_e)]

    exp_cols <- c("exp_y","exp_x1","exp_x2","exp_x3","exp_x4",
                  "exp_x5","exp_x6","exp_x7","exp_x8","exp_x9")
    exp_mom <- focal_vals[,
      lapply(.SD, weighted.mean, w = bw_w_e, na.rm = TRUE),
      .SDcols = exp_cols
    ]

    exp_Y_list[[length(exp_Y_list) + 1L]] <- exp_mom$exp_y
    exp_X_list[[length(exp_X_list) + 1L]] <- as.numeric(
      exp_mom[, .(exp_x1,exp_x2,exp_x3,exp_x4,exp_x5,
                  exp_x6,exp_x7,exp_x8,exp_x9)]
    )
    exp_jmap_vec <- c(exp_jmap_vec, j_idx + 2L)

    # --- Per-ref-destination (sig_V, gam_V) lookup with fallback ---
    # g_code is the product for this product-level estimation call.
    # all_dt$good[1] holds it (all rows share a single product).
    g_code <- all_dt$good[1]

    sig_V_val <- cfg$sigma_V_default
    if (!is.null(sig_V_lkp)) {
      row_s <- sig_V_lkp[importer == ref_dest & good == g_code]
      if (nrow(row_s) > 0L && !is.na(row_s$sigma[1]) && row_s$sigma[1] > 1) {
        sig_V_val <- row_s$sigma[1]
      }
    }

    gam_V_val <- cfg$gamma_V_default
    if (!is.null(gam_V_lkp)) {
      row_g <- gam_V_lkp[importer == ref_dest & good == g_code]
      if (nrow(row_g) > 0L && !is.na(row_g$gamma[1]) && row_g$gamma[1] > 0) {
        gam_V_val <- row_g$gamma[1]
      }
    }

    sig_V_vec    <- c(sig_V_vec, sig_V_val)
    gam_V_vec    <- c(gam_V_vec, gam_V_val)
  }

  M <- length(exp_Y_list)

  if (M > 0L) {
    list(
      exp_Y  = unlist(exp_Y_list),
      exp_X  = matrix(unlist(exp_X_list), nrow = M, ncol = 9, byrow = TRUE),
      jmap   = exp_jmap_vec,
      sig_V  = sig_V_vec,
      gam_V  = gam_V_vec,
      wt_exp = rep(1, M),
      M      = M
    )
  } else {
    list(
      exp_Y  = numeric(0),
      exp_X  = matrix(nrow = 0, ncol = 9),
      jmap   = integer(0),
      sig_V  = numeric(0),
      gam_V  = numeric(0),
      wt_exp = numeric(0),
      M      = 0L
    )
  }
}


#' Pre-split a product-level data.table by exporter.
#'
#' Constructs a named list of per-exporter data.table slices. Passing this
#' to build_export_moments converts O(N_exporters) repeated filtering of
#' all_dt into O(1) list indexing. At HS6 scale (5,300 products, many
#' exporters per product) this is a meaningful performance win.
#'
#' Called ONCE per product (not per importer) in estimate_product /
#' estimate_product_fixed_sigma.
#'
#' @param dt_g Product-level data.table.
#' @return Named list; names are exporter codes, values are their subsets.
compute_exporter_lookup <- function(dt_g) {
  split(dt_g, by = "exporter", keep.by = TRUE)
}


#' Estimate one (importer, product) cell.
#'
#' Returns a data.table on success, or a cell_failure object with
#' a diagnostic reason on failure. The cell_failure is collected
#' by estimate_product for the failure log.
#'
#' @param exp_lookup Optional pre-split per-exporter lookup (see
#'   compute_exporter_lookup). Passed through to build_export_moments.
estimate_importer_product <- function(imp_dt, focal_importer, all_dt, cfg,
                                      exp_lookup = NULL) {

  dt <- imp_dt[importer == focal_importer]

  n_exp <- uniqueN(dt$exporter)
  max_pd <- dt[, max(period_count)]

  if (n_exp < cfg$min_exporters || max_pd < cfg$min_periods) {
    return(cell_failure("insufficient_data"))
  }

  # --- Choose reference exporter k ---
  ref_exporter <- choose_reference(dt)

  # --- Reference-differencing (import side) ---
  ref_vals <- dt[exporter == ref_exporter,
                 .(t, ref_ls_dif = ls_imp_dif, ref_lp_dif = lp_dif)]

  dt <- ref_vals[dt, on = "t"]
  dt <- dt[!is.na(ref_ls_dif) & !is.na(ref_lp_dif)]

  dt[, `:=`(Dk_lp = lp_dif - ref_lp_dif,
            Dk_ls = ls_imp_dif - ref_ls_dif)]

  # --- Import-side moment variables (eq 10) ---
  dt[, `:=`(
    imp_y  = Dk_lp^2,
    imp_x1 = Dk_ls^2,
    imp_x2 = Dk_ls * Dk_lp,
    imp_x3 = Dk_ls * lp_dif,
    imp_x4 = Dk_ls * ref_lp_dif,
    imp_x5 = Dk_lp * ref_lp_dif
  )]

  dt <- dt[!is.na(imp_y) & !is.na(imp_x1) & !is.na(imp_x2)]

  dt_nonref <- dt[exporter != ref_exporter]
  if (nrow(dt_nonref) == 0L) return(cell_failure("no_nonref_exporters"))

  # --- BW weights ---
  setorder(dt_nonref, exporter, t)
  dt_nonref[, cusval_lag := shift(cusval, 1L), by = exporter]
  dt_nonref[, bw_w := bw_weight(cusval, cusval_lag, period_count)]

  # --- Time-average with BW weights ---
  imp_moments <- dt_nonref[,
    lapply(.SD, weighted.mean, w = bw_w, na.rm = TRUE),
    by = exporter,
    .SDcols = c("imp_y", "imp_x1", "imp_x2", "imp_x3", "imp_x4", "imp_x5")
  ]
  setorder(imp_moments, exporter)

  J <- nrow(imp_moments)
  if (J < 1L) return(cell_failure("no_valid_moments"))

  exporter_order <- imp_moments$exporter
  imp_Y_vec <- imp_moments$imp_y
  imp_X_mat <- as.matrix(imp_moments[, .(imp_x1, imp_x2, imp_x3,
                                          imp_x4, imp_x5)])
  wt_imp_vec <- rep(1, J)


  # =========================================================
  #  EXPORT-SIDE MOMENT VARIABLES (eq 11)
  # =========================================================

  exp_mom <- build_export_moments(exporter_order, focal_importer, all_dt, cfg,
                                   exp_lookup = exp_lookup)
  exp_Y_vec    <- exp_mom$exp_Y
  exp_X_mat    <- exp_mom$exp_X
  exp_jmap_vec <- exp_mom$jmap
  sig_V_vec    <- exp_mom$sig_V
  gam_V_vec    <- exp_mom$gam_V
  wt_exp_vec   <- exp_mom$wt_exp


  # =========================================================
  #  JOINT NONLINEAR SUR ESTIMATION
  # =========================================================

  # Per-cell starting values from regional estimates if available
  sig_init <- cfg$sigma_start
  gam_init <- cfg$gamma_start

  if (!is.null(cfg$regional_starts)) {
    g_code <- imp_dt$good[1]
    imp_region <- focal_importer
    # If running country-level, map country code to region
    if (!is.null(cfg$regional_starts_rmap)) {
      imp_region <- assign_regions(focal_importer, cfg$regional_starts_rmap)
    }
    match_row <- cfg$regional_starts[region == imp_region & good == g_code]
    if (nrow(match_row) > 0L) {
      sig_init <- match_row$sigma[1]
      gam_init <- match_row$gamma[1]
    }
  }

  d_start <- c(sig_init, gam_init, rep(gam_init, J))
  lower_bounds <- c(1 + 1e-6, rep(1e-6, J + 1))

  result <- tryCatch(
    optim(
      par    = d_start,
      fn     = het_obj,
      method = "L-BFGS-B",
      lower  = lower_bounds,
      upper  = rep(Inf, J + 2),
      imp_Y  = imp_Y_vec,  imp_X = imp_X_mat,
      exp_Y  = exp_Y_vec,  exp_X = exp_X_mat,
      exp_jmap  = exp_jmap_vec,
      exp_sig_V = sig_V_vec, exp_gam_V = gam_V_vec,
      wt_imp = wt_imp_vec, wt_exp = wt_exp_vec,
      control = list(maxit = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(result) || result$convergence != 0) {
    result <- tryCatch(
      optim(
        par = d_start, fn = het_obj, method = "Nelder-Mead",
        imp_Y = imp_Y_vec, imp_X = imp_X_mat,
        exp_Y = exp_Y_vec, exp_X = exp_X_mat,
        exp_jmap = exp_jmap_vec,
        exp_sig_V = sig_V_vec, exp_gam_V = gam_V_vec,
        wt_imp = wt_imp_vec, wt_exp = wt_exp_vec,
        control = list(maxit = 1000)
      ),
      error = function(e) NULL
    )
  }

  if (is.null(result)) return(cell_failure("optimizer_failed"))

  d_hat <- result$par
  sigma_hat   <- d_hat[1]
  gamma_k_hat <- d_hat[2]
  gamma_j_hat <- d_hat[3:(J + 2)]

  if (sigma_hat <= 1) sigma_hat <- NA_real_
  gamma_k_hat <- max(gamma_k_hat, 0)
  gamma_j_hat <- pmax(gamma_j_hat, 0)

  data.table(
    importer     = focal_importer,
    exporter     = c(exporter_order, ref_exporter),
    sigma        = sigma_hat,
    gamma        = c(gamma_j_hat, gamma_k_hat),
    ref_exporter = ref_exporter,
    convergence  = result$convergence,
    obj_value    = result$value
  )
}


#' Optimal non-cooperative tariff (Proposition 1, paper p. 56).
optimal_tariff <- function(gamma, sigma, trade_values = NULL) {
  if (is.null(trade_values)) trade_values <- rep(1, length(gamma))
  ok <- gamma > 0 & !is.na(gamma) & trade_values > 0
  if (sum(ok) == 0L) return(NA_real_)
  g <- gamma[ok]; w <- trade_values[ok]
  num <- sum(w * g / (1 + g * sigma))
  den <- sum(w / (1 + g * sigma))
  if (den == 0) NA_real_ else num / den
}




# ===========================================================================
#  ESTIMATE ONE PRODUCT (all importers within a product)
# ===========================================================================

estimate_product <- function(g, dt_g, cfg) {
  t0 <- proc.time()["elapsed"]
  importers <- unique(dt_g$importer)
  results_g <- list()
  failures_g <- list()
  n_cells <- 0L; n_ok <- 0L; n_skipped <- 0L

  # --- Pre-filter: skip importers that cannot meet minimum requirements ---
  imp_stats <- dt_g[, .(n_exp = uniqueN(exporter),
                         max_pd = max(period_count)),
                     by = importer]
  viable_importers <- imp_stats[n_exp >= cfg$min_exporters &
                                max_pd >= cfg$min_periods, importer]
  n_skipped <- length(importers) - length(viable_importers)

  # --- Pre-split per-exporter slices ONCE per product (HS6 perf) ---
  exp_lookup <- compute_exporter_lookup(dt_g)

  for (imp in viable_importers) {
    n_cells <- n_cells + 1L
    est <- tryCatch(
      estimate_importer_product(dt_g, imp, dt_g, cfg, exp_lookup = exp_lookup),
      error = function(e) cell_failure(paste0("error: ", conditionMessage(e)))
    )

    if (inherits(est, "cell_failure")) {
      failures_g[[length(failures_g) + 1L]] <- list(
        importer = imp, good = g, reason = est$reason)
      next
    }

    if (!is.null(est)) {
      est[, good := g]
      trade_wt <- dt_g[importer == imp,
                       .(avg_trade = mean(cusval, na.rm = TRUE)), by = exporter]
      est <- trade_wt[est, on = "exporter"]
      ot <- optimal_tariff(est$gamma, est$sigma[1], est$avg_trade)
      est[, opt_tariff := ot]
      results_g[[length(results_g) + 1L]] <- est
      n_ok <- n_ok + 1L
    }
  }

  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  if (length(results_g) > 0L) {
    out <- rbindlist(results_g)
    attr(out, "timing") <- list(product = g, seconds = elapsed,
                                cells = n_cells, succeeded = n_ok,
                                skipped = n_skipped)
    attr(out, "failures") <- failures_g
    out
  } else { NULL }
}


# ===========================================================================
#  DATA PREPARATION PIPELINE
# ===========================================================================

#' Load and clean BACI data through HS4 aggregation (cacheable).
#'
#' Performs steps 1-5 of the data pipeline: load, clean, year filter,
#' HS4 aggregation. The result can be cached and reused across stages
#' to avoid reloading 270M+ rows multiple times.
#'
#' @param cfg Config list.
#' @return data.table with columns: year, good, importer, exporter, cusval, quantity.
prepare_raw_data <- function(cfg) {
  cat("Loading BACI data (raw prep for caching)...\n")
  raw <- load_baci(cfg$filepath)

  setnames(raw,
    old = c(cfg$value, cfg$quan, cfg$good, cfg$importer, cfg$exporter, cfg$time),
    new = c("cusval", "quantity", "good", "importer", "exporter", "year"),
    skip_absent = TRUE)

  raw <- raw[cusval > 0 & quantity > 0 & !is.na(cusval) & !is.na(quantity)]
  raw <- raw[year >= cfg$minyear]
  if (!is.null(cfg$maxyear) && !is.na(cfg$maxyear)) raw <- raw[year <= cfg$maxyear]

  raw[, good := pad_hs6(good)]
  raw <- raw[nchar(good) == 6L]

  if (cfg$agg_level == "hs4") {
    raw[, good := hs6_to_hs4(good)]
    raw <- raw[, .(cusval = sum(cusval), quantity = sum(quantity)),
               by = .(year, good, importer, exporter)]
  }

  cat(sprintf("  Raw data cached: %s obs, %d products, years %d-%d\n",
              format(nrow(raw), big.mark = ","), uniqueN(raw$good),
              min(raw$year), max(raw$year)))
  raw
}


#' Full data preparation pipeline.
#'
#' @param cfg Config list.
#' @param raw_cache Optional pre-loaded data.table from prepare_raw_data().
#'   If provided, skips the expensive loading/cleaning/aggregation steps.
#' @return Named list with dt (estimation data) and qlog (quality log).
prepare_data <- function(cfg, raw_cache = NULL) {

  validate_config(cfg)

  qlog <- new_quality_log()

  if (!is.null(raw_cache)) {
    cat("Using cached raw data (skipping load/clean/aggregate)...\n")
    raw <- copy(raw_cache)
    qlog$add("Raw data (from cache)", n_obs = nrow(raw),
             trade_value = sum(raw$cusval, na.rm = TRUE),
             detail = sprintf("%d products, %d importers, %d exporters",
                              uniqueN(raw$good), uniqueN(raw$importer),
                              uniqueN(raw$exporter)))
  } else {
    cat("Loading BACI data...\n")
    raw <- load_baci(cfg$filepath)

  setnames(raw,
    old = c(cfg$value, cfg$quan, cfg$good, cfg$importer, cfg$exporter, cfg$time),
    new = c("cusval", "quantity", "good", "importer", "exporter", "year"),
    skip_absent = TRUE)

  qlog$add("Raw data loaded", n_obs = nrow(raw),
           trade_value = sum(raw$cusval, na.rm = TRUE),
           detail = sprintf("%d products (HS6), %d importers, %d exporters, years %d-%d",
                            uniqueN(raw$good), uniqueN(raw$importer),
                            uniqueN(raw$exporter), min(raw$year), max(raw$year)))

  n_before <- nrow(raw)
  raw <- raw[cusval > 0 & quantity > 0 & !is.na(cusval) & !is.na(quantity)]
  qlog$add("Drop zero/missing value or quantity",
           n_obs = nrow(raw), n_dropped = n_before - nrow(raw),
           trade_value = sum(raw$cusval, na.rm = TRUE))

  n_before <- nrow(raw)
  raw <- raw[year >= cfg$minyear]
  if (!is.null(cfg$maxyear) && !is.na(cfg$maxyear)) {
    raw <- raw[year <= cfg$maxyear]
  }
  actual_maxyear <- if (!is.null(cfg$maxyear) && !is.na(cfg$maxyear)) cfg$maxyear else max(raw$year)
  qlog$add(sprintf("Keep years %d-%d", cfg$minyear, actual_maxyear),
           n_obs = nrow(raw), n_dropped = n_before - nrow(raw),
           trade_value = sum(raw$cusval, na.rm = TRUE))

  # HS6 -> HS4
  raw[, good := pad_hs6(good)]
  bad_len <- sum(nchar(raw$good) != 6L)
  if (bad_len > 0L) {
    warning(sprintf("%d HS6 codes not 6 digits after padding; dropping.", bad_len))
    raw <- raw[nchar(good) == 6L]
  }

  if (cfg$agg_level == "hs4") {
    raw[, good := hs6_to_hs4(good)]
    n_before <- nrow(raw)
    raw <- raw[, .(cusval = sum(cusval), quantity = sum(quantity)),
               by = .(year, good, importer, exporter)]
    qlog$add("Aggregate HS6 to HS4", n_obs = nrow(raw),
             n_dropped = n_before - nrow(raw),
             trade_value = sum(raw$cusval, na.rm = TRUE),
             detail = sprintf("%d unique HS4 products", uniqueN(raw$good)))
  }
  } # end if/else raw_cache

  # Regional aggregation
  if (cfg$use_regions) {
    rmap <- if (!is.null(cfg$custom_region_map)) cfg$custom_region_map else build_region_map()
    rmap[, cty_code := as.integer(cty_code)]

    imp_merged <- rmap[data.table(cty_code = as.integer(raw$importer)), on = "cty_code"]
    imp_merged[is.na(region), region := "OTHER"]
    raw[, importer := imp_merged$region]

    exp_merged <- rmap[data.table(cty_code = as.integer(raw$exporter)), on = "cty_code"]
    exp_merged[is.na(region), region := "OTHER"]
    raw[, exporter := exp_merged$region]

    n_unmapped_imp <- sum(raw$importer == "OTHER")
    n_unmapped_exp <- sum(raw$exporter == "OTHER")

    n_before <- nrow(raw)
    raw <- raw[, .(cusval = sum(cusval), quantity = sum(quantity)),
               by = .(year, good, importer, exporter)]
    qlog$add("Regional aggregation (Soderbery Table 1)", n_obs = nrow(raw),
             n_dropped = n_before - nrow(raw),
             trade_value = sum(raw$cusval, na.rm = TRUE),
             detail = sprintf("%d importers, %d exporters; %s unmapped imp, %s unmapped exp",
                              uniqueN(raw$importer), uniqueN(raw$exporter),
                              format(n_unmapped_imp, big.mark = ","),
                              format(n_unmapped_exp, big.mark = ",")))
  } else {
    raw[, `:=`(importer = as.character(importer), exporter = as.character(exporter))]
    qlog$add("No regional aggregation (individual countries)", n_obs = nrow(raw),
             trade_value = sum(raw$cusval, na.rm = TRUE),
             detail = sprintf("%d importers, %d exporters",
                              uniqueN(raw$importer), uniqueN(raw$exporter)))
  }

  # Prices, shares, first-differences
  dt <- copy(raw)
  dt[, `:=`(t = year - cfg$minyear + 1L, lp = log(cusval / quantity))]
  dt[, imp_total := sum(cusval), by = .(t, importer, good)]
  dt[, `:=`(s_imp = cusval / imp_total, ls_imp = log(cusval / imp_total))]
  dt[, exp_total := sum(cusval), by = .(t, exporter, good)]
  dt[, `:=`(s_exp = cusval / exp_total, ls_exp = log(cusval / exp_total))]

  setorder(dt, importer, exporter, good, t)
  dt[, `:=`(lp_dif = lp - shift(lp, 1L),
            ls_imp_dif = ls_imp - shift(ls_imp, 1L),
            ls_exp_dif = ls_exp - shift(ls_exp, 1L),
            period_count = .N), by = .(importer, exporter, good)]

  n_before <- nrow(dt)
  dt <- dt[!is.na(lp_dif)]
  qlog$add("First-differencing (drop first obs per panel)",
           n_obs = nrow(dt), n_dropped = n_before - nrow(dt),
           trade_value = sum(dt$cusval, na.rm = TRUE))

  if (!is.na(cfg$uv_outlier_threshold) && cfg$uv_outlier_threshold > 0) {
    thresh <- cfg$uv_outlier_threshold
    n_before <- nrow(dt)
    dt <- dt[abs(lp_dif) < thresh]
    qlog$add(sprintf("Unit value outlier filter |d ln(p)| < %.1f", thresh),
             n_obs = nrow(dt), n_dropped = n_before - nrow(dt),
             trade_value = sum(dt$cusval, na.rm = TRUE),
             detail = sprintf("Threshold = factor of ~%.1f in levels", exp(thresh)))
  }

  qlog$add("Data entering estimation", n_obs = nrow(dt),
           trade_value = sum(dt$cusval, na.rm = TRUE),
           detail = sprintf("%d products, %d importers, %d exporters",
                            uniqueN(dt$good), uniqueN(dt$importer),
                            uniqueN(dt$exporter)))

  cat(sprintf("\nEstimation sample: %s obs, %d products, %d importers, %d exporters\n",
              format(nrow(dt), big.mark = ","),
              uniqueN(dt$good), uniqueN(dt$importer), uniqueN(dt$exporter)))

  cell_stats <- dt[, .(n_exp = uniqueN(exporter), n_per = uniqueN(t)),
                   by = .(importer, good)]
  cat(sprintf("  Cells: %s | Exporters/cell: median=%d, mean=%.1f | Periods/cell: median=%d, mean=%.1f\n\n",
              format(nrow(cell_stats), big.mark = ","),
              median(cell_stats$n_exp), mean(cell_stats$n_exp),
              median(cell_stats$n_per), mean(cell_stats$n_per)))

  list(dt = dt, qlog = qlog)
}


# ===========================================================================
#  PARALLEL ESTIMATION ENGINE
# ===========================================================================

#' @param cfg Config list.
#' @param ncores Number of CPU cores. NULL = detectCores() - 2.
#' @return data.table of estimates.
estimate_all_parallel <- function(cfg, ncores = NULL) {

  if (is.null(ncores)) ncores <- max(1L, detectCores() - 2L)
  ncores <- min(ncores, detectCores())
  cat(sprintf("Estimation: %d cores, OS = %s\n\n", ncores, .Platform$OS.type))

  prep <- prepare_data(cfg)
  dt <- prep$dt; qlog <- prep$qlog

  products <- unique(dt$good)
  n_products <- length(products)
  dt_by_product <- split(dt, by = "good", keep.by = TRUE)

  worker_fns <- c("estimate_product", "estimate_importer_product",
                   "choose_reference", "bw_weight", "optimal_tariff",
                   "assign_regions", "build_region_map",
                   "build_export_moments", "compute_exporter_lookup",
                   "compute_exporter_weights", "cell_failure")

  is_windows <- .Platform$OS.type == "windows"
  t_start <- Sys.time()

  if (ncores == 1L) {
    cat("Running serially (ncores = 1) with incremental checkpoints...\n")

    # --- Checkpoint configuration ---
    checkpoint_every <- 50L  # Save every N products
    checkpoint_file  <- paste0(build_output_prefix(cfg), "_checkpoint.rds")

    # --- Resume from checkpoint if available ---
    results_list <- list()
    start_idx <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_idx <- ckpt$next_idx
      cat(sprintf("  Resuming from checkpoint: %d/%d products already done (%s)\n",
                  start_idx - 1L, n_products, checkpoint_file))
    }

    if (start_idx <= n_products) {
      for (idx in start_idx:n_products) {
        g <- products[idx]
        if (idx %% 10 == 0 || idx == start_idx) {
          el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
          done <- idx - start_idx
          rate <- if (done > 0L) el / done else NA
          eta  <- if (done > 0L) rate * (n_products - idx) else NA
          cat(sprintf("  [%d/%d] %.1f min elapsed%s\n",
                      idx, n_products, el,
                      if (!is.na(eta)) sprintf(", ~%.1f min remaining", eta) else ""))
        }
        res <- estimate_product(g, dt_by_product[[g]], cfg)
        results_list[[idx]] <- res

        # --- Checkpoint save ---
        if (idx %% checkpoint_every == 0L) {
          saveRDS(list(results = results_list, next_idx = idx + 1L),
                  checkpoint_file)
          cat(sprintf("    >> Checkpoint saved at product %d/%d (%s)\n",
                      idx, n_products, checkpoint_file))
        }
      }
    } # end if (start_idx <= n_products)

    # Final checkpoint (in case n_products is not a multiple of checkpoint_every)
    if (n_products %% checkpoint_every != 0L) {
      saveRDS(list(results = results_list, next_idx = n_products + 1L),
              checkpoint_file)
    }

    # Clean up checkpoint after successful completion
    # (kept until post-processing succeeds — deleted at the end of the function)

  } else if (is_windows) {
    cat("Starting socket cluster (Windows)...\n")

    # --- Write product slices to temp files ---
    # Windows socket clusters serialize everything sent to workers.
    # With 90M+ rows, sending data through parLapply causes memory
    # failures. Instead, we write each product's data to a small
    # temp RDS file and have workers read only their own slice.
    tmp_dir <- file.path(tempdir(), "het_products")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

    cat(sprintf("  Writing %d product slices to temp dir...\n", n_products))
    for (g in products) {
      saveRDS(dt_by_product[[g]], file.path(tmp_dir, paste0(g, ".rds")))
    }
    # Free memory — workers will read from disk
    rm(dt_by_product); gc()
    cat("  Product slices written. Starting workers...\n")

    cl <- makeCluster(ncores)
    on.exit({
      tryCatch(stopCluster(cl), error = function(e) NULL)
      unlink(tmp_dir, recursive = TRUE)
    }, add = TRUE)

    clusterExport(cl, varlist = worker_fns, envir = environment())
    clusterExport(cl, varlist = c("cfg", "tmp_dir"), envir = environment())

    # Load het_obj on each worker: try Rcpp, fall back to R
    clusterEvalQ(cl, {
      library(data.table)
      .loaded_rcpp <- FALSE
      tryCatch({
        if (requireNamespace("Rcpp", quietly = TRUE) &&
            file.exists("het_obj_rcpp.cpp")) {
          Rcpp::sourceCpp("het_obj_rcpp.cpp")
          het_obj <- het_obj_rcpp
          .loaded_rcpp <- TRUE
        }
      }, error = function(e) NULL)
      if (!.loaded_rcpp) source("het_obj.R")
    })

    batch_size <- ncores * 2L
    n_batches <- ceiling(n_products / batch_size)
    checkpoint_file <- paste0(build_output_prefix(cfg), "_checkpoint.rds")

    # --- Resume from checkpoint if available ---
    results_list <- list()
    start_batch <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_batch <- ckpt$next_batch
      cat(sprintf("  Resuming from checkpoint: %d/%d batches already done (%s)\n",
                  start_batch - 1L, n_batches, checkpoint_file))
    }

    cat(sprintf("Estimating %d products in %d batches of ~%d...\n\n",
                n_products, n_batches, batch_size))

    for (b in seq(start_batch, n_batches)) {
      if (b > n_batches) break
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      # Send only product names — workers read their own data from disk
      batch_res <- parLapply(cl, batch_products, function(g) {
        dt_g <- readRDS(file.path(tmp_dir, paste0(g, ".rds")))
        estimate_product(g, dt_g, cfg)
      })
      results_list <- c(results_list, batch_res)

      # --- Checkpoint save after every batch ---
      saveRDS(list(results = results_list, next_batch = b + 1L),
              checkpoint_file)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d of %d (%.0f%%) | %.1f min, ~%.1f min left\n",
                  b, n_batches, idx_s, idx_e, n_products,
                  100 * idx_e / n_products, el, el / idx_e * (n_products - idx_e)))
    }

  } else {
    cat("Using forked parallelism (Unix/Mac)...\n")
    # Batch checkpointing for the joint estimator.
    checkpoint_file <- paste0(build_output_prefix(cfg), "_checkpoint.rds")
    batch_size <- max(ncores * 4L, 50L)
    n_batches <- ceiling(n_products / batch_size)

    results_list <- list()
    start_batch <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_batch <- ckpt$next_batch
      cat(sprintf("  Resuming joint from checkpoint: %d/%d batches done\n",
                  start_batch - 1L, n_batches))
    }

    cat(sprintf("  Forked parallel: %d products in %d batches of ~%d\n\n",
                n_products, n_batches, batch_size))

    for (b in seq(start_batch, n_batches)) {
      if (b > n_batches) break
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      batch_res <- mclapply(batch_products, function(g)
        estimate_product(g, dt_by_product[[g]], cfg), mc.cores = ncores)
      results_list <- c(results_list, batch_res)

      saveRDS(list(results = results_list, next_batch = b + 1L),
              checkpoint_file)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d (%.0f%%) | %.1f min\n",
                  b, n_batches, idx_s, idx_e, 100*idx_e/n_products, el))
    }
  }

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

  # --- Combine and truncate ---
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0L) { cat("\nNo estimates produced.\n"); return(NULL) }

  # Extract timing and failures before rbindlist strips attributes
  timing_info <- lapply(results_list, function(r) attr(r, "timing"))
  timing_info <- timing_info[!sapply(timing_info, is.null)]

  failure_info <- unlist(lapply(results_list, function(r) attr(r, "failures")),
                         recursive = FALSE)
  if (is.null(failure_info)) failure_info <- list()

  output <- rbindlist(results_list)
  n_raw <- nrow(output)
  n_succeeded <- length(results_list)

  # -------------------------------------------------------------------
  # OUTPUT COLUMN INTERPRETATION
  #
  # sigma  Elasticity of substitution (positive, > 1). NOT the price
  #        elasticity of demand. Price elasticity = -sigma.
  # gamma  Inverse export supply elasticity (positive). Supply
  #        elasticity = 1/gamma. Pass-through = 1/(1+gamma).
  # opt_tariff  Heterogeneous formula (Eq. 15). Do NOT use 1/gamma.
  # -------------------------------------------------------------------

  # --- Post-estimation trimming ---
  # Symmetric percentile trim applied to both sigma and gamma.
  # Drops the top and bottom tail_trim_pct of each distribution.
  # This is the approach Soderbery uses for gamma (paper Table 2
  # footnote); we extend it to sigma for consistency rather than
  # hardcoding his 99.5th percentile (131.05), which is specific
  # to his Comtrade sample.

  n_trim_total <- 0L
  trim_pct <- cfg$tail_trim_pct

  if (!is.na(trim_pct) && trim_pct > 0) {
    n_b <- nrow(output)
    output <- output[!is.na(sigma) & !is.na(gamma)]

    sig_lo <- quantile(output$sigma, trim_pct, na.rm = TRUE)
    sig_hi <- quantile(output$sigma, 1 - trim_pct, na.rm = TRUE)
    gam_lo <- quantile(output$gamma, trim_pct, na.rm = TRUE)
    gam_hi <- quantile(output$gamma, 1 - trim_pct, na.rm = TRUE)

    output <- output[sigma >= sig_lo & sigma <= sig_hi &
                     gamma >= gam_lo & gamma <= gam_hi]

    n_trim_total <- n_b - nrow(output)

    cat(sprintf("  Trimmed: %s rows (%.1f%% each tail)\n",
                format(n_trim_total, big.mark = ","), trim_pct * 100))
    cat(sprintf("    Sigma kept: [%.2f, %.2f]  Gamma kept: [%.3f, %.3f]\n",
                sig_lo, sig_hi, gam_lo, gam_hi))
  }

  output <- output[, .(importer, exporter, good, sigma, gamma,
                       ref_exporter, opt_tariff, convergence, obj_value)]

  # --- Quality log ---
  qlog$add("Estimation results (pre-trim)", n_obs = n_raw,
           detail = sprintf("%d products succeeded, %d failed, %.1f min",
                            n_succeeded, n_products - n_succeeded, t_elapsed))
  qlog$add(sprintf("Symmetric tail trim (%.1f%% each tail, both sigma and gamma)",
                   trim_pct * 100),
           n_obs = nrow(output), n_dropped = n_trim_total)

  # --- Final summary ---
  cat("\n=============================================\n")
  cat("ESTIMATION COMPLETE\n")
  cat("=============================================\n\n")

  if (length(timing_info) > 0L) {
    ps <- sapply(timing_info, function(x) x$seconds)
    pc <- sapply(timing_info, function(x) x$cells)
    po <- sapply(timing_info, function(x) x$succeeded)
    cat(sprintf("  Time: %.1f min (%.2f hours) | %.1f products/min\n",
                t_elapsed, t_elapsed / 60, length(timing_info) / t_elapsed))
    cat(sprintf("  Per product: median=%.1fs, mean=%.1fs, max=%.1fs\n",
                median(ps), mean(ps), max(ps)))
    cat(sprintf("  Cell success: %d/%d (%.1f%%)\n\n",
                sum(po), sum(pc), 100 * sum(po) / sum(pc)))
  }

  # --- Failure summary ---
  if (length(failure_info) > 0L) {
    fail_dt <- rbindlist(failure_info)
    fail_counts <- fail_dt[, .N, by = reason]
    setorder(fail_counts, -N)
    cat(sprintf("  Cell failures: %d total\n", nrow(fail_dt)))
    for (i in seq_len(nrow(fail_counts))) {
      cat(sprintf("    %-30s %d\n", fail_counts$reason[i], fail_counts$N[i]))
    }
    cat("\n")
  }

  cat(sprintf("  Estimates: %s rows\n", format(nrow(output), big.mark = ",")))
  cat(sprintf("  Sigma:  median=%.3f, mean=%.3f, IQR=[%.3f, %.3f]\n",
              median(output$sigma, na.rm = TRUE), mean(output$sigma, na.rm = TRUE),
              quantile(output$sigma, .25, na.rm = TRUE),
              quantile(output$sigma, .75, na.rm = TRUE)))
  cat(sprintf("  Gamma:  median=%.3f, mean=%.3f, MAD=%.3f\n",
              median(output$gamma, na.rm = TRUE), mean(output$gamma, na.rm = TRUE),
              mad(output$gamma, na.rm = TRUE)))
  cat(sprintf("  Opt tariff: median=%.3f, mean=%.3f\n",
              median(output$opt_tariff, na.rm = TRUE),
              mean(output$opt_tariff, na.rm = TRUE)))

  print_quality_log(qlog)

  # --- Attach metadata for summary generation ---
  attr(output, "run_meta") <- list(
    qlog         = qlog,
    timing_info  = timing_info,
    failure_info = failure_info,
    n_products   = n_products,
    n_succeeded  = n_succeeded,
    n_failed     = n_products - n_succeeded,
    t_elapsed    = t_elapsed,
    ncores       = ncores,
    rcpp_loaded  = .het_obj_rcpp_loaded,
    trim_pct     = trim_pct,
    trim_bounds  = if (exists("sig_lo")) list(
      sig_lo = sig_lo, sig_hi = sig_hi, gam_lo = gam_lo, gam_hi = gam_hi
    ) else NULL,
    n_pre_trim   = n_raw,
    n_trimmed    = n_trim_total,
    timestamp    = Sys.time()
  )

  # --- Clean up checkpoint file after successful completion ---
  checkpoint_file_cleanup <- paste0(build_output_prefix(cfg), "_checkpoint.rds")
  if (file.exists(checkpoint_file_cleanup)) {
    file.remove(checkpoint_file_cleanup)
    cat(sprintf("  Checkpoint file removed: %s\n", checkpoint_file_cleanup))
  }

  output
}


# ===========================================================================
#  ITERATION AND STARTING VALUE HELPERS
# ===========================================================================

#' Update defaults from completed results for iterative refinement.
#'
#' NOTE: currently unused by run_estimation.R (three-stage pipeline uses
#' explicit Stage 1 -> Stage 2a -> Stage 2b handoff with its own default
#' updates). Retained for interactive use / iterative robustness checks.
update_defaults_from_results <- function(cfg, results, pass = 1L) {
  ms <- median(results$sigma, na.rm = TRUE)
  mg <- median(results$gamma, na.rm = TRUE)
  cfg$sigma_start <- ms; cfg$gamma_start <- mg
  cfg$sigma_V_default <- ms; cfg$gamma_V_default <- mg
  cat(sprintf("  Defaults updated from pass %d: sigma=%.3f, gamma=%.3f\n", pass, ms, mg))
  cfg
}


#' Initialize country-level starting values from regional estimates.
#'
#' Creates a lookup table mapping (region, product) -> (sigma, gamma)
#' from regional estimation results. When set on config, the estimator
#' uses per-cell starting values instead of global defaults, which
#' speeds convergence and reduces failures.
#'
#' @param cfg Config list (for country-level estimation).
#' @param regional_results data.table from a regional estimation run.
#' @param custom_map Optional region map. NULL uses build_region_map().
#' @return Config list with regional_starts and regional_starts_rmap set.
init_from_regional <- function(cfg, regional_results, custom_map = NULL) {

  # Build per-(region, product) medians from regional results
  starts <- regional_results[!is.na(sigma) & sigma > 1 & gamma > 0,
                              .(sigma = median(sigma),
                                gamma = median(gamma)),
                              by = .(region = importer, good)]

  cfg$regional_starts <- starts
  cfg$regional_starts_rmap <- if (!is.null(custom_map)) custom_map else build_region_map()

  cat(sprintf("  Regional starting values loaded: %d (region x product) cells\n",
              nrow(starts)))
  cat(sprintf("  Coverage: sigma range [%.2f, %.2f], gamma range [%.2f, %.2f]\n",
              min(starts$sigma), max(starts$sigma),
              min(starts$gamma), max(starts$gamma)))
  cfg
}


# ===========================================================================
#  LAMBDA CALIBRATION DIAGNOSTIC
#
#  Computes the four criteria from the README for choosing shrinkage
#  lambda empirically:
#    1. Within-pair gamma MAD — should approach Soderbery's 0.125
#    2. Cross-cell variance decomposition R-squared — should approach 0.72
#    3. Share of gamma estimates on the plateau — should be near 0
#    4. Distinctness of Tier 1 vs Tier 3 distributions (KS-style gap)
#
#  Designed to be called on the output of estimate_all_fixed_sigma()
#  for a single lambda value. Typical use:
#
#    for (lam in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
#      cfg_test$shrinkage_lambda <- lam
#      res <- estimate_all_fixed_sigma(cfg_test, prepared_dt = dt_test)
#      diag <- lambda_calibration_diagnostic(res, lambda = lam)
#      print(diag)
#    }
#
#  Targets listed in the return value are the Soderbery benchmarks.
#  Interpretation:
#    - Very low lambda (0.01): weak regularization, high plateau share
#    - Very high lambda (0.5): within-pair MAD collapses (overshrinkage)
#    - Sweet spot: plateau share ~0, within-pair MAD near 0.125,
#      Tier 1 / Tier 3 distributions clearly distinct
# ===========================================================================

#' Compute lambda calibration diagnostic for a single fixed-sigma run.
#'
#' @param results Output of estimate_all_fixed_sigma().
#' @param lambda Scalar lambda value used for the run (for labeling).
#' @param plateau_gamma_cutoff gamma value above which estimates are
#'   considered on the plateau. Default 20 (gamma/(1+gamma) > 0.95).
#' @return data.table with one row of diagnostics. Columns:
#'   lambda, n_estimates, n_pairs, within_pair_mad,
#'   r_squared, plateau_share, tier1_gamma_median, tier3_gamma_median,
#'   tier_distinctness
lambda_calibration_diagnostic <- function(results, lambda = NA_real_,
                                           plateau_gamma_cutoff = 20) {
  dt <- results[!is.na(gamma) & gamma > 0]
  if (nrow(dt) == 0L) return(NULL)

  # --- 1. Within-pair gamma MAD ---
  pair_mad <- dt[, .(
    n_exp = .N,
    mad_g = if (.N >= 2L) mad(gamma, constant = 1) else NA_real_
  ), by = .(importer, good)]
  pair_mad <- pair_mad[n_exp >= 2L]
  within_pair_mad <- median(pair_mad$mad_g, na.rm = TRUE)
  n_pairs <- nrow(pair_mad)

  # --- 2. Cross-cell variance decomposition R² ---
  vd <- tryCatch(variance_decomposition(dt), error = function(e) NULL)
  r_sq <- if (!is.null(vd)) vd$r_squared else NA_real_

  # --- 3. Plateau share ---
  plateau_share <- mean(dt$gamma > plateau_gamma_cutoff, na.rm = TRUE)

  # --- 4. Tier distinctness ---
  tier1_med <- NA_real_; tier3_med <- NA_real_; distinctness <- NA_real_
  if ("tier" %in% names(dt)) {
    t1 <- dt[tier == 1L, gamma]
    t3 <- dt[tier == 3L, gamma]
    if (length(t1) > 10L && length(t3) > 10L) {
      tier1_med <- median(t1)
      tier3_med <- median(t3)
      # "Distinctness" = absolute log-median gap. Larger = more distinct.
      # When shrinkage is too aggressive, Tier 1 distributions collapse
      # toward the Tier 3 prior value and this approaches 0.
      distinctness <- abs(log(tier1_med) - log(tier3_med))
    }
  }

  data.table(
    lambda              = lambda,
    n_estimates         = nrow(dt),
    n_pairs             = n_pairs,
    within_pair_mad     = within_pair_mad,
    mad_target          = 0.125,
    r_squared           = r_sq,
    r_squared_target    = 0.72,
    plateau_share       = plateau_share,
    tier1_gamma_median  = tier1_med,
    tier3_gamma_median  = tier3_med,
    tier_distinctness   = distinctness
  )
}


#' Print a lambda calibration diagnostic in a human-readable format.
#'
#' @param diag data.table from lambda_calibration_diagnostic().
print_lambda_diagnostic <- function(diag) {
  if (is.null(diag) || nrow(diag) == 0L) {
    cat("  (no diagnostic data)\n")
    return(invisible())
  }
  cat(sprintf("  Lambda = %g\n", diag$lambda[1]))
  cat(sprintf("    Estimates:             %s (%s pairs)\n",
              format(diag$n_estimates[1], big.mark = ","),
              format(diag$n_pairs[1], big.mark = ",")))
  cat(sprintf("    Within-pair gamma MAD: %.3f  (target: %.3f)\n",
              diag$within_pair_mad[1], diag$mad_target[1]))
  cat(sprintf("    R-squared (imp x good FE): %.3f  (target: %.3f)\n",
              diag$r_squared[1], diag$r_squared_target[1]))
  cat(sprintf("    Plateau share (g>20):  %.4f  (target: ~0)\n",
              diag$plateau_share[1]))
  if (!is.na(diag$tier_distinctness[1])) {
    cat(sprintf("    Tier distinctness:     %.3f (log gap |log T1 - log T3|)\n",
                diag$tier_distinctness[1]))
    cat(sprintf("    Tier 1 gamma median:   %.3f\n",
                diag$tier1_gamma_median[1]))
    cat(sprintf("    Tier 3 gamma median:   %.3f\n",
                diag$tier3_gamma_median[1]))
  }
  invisible()
}


#' Run a lambda sweep and return combined diagnostics.
#'
#' Intended as a one-shot driver for calibrating shrinkage lambda on a
#' subset of products. Takes a prepared data.table, runs Stage 2 at each
#' lambda, and returns a single combined diagnostic table.
#'
#' @param cfg_base Base config list (sigma_lookup, shrinkage_priors, etc.
#'   all set). shrinkage_lambda is overridden per run.
#' @param prepared_dt Prepared data.table (as from prepare_data()$dt).
#' @param lambda_grid Numeric vector of lambdas to sweep.
#' @param ncores Number of cores to use.
#' @return data.table of diagnostics, one row per lambda.
lambda_calibration_sweep <- function(cfg_base, prepared_dt,
                                      lambda_grid = c(0.01, 0.05, 0.1, 0.2, 0.5),
                                      ncores = NULL) {
  diags <- list()
  for (lam in lambda_grid) {
    cat(sprintf("\n=== Lambda sweep: lambda = %g ===\n", lam))
    cfg_lam <- cfg_base
    cfg_lam$shrinkage_lambda <- lam
    res <- estimate_all_fixed_sigma(cfg_lam, ncores = ncores,
                                     prepared_dt = prepared_dt)
    diag <- lambda_calibration_diagnostic(res, lambda = lam)
    print_lambda_diagnostic(diag)
    diags[[length(diags) + 1L]] <- diag
  }
  rbindlist(diags)
}


# ===========================================================================
#  ESTIMATION SUMMARY
#
#  Builds a comprehensive summary from estimation results, modeled on
#  the diagnostics reported in Soderbery (2018):
#    - Table 2-style per-importer distributional statistics
#    - Variance decomposition of gamma (importer×product FE R²)
#    - Within importer-product heterogeneity statistics
#    - Trade value coverage
#    - Performance and convergence diagnostics
#    - Config provenance
#
#  Output: both a machine-readable RDS list and a formatted text report.
# ===========================================================================


#' Build a per-importer summary table in the style of Soderbery Table 2.
#'
#' For each importer, reports: observation count, mean/median/MAD of
#' sigma and gamma. The "World" row gives pooled statistics.
#'
#' @param results data.table of estimation results.
#' @return data.table with one row per importer plus a "World" row.
build_table2 <- function(results) {
  dt <- results[!is.na(sigma) & !is.na(gamma)]
  by_imp <- dt[, .(
    obs          = .N,
    sigma_mean   = mean(sigma),
    sigma_median = median(sigma),
    sigma_mad    = mad(sigma, constant = 1),
    gamma_mean   = mean(gamma),
    gamma_median = median(gamma),
    gamma_mad    = mad(gamma, constant = 1)
  ), by = importer]
  setorder(by_imp, importer)

  world <- dt[, .(
    importer     = "World",
    obs          = .N,
    sigma_mean   = mean(sigma),
    sigma_median = median(sigma),
    sigma_mad    = mad(sigma, constant = 1),
    gamma_mean   = mean(gamma),
    gamma_median = median(gamma),
    gamma_mad    = mad(gamma, constant = 1)
  )]

  rbindlist(list(by_imp, world))
}


#' Compute within importer-product heterogeneity statistics.
#'
#' For each (importer, product), computes the within-pair median,
#' SD, and MAD of gamma. Then reports the median of these statistics
#' across all importer-product pairs, plus 25th/75th percentiles.
#' (Corresponds to Soderbery p. 52 discussion.)
#'
#' @param results data.table of estimation results.
#' @return Named list of within-pair statistics.
within_pair_stats <- function(results) {
  dt <- results[!is.na(gamma) & gamma > 0]
  # Need at least 2 exporters per (importer, good) to have within-pair variation
  pair_stats <- dt[, .(
    n_exporters  = .N,
    within_med   = median(gamma),
    within_sd    = if (.N > 1L) sd(gamma) else NA_real_,
    within_mad   = mad(gamma, constant = 1)
  ), by = .(importer, good)]

  pair_stats <- pair_stats[n_exporters >= 2L]

  if (nrow(pair_stats) == 0L) return(NULL)

  list(
    n_pairs                = nrow(pair_stats),
    median_of_within_med   = median(pair_stats$within_med, na.rm = TRUE),
    median_of_within_sd    = median(pair_stats$within_sd, na.rm = TRUE),
    median_of_within_mad   = median(pair_stats$within_mad, na.rm = TRUE),
    q25_within_med         = quantile(pair_stats$within_med, 0.25, na.rm = TRUE),
    q75_within_med         = quantile(pair_stats$within_med, 0.75, na.rm = TRUE),
    q25_within_sd          = quantile(pair_stats$within_sd, 0.25, na.rm = TRUE),
    q75_within_sd          = quantile(pair_stats$within_sd, 0.75, na.rm = TRUE),
    q25_within_mad         = quantile(pair_stats$within_mad, 0.25, na.rm = TRUE),
    q75_within_mad         = quantile(pair_stats$within_mad, 0.75, na.rm = TRUE)
  )
}


#' Compute variance decomposition of gamma.
#'
#' Regresses log(gamma) on importer×product fixed effects and reports R².
#' This measures how much of the variation in export supply elasticities
#' is explained by the import market vs. exporter heterogeneity within
#' markets. Soderbery reports 72% (p. 52).
#'
#' @param results data.table of estimation results.
#' @return Named list with r_squared and n_obs, or NULL if insufficient data.
variance_decomposition <- function(results) {
  dt <- results[!is.na(gamma) & gamma > 0]
  dt[, log_gamma := log(gamma)]
  dt[, imp_good := paste(importer, good, sep = "_")]

  # Need enough variation — at least some groups with > 1 obs
  grp <- dt[, .N, by = imp_good]
  if (sum(grp$N > 1L) < 10L) return(NULL)

  # Use within-group SS / total SS for R² (equivalent to FE regression)
  grand_mean <- mean(dt$log_gamma, na.rm = TRUE)
  ss_total   <- sum((dt$log_gamma - grand_mean)^2, na.rm = TRUE)

  group_means <- dt[, .(gm = mean(log_gamma, na.rm = TRUE)), by = imp_good]
  dt <- group_means[dt, on = "imp_good"]
  ss_within  <- sum((dt$log_gamma - dt$gm)^2, na.rm = TRUE)

  r_sq <- 1 - ss_within / ss_total

  list(
    r_squared = r_sq,
    n_obs     = nrow(dt),
    n_groups  = nrow(group_means)
  )
}


#' Build the full estimation summary.
#'
#' @param results data.table of final (Step 2) results.
#' @param cfg Config list used for the estimation.
#' @param step1_results Optional data.table of Step 1 results for comparison.
#' @param scope Character: "regional" or "country".
#' @return Named list containing all summary components.
build_summary <- function(results, cfg, step1_results = NULL, scope = NULL) {

  if (is.null(scope)) {
    scope <- if (isTRUE(cfg$use_regions)) "regional" else "country"
  }

  meta <- attr(results, "run_meta")
  dt <- results[!is.na(sigma) & !is.na(gamma)]

  # --- Config provenance ---
  provenance <- list(
    baci_source       = parse_baci_source(cfg$filepath),
    filepath          = cfg$filepath,
    scope             = scope,
    agg_level         = cfg$agg_level,
    minyear           = cfg$minyear,
    maxyear           = cfg$maxyear,
    min_exporters     = cfg$min_exporters,
    min_destinations  = cfg$min_destinations,
    min_periods       = cfg$min_periods,
    uv_outlier_thresh = cfg$uv_outlier_threshold,
    tail_trim_pct     = cfg$tail_trim_pct,
    sigma_V_default   = cfg$sigma_V_default,
    gamma_V_default   = cfg$gamma_V_default,
    sigma_start       = cfg$sigma_start,
    gamma_start       = cfg$gamma_start,
    timestamp         = if (!is.null(meta)) meta$timestamp else Sys.time()
  )

  # --- Performance ---
  perf <- list(ncores = NULL, rcpp = NULL, elapsed_min = NULL,
               products_per_min = NULL, per_product_median_s = NULL,
               per_product_mean_s = NULL, per_product_max_s = NULL,
               cells_attempted = NULL, cells_succeeded = NULL,
               cell_success_rate = NULL)
  if (!is.null(meta)) {
    perf$ncores            <- meta$ncores
    perf$rcpp              <- meta$rcpp_loaded
    perf$elapsed_min       <- meta$t_elapsed
    perf$elapsed_hours     <- meta$t_elapsed / 60
    perf$n_products        <- meta$n_products
    perf$n_succeeded       <- meta$n_succeeded
    perf$n_failed          <- meta$n_failed
    perf$products_per_min  <- meta$n_succeeded / max(meta$t_elapsed, 0.01)
    if (length(meta$timing_info) > 0L) {
      ps <- sapply(meta$timing_info, function(x) x$seconds)
      pc <- sapply(meta$timing_info, function(x) x$cells)
      po <- sapply(meta$timing_info, function(x) x$succeeded)
      perf$per_product_median_s <- median(ps)
      perf$per_product_mean_s   <- mean(ps)
      perf$per_product_max_s    <- max(ps)
      perf$cells_attempted      <- sum(pc)
      perf$cells_succeeded      <- sum(po)
      perf$cell_success_rate    <- 100 * sum(po) / max(sum(pc), 1)
    }
  }

  # --- Data quality (from qlog) ---
  quality <- NULL
  if (!is.null(meta) && !is.null(meta$qlog)) {
    ql <- meta$qlog
    quality <- data.table(
      stage       = sapply(ql$steps, function(s) s$stage),
      n_obs       = sapply(ql$steps, function(s) s$n_obs),
      n_dropped   = sapply(ql$steps, function(s) s$n_dropped),
      trade_value = sapply(ql$steps, function(s) s$trade_value),
      detail      = sapply(ql$steps, function(s) s$detail)
    )
  }

  # --- Global distributional statistics ---
  global <- list(
    n_estimates    = nrow(dt),
    n_products     = uniqueN(dt$good),
    n_importers    = uniqueN(dt$importer),
    n_exporters    = uniqueN(dt$exporter),
    n_sigma        = nrow(unique(dt[, .(importer, good)])),
    n_gamma        = nrow(dt),
    sigma_median   = median(dt$sigma),
    sigma_mean     = mean(dt$sigma),
    sigma_sd       = sd(dt$sigma),
    sigma_q25      = as.numeric(quantile(dt$sigma, 0.25)),
    sigma_q75      = as.numeric(quantile(dt$sigma, 0.75)),
    gamma_median   = median(dt$gamma),
    gamma_mean     = mean(dt$gamma),
    gamma_sd       = sd(dt$gamma),
    gamma_mad      = mad(dt$gamma, constant = 1),
    gamma_q25      = as.numeric(quantile(dt$gamma, 0.25)),
    gamma_q75      = as.numeric(quantile(dt$gamma, 0.75)),
    opt_tariff_median = median(dt$opt_tariff, na.rm = TRUE),
    opt_tariff_mean   = mean(dt$opt_tariff, na.rm = TRUE)
  )
  # Trim bounds
  if (!is.null(meta) && !is.null(meta$trim_bounds)) {
    global$trim <- meta$trim_bounds
    global$n_pre_trim <- meta$n_pre_trim
    global$n_trimmed  <- meta$n_trimmed
  }

  # --- Table 2 (per-importer) ---
  table2 <- build_table2(results)

  # --- Within importer-product heterogeneity ---
  within_het <- within_pair_stats(results)

  # --- Variance decomposition ---
  var_decomp <- variance_decomposition(results)

  # --- Two-step comparison ---
  step_comparison <- NULL
  if (!is.null(step1_results)) {
    s1 <- step1_results[!is.na(sigma) & !is.na(gamma)]
    med_s1_sig <- median(s1$sigma, na.rm = TRUE)
    med_s1_gam <- median(s1$gamma, na.rm = TRUE)
    med_s2_sig <- global$sigma_median
    med_s2_gam <- global$gamma_median

    merged <- merge(
      s1[, .(importer, exporter, good, sigma_1 = sigma, gamma_1 = gamma)],
      dt[, .(importer, exporter, good, sigma_2 = sigma, gamma_2 = gamma)],
      by = c("importer", "exporter", "good"))

    step_comparison <- list(
      step1_sigma_median = med_s1_sig,
      step1_gamma_median = med_s1_gam,
      step2_sigma_median = med_s2_sig,
      step2_gamma_median = med_s2_gam,
      sigma_shift        = med_s2_sig - med_s1_sig,
      gamma_shift        = med_s2_gam - med_s1_gam,
      sigma_shift_pct    = 100 * (med_s2_sig - med_s1_sig) / med_s1_sig,
      gamma_shift_pct    = 100 * (med_s2_gam - med_s1_gam) / med_s1_gam,
      n_matched          = nrow(merged)
    )
    if (nrow(merged) > 10L) {
      step_comparison$cor_sigma <- cor(merged$sigma_1, merged$sigma_2, use = "complete.obs")
      step_comparison$cor_gamma <- cor(merged$gamma_1, merged$gamma_2, use = "complete.obs")
    }
  }

  # --- Failure diagnostics ---
  fail_summary <- NULL
  if (!is.null(meta) && length(meta$failure_info) > 0L) {
    fail_dt <- rbindlist(meta$failure_info)
    fail_summary <- fail_dt[, .N, by = reason]
    setorder(fail_summary, -N)
    setnames(fail_summary, c("reason", "count"))
  }

  list(
    provenance      = provenance,
    performance     = perf,
    quality         = quality,
    global          = global,
    table2          = table2,
    within_het      = within_het,
    var_decomp      = var_decomp,
    step_comparison = step_comparison,
    failures        = fail_summary
  )
}


#' Write summary to a formatted text file.
#'
#' @param summary List from build_summary().
#' @param filepath Path to write the text file.
write_summary_text <- function(summary, filepath) {

  lines <- character(0)
  a <- function(...) lines <<- c(lines, sprintf(...))
  rule <- function() a(paste(rep("=", 72), collapse = ""))
  dash <- function() a(paste(rep("-", 72), collapse = ""))

  prov <- summary$provenance
  perf <- summary$performance
  glob <- summary$global
  tb2  <- summary$table2
  wh   <- summary$within_het
  vd   <- summary$var_decomp
  sc   <- summary$step_comparison
  fl   <- summary$failures

  rule()
  a("  ESTIMATION SUMMARY REPORT")
  a("  Soderbery (2018) Heterogeneous Elasticity Estimator")
  rule()
  a("")

  # --- Provenance ---
  a("CONFIGURATION")
  dash()
  a("  BACI source:          %s", prov$baci_source)
  a("  Data path:            %s", prov$filepath)
  a("  Scope:                %s", prov$scope)
  a("  Aggregation:          %s", prov$agg_level)
  a("  Min year:             %d", prov$minyear)
  a("  Max year:             %s", if (is.null(prov$maxyear) || is.na(prov$maxyear))
                                    "all available" else as.character(prov$maxyear))
  a("  Structural defaults:  sigma_V=%.3f, gamma_V=%.3f", prov$sigma_V_default, prov$gamma_V_default)
  a("  Starting values:      sigma=%.3f, gamma=%.3f", prov$sigma_start, prov$gamma_start)
  a("  Min exporters:        %d", prov$min_exporters)
  a("  Min destinations:     %d", prov$min_destinations)
  a("  Min periods:          %d", prov$min_periods)
  a("  UV outlier threshold: %s", if (is.na(prov$uv_outlier_thresh)) "disabled"
                                  else sprintf("%.1f", prov$uv_outlier_thresh))
  a("  Tail trim:            %.1f%% each tail", prov$tail_trim_pct * 100)
  a("  Timestamp:            %s", format(prov$timestamp, "%Y-%m-%d %H:%M:%S"))
  a("")

  # --- Performance ---
  if (!is.null(perf$ncores)) {
    a("PERFORMANCE")
    dash()
    a("  Cores:                %d", perf$ncores)
    a("  Objective function:   %s", if (isTRUE(perf$rcpp)) "Rcpp (C++)" else "pure R")
    a("  Wall clock time:      %.1f min (%.2f hours)", perf$elapsed_min, perf$elapsed_hours)
    a("  Products:             %d total, %d succeeded, %d failed",
      perf$n_products, perf$n_succeeded, perf$n_failed)
    a("  Throughput:           %.1f products/min", perf$products_per_min)
    if (!is.null(perf$per_product_median_s)) {
      a("  Per product (sec):    median=%.1f, mean=%.1f, max=%.1f",
        perf$per_product_median_s, perf$per_product_mean_s, perf$per_product_max_s)
    }
    if (!is.null(perf$cells_attempted)) {
      a("  Cell convergence:     %d / %d (%.1f%%)",
        perf$cells_succeeded, perf$cells_attempted, perf$cell_success_rate)
    }
    a("")
  }

  # --- Data quality ---
  if (!is.null(summary$quality)) {
    ql <- summary$quality
    a("DATA PIPELINE")
    dash()
    has_tv <- any(!is.na(ql$trade_value))
    if (has_tv) {
      a("  %-38s %12s %12s %14s", "Stage", "Obs", "Dropped", "Trade Val ($B)")
      a("  %-38s %12s %12s %14s",
        paste(rep("-", 38), collapse = ""),
        paste(rep("-", 12), collapse = ""),
        paste(rep("-", 12), collapse = ""),
        paste(rep("-", 14), collapse = ""))
    } else {
      a("  %-40s %12s %12s", "Stage", "Obs", "Dropped")
      a("  %-40s %12s %12s",
        paste(rep("-", 40), collapse = ""),
        paste(rep("-", 12), collapse = ""),
        paste(rep("-", 12), collapse = ""))
    }
    for (i in seq_len(nrow(ql))) {
      n_str <- format(ql$n_obs[i], big.mark = ",")
      d_str <- if (is.na(ql$n_dropped[i])) "" else format(ql$n_dropped[i], big.mark = ",")
      if (has_tv) {
        tv_str <- if (is.na(ql$trade_value[i])) "" else sprintf("%.1f", ql$trade_value[i] / 1e6)
        a("  %-38s %12s %12s %14s", ql$stage[i], n_str, d_str, tv_str)
      } else {
        a("  %-40s %12s %12s", ql$stage[i], n_str, d_str)
      }
      if (nchar(ql$detail[i]) > 0) a("    %s", ql$detail[i])
    }
    # Trade value retention
    if (has_tv) {
      tv_first <- ql$trade_value[which(!is.na(ql$trade_value))[1]]
      tv_last  <- ql$trade_value[max(which(!is.na(ql$trade_value)))]
      if (!is.na(tv_first) && !is.na(tv_last) && tv_first > 0) {
        a("")
        a("  Trade value retention: $%.1fB / $%.1fB (%.1f%%)",
          tv_last / 1e6, tv_first / 1e6, 100 * tv_last / tv_first)
      }
    }
    a("")
  }

  # --- Global distributional statistics ---
  a("GLOBAL DISTRIBUTIONAL STATISTICS")
  dash()
  a("  Total estimate rows:  %s", format(glob$n_estimates, big.mark = ","))
  a("  Products:             %d", glob$n_products)
  a("  Importers:            %d", glob$n_importers)
  a("  Exporters:            %d", glob$n_exporters)
  a("  Unique sigma (imp x good): %s", format(glob$n_sigma, big.mark = ","))
  a("  Unique gamma (imp x exp x good): %s", format(glob$n_gamma, big.mark = ","))
  a("")
  a("  Sigma (elasticity of substitution):")
  a("    Median=%.3f  Mean=%.3f  SD=%.3f  IQR=[%.3f, %.3f]",
    glob$sigma_median, glob$sigma_mean, glob$sigma_sd, glob$sigma_q25, glob$sigma_q75)
  a("")
  a("  Gamma (inverse export supply elasticity):")
  a("    Median=%.3f  Mean=%.3f  SD=%.3f  MAD=%.3f  IQR=[%.3f, %.3f]",
    glob$gamma_median, glob$gamma_mean, glob$gamma_sd, glob$gamma_mad,
    glob$gamma_q25, glob$gamma_q75)
  a("")
  a("  Optimal tariff (Proposition 1):")
  a("    Median=%.3f  Mean=%.3f", glob$opt_tariff_median, glob$opt_tariff_mean)
  if (!is.null(glob$trim)) {
    a("")
    a("  Trim bounds applied:")
    a("    Sigma: [%.3f, %.3f]  Gamma: [%.4f, %.4f]",
      glob$trim$sig_lo, glob$trim$sig_hi, glob$trim$gam_lo, glob$trim$gam_hi)
    a("    Pre-trim: %s  Trimmed: %s",
      format(glob$n_pre_trim, big.mark = ","), format(glob$n_trimmed, big.mark = ","))
  }
  a("")

  # --- Table 2 (per-importer) ---
  a("PER-IMPORTER SUMMARY (cf. Soderbery Table 2)")
  dash()
  a("  %-20s %7s %8s %8s %8s %8s %8s %8s",
    "Importer", "Obs", "sig_mn", "sig_md", "sig_MAD", "gam_mn", "gam_md", "gam_MAD")
  a("  %-20s %7s %8s %8s %8s %8s %8s %8s",
    paste(rep("-", 20), collapse = ""),
    paste(rep("-", 7), collapse = ""),
    paste(rep("-", 8), collapse = ""),
    paste(rep("-", 8), collapse = ""),
    paste(rep("-", 8), collapse = ""),
    paste(rep("-", 8), collapse = ""),
    paste(rep("-", 8), collapse = ""),
    paste(rep("-", 8), collapse = ""))
  for (i in seq_len(nrow(tb2))) {
    r <- tb2[i]
    a("  %-20s %7s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",
      r$importer, format(r$obs, big.mark = ","),
      r$sigma_mean, r$sigma_median, r$sigma_mad,
      r$gamma_mean, r$gamma_median, r$gamma_mad)
  }
  a("")

  # --- Within importer-product heterogeneity ---
  if (!is.null(wh)) {
    a("WITHIN IMPORTER-PRODUCT HETEROGENEITY (cf. Soderbery p. 52)")
    dash()
    a("  Importer-product pairs with 2+ exporters: %s", format(wh$n_pairs, big.mark = ","))
    a("")
    a("  Median gamma within pair:  %.3f  [Q25=%.3f, Q75=%.3f]",
      wh$median_of_within_med, wh$q25_within_med, wh$q75_within_med)
    a("  SD of gamma within pair:   %.3f  [Q25=%.3f, Q75=%.3f]",
      wh$median_of_within_sd, wh$q25_within_sd, wh$q75_within_sd)
    a("  MAD of gamma within pair:  %.3f  [Q25=%.3f, Q75=%.3f]",
      wh$median_of_within_mad, wh$q25_within_mad, wh$q75_within_mad)
    a("")
    a("  Interpretation: For the median importer-product, exporters have")
    a("  a median inverse supply elasticity of %.3f with within-pair", wh$median_of_within_med)
    a("  SD of %.3f. Soderbery reports: median=0.712, SD=0.625, MAD=0.125.", wh$median_of_within_sd)
    a("")
  }

  # --- Variance decomposition ---
  if (!is.null(vd)) {
    a("VARIANCE DECOMPOSITION OF GAMMA (cf. Soderbery p. 52)")
    dash()
    a("  R-squared (importer x product FEs on log gamma): %.3f", vd$r_squared)
    a("  Observations: %s  Groups: %s",
      format(vd$n_obs, big.mark = ","), format(vd$n_groups, big.mark = ","))
    a("")
    a("  Interpretation: %.0f%% of variation in log(gamma) is explained by",
      vd$r_squared * 100)
    a("  importer-product fixed effects. The remaining %.0f%% comes from",
      (1 - vd$r_squared) * 100)
    a("  exporter heterogeneity within importer-product pairs.")
    a("  Soderbery reports 72%% explained by importer-product FEs.")
    a("")
  }

  # --- Cell failure diagnostics ---
  if (!is.null(fl) && nrow(fl) > 0L) {
    a("CELL FAILURE DIAGNOSTICS")
    dash()
    total_fails <- sum(fl$count)
    a("  Total failed cells: %s", format(total_fails, big.mark = ","))
    a("")
    a("  %-35s %8s", "Reason", "Count")
    a("  %-35s %8s", paste(rep("-", 35), collapse = ""),
      paste(rep("-", 8), collapse = ""))
    for (i in seq_len(nrow(fl))) {
      a("  %-35s %8s", fl$reason[i], format(fl$count[i], big.mark = ","))
    }
    a("")
  }

  # --- Two-step comparison ---
  if (!is.null(sc)) {
    a("TWO-STEP SENSITIVITY COMPARISON")
    dash()
    a("  Step 1 medians (Comtrade defaults): sigma=%.3f  gamma=%.3f",
      sc$step1_sigma_median, sc$step1_gamma_median)
    a("  Step 2 medians (BACI defaults):     sigma=%.3f  gamma=%.3f",
      sc$step2_sigma_median, sc$step2_gamma_median)
    a("  Shift: sigma %+.3f (%.1f%%)  gamma %+.3f (%.1f%%)",
      sc$sigma_shift, sc$sigma_shift_pct, sc$gamma_shift, sc$gamma_shift_pct)
    if (!is.null(sc$cor_sigma)) {
      a("  Cross-step correlation: sigma=%.4f  gamma=%.4f", sc$cor_sigma, sc$cor_gamma)
    }
    a("  Matched cells: %s", format(sc$n_matched, big.mark = ","))
    a("")
  }

  rule()
  a("  END OF REPORT")
  rule()

  writeLines(lines, filepath)
  cat(sprintf("  Summary text: %s\n", filepath))
}


#' Write full summary (RDS + text).
#'
#' @param results Final (Step 2) results data.table.
#' @param cfg Config list.
#' @param out_prefix Output file prefix from build_output_prefix().
#' @param step1_results Optional Step 1 results for comparison.
#' @param scope Character: "regional" or "country".
write_estimation_summary <- function(results, cfg, out_prefix,
                                     step1_results = NULL, scope = NULL) {
  cat("\nBuilding estimation summary...\n")
  summary <- build_summary(results, cfg, step1_results, scope)

  rds_path  <- paste0(out_prefix, "_summary.rds")
  text_path <- paste0(out_prefix, "_summary.txt")

  saveRDS(summary, rds_path)
  cat(sprintf("  Summary RDS:  %s\n", rds_path))

  write_summary_text(summary, text_path)
}
#
#  Generates a file prefix from config parameters so output files
#  reflect the data source and estimation settings.
#
#  Example outputs:
#    baci_hs92_v202601_elast_regional_hs4
#    baci_hs07_v202601_elast_country_hs6
#
#  The prefix is built from:
#    - baci_source:  e.g., "baci_hs92_v202601" (parsed from filepath)
#    - scope:        "regional" or "country" (from use_regions)
#    - agg_level:    "hs4" or "hs6" (from config)
# ===========================================================================

#' Parse BACI source identifier from the filepath.
#'
#' Extracts the HS revision and version from the BACI directory or file name.
#' Falls back to "baci" if the pattern is not recognized.
#'
#' Examples:
#'   "BACI_HS92_V202601/"       -> "baci_hs92_v202601"
#'   "BACI_HS07_V202601/"       -> "baci_hs07_v202601"
#'   "data/BACI_HS17_V202501/"  -> "baci_hs17_v202501"
#'   "my_trade_data.csv"        -> "baci"
#'
#' @param filepath The BACI data path from config.
#' @return Character string identifying the BACI source.
parse_baci_source <- function(filepath) {
  # Match patterns like BACI_HS92_V202601 anywhere in the path
  m <- regmatches(filepath,
                  regexpr("BACI_HS\\d{2}_V\\d{6}", filepath, ignore.case = TRUE))
  if (length(m) == 1L && nchar(m) > 0L) {
    return(tolower(m))  # e.g., "baci_hs92_v202601"
  }
  "baci"
}


#' Build output file prefix from config.
#'
#' @param cfg Config list with filepath, use_regions, and agg_level.
#' @param scope Character: "regional" or "country". If NULL, inferred
#'   from cfg$use_regions.
#' @return Character string like "baci_hs92_v202601_elast_regional_hs4".
build_output_prefix <- function(cfg, scope = NULL) {
  source_id <- parse_baci_source(cfg$filepath)
  if (is.null(scope)) {
    scope <- if (isTRUE(cfg$use_regions)) "regional" else "country"
  }
  agg <- cfg$agg_level
  paste(source_id, "elast", scope, agg, sep = "_")
}




# ===========================================================================
#  FIXED-SIGMA OBJECTIVE FUNCTION LOADING
# ===========================================================================

.het_obj_fs_rcpp_loaded <- FALSE

tryCatch({
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    cpp_file <- file.path(getwd(), "het_obj_fixed_sigma_rcpp.cpp")
    if (file.exists(cpp_file)) {
      Rcpp::sourceCpp(cpp_file)
      het_obj_fixed_sigma <- het_obj_fixed_sigma_rcpp
      .het_obj_fs_rcpp_loaded <- TRUE
      cat("  Fixed-sigma objective: Rcpp (compiled C++)\n")
    }
  }
}, error = function(e) {
  cat(sprintf("  Fixed-sigma Rcpp compilation failed: %s\n", conditionMessage(e)))
})

if (!.het_obj_fs_rcpp_loaded) {
  het_obj_fixed_sigma <- function(d, sigma, imp_Y, imp_X, exp_Y, exp_X,
                                  exp_jmap, exp_sig_V, exp_gam_V,
                                  wt_imp, wt_exp,
                                  ln_gamma_prior, shrinkage_lambda) {
    d_full <- c(sigma, d)
    ssr <- het_obj(d_full, imp_Y, imp_X, exp_Y, exp_X,
                   exp_jmap, exp_sig_V, exp_gam_V, wt_imp, wt_exp)
    if (ssr >= 1e12) return(1e12)
    if (shrinkage_lambda > 0 && !is.na(ln_gamma_prior)) {
      gam_vals <- d[d > 1e-5]
      if (length(gam_vals) > 0L) {
        ssr <- ssr + shrinkage_lambda * sum((log(gam_vals) - ln_gamma_prior)^2)
      }
    }
    ssr
  }
  cat("  Fixed-sigma objective: pure R (wrapper)\n")
}


# ===========================================================================
#  FEENSTRA (1994) SIGMA ESTIMATION
#
#  Import-side-only, 2-parameter (sigma, gamma_common) objective.
#  Under homogeneity gamma_j = gamma_k, Eq (10) simplifies.
# ===========================================================================

feenstra_sigma_obj <- function(d, imp_Y, imp_X, wt_imp) {
  sig <- d[1]
  gam <- d[2]
  if (sig <= 1 || gam <= 0) return(1e12)
  sm1    <- sig - 1
  g_frac <- gam / (1 + gam)
  pred <- (g_frac / sm1)   * imp_X[, 1] +
          g_frac           * imp_X[, 2] +
          (-1 / sm1)       * imp_X[, 3] +
          (g_frac^2 / sm1) * imp_X[, 4]
  sum(wt_imp * (imp_Y - pred)^2)
}


estimate_feenstra_sigma_cell <- function(imp_dt, focal_importer, cfg) {
  dt <- imp_dt[importer == focal_importer]
  n_exp <- uniqueN(dt$exporter)
  max_pd <- dt[, max(period_count)]
  if (n_exp < cfg$min_exporters || max_pd < cfg$min_periods)
    return(cell_failure("insufficient_data"))

  ref_exporter <- choose_reference(dt)
  ref_vals <- dt[exporter == ref_exporter,
                 .(t, ref_ls_dif = ls_imp_dif, ref_lp_dif = lp_dif)]
  dt <- ref_vals[dt, on = "t"]
  dt <- dt[!is.na(ref_ls_dif) & !is.na(ref_lp_dif)]
  dt[, `:=`(Dk_lp = lp_dif - ref_lp_dif, Dk_ls = ls_imp_dif - ref_ls_dif)]
  dt[, `:=`(imp_y = Dk_lp^2, imp_x1 = Dk_ls^2, imp_x2 = Dk_ls * Dk_lp,
            imp_x3 = Dk_ls * lp_dif, imp_x4 = Dk_ls * ref_lp_dif,
            imp_x5 = Dk_lp * ref_lp_dif)]
  dt <- dt[!is.na(imp_y) & !is.na(imp_x1)]
  dt_nonref <- dt[exporter != ref_exporter]
  if (nrow(dt_nonref) == 0L) return(cell_failure("no_nonref_exporters"))

  setorder(dt_nonref, exporter, t)
  dt_nonref[, cusval_lag := shift(cusval, 1L), by = exporter]
  dt_nonref[, bw_w := bw_weight(cusval, cusval_lag, period_count)]
  imp_moments <- dt_nonref[,
    lapply(.SD, weighted.mean, w = bw_w, na.rm = TRUE),
    by = exporter, .SDcols = c("imp_y","imp_x1","imp_x2","imp_x3","imp_x4","imp_x5")]
  J <- nrow(imp_moments)
  if (J < 1L) return(cell_failure("no_valid_moments"))

  imp_Y_vec  <- imp_moments$imp_y
  imp_X_mat  <- as.matrix(imp_moments[, .(imp_x1,imp_x2,imp_x3,imp_x4,imp_x5)])
  wt_imp_vec <- compute_exporter_weights(dt_nonref, imp_moments$exporter, cfg)

  result <- tryCatch(
    optim(par = c(cfg$sigma_start, cfg$gamma_start),
          fn = feenstra_sigma_obj, method = "L-BFGS-B",
          lower = c(1 + 1e-6, 1e-6), upper = c(Inf, Inf),
          imp_Y = imp_Y_vec, imp_X = imp_X_mat, wt_imp = wt_imp_vec,
          control = list(maxit = 500)),
    error = function(e) NULL)

  if (is.null(result) || result$convergence != 0) {
    result <- tryCatch(
      optim(par = c(cfg$sigma_start, cfg$gamma_start),
            fn = feenstra_sigma_obj, method = "Nelder-Mead",
            imp_Y = imp_Y_vec, imp_X = imp_X_mat, wt_imp = wt_imp_vec,
            control = list(maxit = 1000)),
      error = function(e) NULL)
  }

  if (is.null(result)) return(cell_failure("optimizer_failed"))
  sig_hat <- result$par[1]; gam_hat <- result$par[2]
  if (sig_hat <= 1) sig_hat <- NA_real_
  if (gam_hat <= 0) gam_hat <- NA_real_
  list(importer = focal_importer, sigma = sig_hat, gamma = gam_hat,
       convergence = result$convergence, obj_value = result$value)
}


estimate_product_feenstra <- function(g, dt_g, cfg) {
  t0 <- proc.time()["elapsed"]
  results_g <- list(); failures_g <- list()
  n_cells <- 0L; n_ok <- 0L

  imp_stats <- dt_g[, .(n_exp = uniqueN(exporter),
                         max_pd = max(period_count)), by = importer]
  viable <- imp_stats[n_exp >= cfg$min_exporters &
                      max_pd >= cfg$min_periods, importer]

  for (imp in viable) {
    n_cells <- n_cells + 1L
    est <- tryCatch(estimate_feenstra_sigma_cell(dt_g, imp, cfg),
                    error = function(e) cell_failure(paste0("error: ", conditionMessage(e))))
    if (inherits(est, "cell_failure")) {
      failures_g[[length(failures_g) + 1L]] <- list(importer=imp, good=g, reason=est$reason)
      next
    }
    if (!is.null(est) && !is.na(est$sigma)) {
      results_g[[length(results_g) + 1L]] <- data.table(
        importer=est$importer, good=g, sigma=est$sigma, gamma=est$gamma,
        convergence=est$convergence, obj_value=est$obj_value)
      n_ok <- n_ok + 1L
    }
  }
  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  if (length(results_g) > 0L) {
    out <- rbindlist(results_g)
    attr(out, "timing") <- list(product=g, seconds=elapsed, cells=n_cells, succeeded=n_ok)
    attr(out, "failures") <- failures_g
    out
  } else NULL
}


#' Run Feenstra sigma estimation.
#' @param cfg Config list.
#' @param ncores Cores.
#' @param prepared_dt Optional pre-prepared data.table from prepare_data().
estimate_all_feenstra_sigma <- function(cfg, ncores = NULL, prepared_dt = NULL) {
  if (is.null(ncores)) ncores <- max(1L, detectCores() - 2L)
  ncores <- min(ncores, detectCores())
  cat(sprintf("FEENSTRA (1994) SIGMA ESTIMATION: %d cores\n\n", ncores))

  if (is.null(prepared_dt)) {
    prep <- prepare_data(cfg)
    dt <- prep$dt
  } else {
    dt <- prepared_dt
  }

  products <- unique(dt$good)
  n_products <- length(products)
  dt_by_product <- split(dt, by = "good", keep.by = TRUE)

  t_start <- Sys.time()
  is_windows <- .Platform$OS.type == "windows"

  if (ncores == 1L) {
    results_list <- list()
    for (idx in seq_along(products)) {
      g <- products[idx]
      if (idx %% 50 == 0 || idx == 1L) {
        el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        cat(sprintf("  [%d/%d] %.1f min elapsed\n", idx, n_products, el))
      }
      results_list[[idx]] <- estimate_product_feenstra(g, dt_by_product[[g]], cfg)
    }
  } else if (is_windows) {
    tmp_dir <- file.path(tempdir(), "het_feenstra")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    for (g in products) saveRDS(dt_by_product[[g]], file.path(tmp_dir, paste0(g, ".rds")))
    rm(dt_by_product); gc()

    cl <- makeCluster(ncores)
    on.exit({ tryCatch(stopCluster(cl), error = function(e) NULL)
              unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

    clusterExport(cl, varlist = c("estimate_product_feenstra",
      "estimate_feenstra_sigma_cell", "choose_reference", "bw_weight",
      "feenstra_sigma_obj", "compute_exporter_weights",
      "cell_failure", "cfg", "tmp_dir"),
      envir = environment())
    clusterEvalQ(cl, library(data.table))

    results_list <- parLapply(cl, products, function(g) {
      dt_g <- readRDS(file.path(tmp_dir, paste0(g, ".rds")))
      estimate_product_feenstra(g, dt_g, cfg)
    })
  } else {
    # Forked parallel path (Unix/Mac) with batch checkpointing.
    checkpoint_file <- paste0(build_output_prefix(cfg), "_feenstra_checkpoint.rds")
    batch_size <- max(ncores * 4L, 50L)
    n_batches <- ceiling(n_products / batch_size)

    results_list <- list()
    start_batch <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_batch <- ckpt$next_batch
      cat(sprintf("  Resuming Feenstra from checkpoint: %d/%d batches done\n",
                  start_batch - 1L, n_batches))
    }

    cat(sprintf("  Forked parallel: %d products in %d batches of ~%d\n\n",
                n_products, n_batches, batch_size))

    for (b in seq(start_batch, n_batches)) {
      if (b > n_batches) break
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      batch_res <- mclapply(batch_products, function(g)
        estimate_product_feenstra(g, dt_by_product[[g]], cfg),
        mc.cores = ncores)
      results_list <- c(results_list, batch_res)

      saveRDS(list(results = results_list, next_batch = b + 1L),
              checkpoint_file)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d (%.0f%%) | %.1f min\n",
                  b, n_batches, idx_s, idx_e, 100*idx_e/n_products, el))
    }
  }

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0L) { cat("\nNo sigma estimates.\n"); return(NULL) }
  output <- rbindlist(results_list)

  # Clean up Feenstra checkpoint after successful completion
  ckpt_cleanup <- paste0(build_output_prefix(cfg), "_feenstra_checkpoint.rds")
  if (file.exists(ckpt_cleanup)) {
    file.remove(ckpt_cleanup)
    cat(sprintf("  Feenstra checkpoint removed: %s\n", ckpt_cleanup))
  }

  cat(sprintf("\nFeenstra sigma: %.1f min, %s cells, sigma median=%.3f\n",
              t_elapsed, format(nrow(output), big.mark = ","),
              median(output$sigma, na.rm = TRUE)))
  output
}


# ===========================================================================
#  EXPORTER TIER CLASSIFICATION
#
#  Tier 1: Dense — full import + export side estimation.
#  Tier 2: Moderate — import-side only, with shrinkage.
#  Tier 3: Sparse — gamma assigned from regional prior.
# ===========================================================================

#' Compute destination counts per exporter for an entire product.
#'
#' Called ONCE per product in estimate_product_fixed_sigma, then
#' passed to each cell. Eliminates redundant scans of all_dt.
#'
#' @param dt_g Full product data.table (all importers).
#' @return data.table with (exporter, n_dests_total).
compute_exporter_dest_counts <- function(dt_g) {
  dt_g[, .(n_dests_total = uniqueN(importer)), by = exporter]
}


#' Classify exporters into estimation tiers.
#'
#' @param dt_nonref Non-reference exporter data for the focal importer.
#' @param exporter_dests Pre-computed data.table with (exporter, n_dests_total).
#' @param cfg Config list.
#' @return data.table with (exporter, tier).
classify_exporter_tiers <- function(dt_nonref, exporter_dests, cfg) {

  tier1_min_periods <- if (!is.null(cfg$tier1_min_periods)) cfg$tier1_min_periods else 8L
  tier1_min_dests   <- if (!is.null(cfg$tier1_min_dests))   cfg$tier1_min_dests   else 3L
  tier2_min_periods <- if (!is.null(cfg$tier2_min_periods)) cfg$tier2_min_periods else 5L

  # Periods in focal market per exporter (cheap — just this importer's data)
  imp_stats <- dt_nonref[, .(n_periods = uniqueN(t)),
                          by = exporter]

  # Merge with pre-computed destination counts
  # Subtract 1 for the focal importer (they ship there, so n_dests excludes it)
  stats <- merge(imp_stats, exporter_dests, by = "exporter", all.x = TRUE)
  stats[is.na(n_dests_total), n_dests_total := 0L]
  stats[, n_dests := pmax(n_dests_total - 1L, 0L)]

  stats[, tier := 3L]
  stats[n_periods >= tier2_min_periods, tier := 2L]
  stats[n_periods >= tier1_min_periods & n_dests >= tier1_min_dests, tier := 1L]

  stats[, .(exporter, tier)]
}


# ===========================================================================
#  FIXED-SIGMA GAMMA ESTIMATION WITH TIERED EXPORT MOMENTS
# ===========================================================================

estimate_importer_product_fixed_sigma <- function(imp_dt, focal_importer,
                                                   all_dt, cfg,
                                                   exporter_dests = NULL,
                                                   exp_lookup = NULL) {

  dt <- imp_dt[importer == focal_importer]
  n_exp <- uniqueN(dt$exporter)
  max_pd <- dt[, max(period_count)]

  if (n_exp < cfg$min_exporters || max_pd < cfg$min_periods)
    return(cell_failure("insufficient_data"))

  # --- Look up pre-estimated sigma ---
  g_code <- imp_dt$good[1]
  sigma_val <- NA_real_
  if (!is.null(cfg$sigma_lookup)) {
    match_row <- cfg$sigma_lookup[importer == focal_importer & good == g_code]
    if (nrow(match_row) > 0L) sigma_val <- match_row$sigma[1]
  }
  if (is.na(sigma_val) || sigma_val <= 1) {
    sigma_val <- if (!is.null(cfg$sigma_fallback)) cfg$sigma_fallback else
      return(cell_failure("no_sigma_estimate"))
  }

  # --- Look up shrinkage prior ---
  ln_gamma_prior <- NA_real_
  shrinkage_lambda <- if (!is.null(cfg$shrinkage_lambda)) cfg$shrinkage_lambda else 0
  if (shrinkage_lambda > 0 && !is.null(cfg$shrinkage_priors)) {
    prior_row <- cfg$shrinkage_priors[good == g_code]
    if (nrow(prior_row) > 0L) ln_gamma_prior <- prior_row$ln_gamma_prior[1]
  }

  # --- Reference-differencing and import-side moments ---
  ref_exporter <- choose_reference(dt)
  ref_vals <- dt[exporter == ref_exporter,
                 .(t, ref_ls_dif = ls_imp_dif, ref_lp_dif = lp_dif)]
  dt <- ref_vals[dt, on = "t"]
  dt <- dt[!is.na(ref_ls_dif) & !is.na(ref_lp_dif)]
  dt[, `:=`(Dk_lp = lp_dif - ref_lp_dif, Dk_ls = ls_imp_dif - ref_ls_dif)]
  dt[, `:=`(imp_y = Dk_lp^2, imp_x1 = Dk_ls^2, imp_x2 = Dk_ls * Dk_lp,
            imp_x3 = Dk_ls * lp_dif, imp_x4 = Dk_ls * ref_lp_dif,
            imp_x5 = Dk_lp * ref_lp_dif)]
  dt <- dt[!is.na(imp_y) & !is.na(imp_x1) & !is.na(imp_x2)]
  dt_nonref <- dt[exporter != ref_exporter]
  if (nrow(dt_nonref) == 0L) return(cell_failure("no_nonref_exporters"))

  # --- Tier classification (using pre-computed destination counts) ---
  if (is.null(exporter_dests)) {
    exporter_dests <- compute_exporter_dest_counts(all_dt)
  }
  tiers <- classify_exporter_tiers(dt_nonref, exporter_dests, cfg)

  tier1_exp <- tiers[tier == 1L, exporter]
  tier2_exp <- tiers[tier == 2L, exporter]
  tier3_exp <- tiers[tier == 3L, exporter]

  # Exporters to estimate: Tier 1 first, then Tier 2
  # (ordering matters for exp_jmap alignment)
  est_exporters <- c(tier1_exp, tier2_exp)
  if (length(est_exporters) == 0L) {
    # All exporters are Tier 3 — assign from prior only
    if (is.na(ln_gamma_prior)) return(cell_failure("all_tier3_no_prior"))
    gamma_prior <- exp(ln_gamma_prior)
    all_exp <- c(tiers$exporter, ref_exporter)
    return(data.table(
      importer = focal_importer, exporter = all_exp,
      sigma = sigma_val, gamma = gamma_prior,
      ref_exporter = ref_exporter, convergence = -1L,
      obj_value = NA_real_, tier = c(tiers$tier, 0L)
    ))
  }

  # --- BW weights ---
  setorder(dt_nonref, exporter, t)
  dt_nonref[, cusval_lag := shift(cusval, 1L), by = exporter]
  dt_nonref[, bw_w := bw_weight(cusval, cusval_lag, period_count)]

  # --- Import-side moments for estimated exporters only ---
  imp_moments <- dt_nonref[exporter %in% est_exporters,
    lapply(.SD, weighted.mean, w = bw_w, na.rm = TRUE),
    by = exporter,
    .SDcols = c("imp_y","imp_x1","imp_x2","imp_x3","imp_x4","imp_x5")]

  # Force ordering: Tier 1 first, then Tier 2
  imp_moments[, tier_order := match(exporter, est_exporters)]
  setorder(imp_moments, tier_order)
  imp_moments[, tier_order := NULL]

  J <- nrow(imp_moments)
  if (J < 1L) return(cell_failure("no_valid_moments"))

  exporter_order <- imp_moments$exporter
  N_tier1 <- sum(exporter_order %in% tier1_exp)

  imp_Y_vec  <- imp_moments$imp_y
  imp_X_mat  <- as.matrix(imp_moments[, .(imp_x1,imp_x2,imp_x3,imp_x4,imp_x5)])
  wt_imp_vec <- compute_exporter_weights(dt_nonref, exporter_order, cfg)

  # --- Export-side moments for Tier 1 only ---
  if (N_tier1 > 0L) {
    exp_mom <- build_export_moments(exporter_order[1:N_tier1],
                                     focal_importer, all_dt, cfg,
                                     exp_lookup = exp_lookup)
    # Align export-side weights with import-side weighting scheme.
    # exp_jmap[m] = j_idx + 2 where j_idx is the position within
    # exporter_order[1:N_tier1], so exp_jmap[m] - 2 is a valid index
    # into wt_imp_vec (Tier 1 exporters occupy positions 1:N_tier1).
    if (exp_mom$M > 0L) {
      exp_mom$wt_exp <- wt_imp_vec[exp_mom$jmap - 2L]
    }
  } else {
    exp_mom <- list(exp_Y = numeric(0), exp_X = matrix(nrow=0, ncol=9),
                    jmap = integer(0), sig_V = numeric(0),
                    gam_V = numeric(0), wt_exp = numeric(0), M = 0L)
  }

  # --- Optimization: gamma only ---
  gam_init <- cfg$gamma_start
  if (!is.null(cfg$regional_starts)) {
    imp_region <- focal_importer
    if (!is.null(cfg$regional_starts_rmap))
      imp_region <- assign_regions(focal_importer, cfg$regional_starts_rmap)
    match_row <- cfg$regional_starts[region == imp_region & good == g_code]
    if (nrow(match_row) > 0L) gam_init <- match_row$gamma[1]
  }

  d_start <- rep(gam_init, J + 1)  # gamma_k + J gamma_j
  lower_bounds <- rep(1e-6, J + 1)

  result <- tryCatch(
    optim(par = d_start, fn = het_obj_fixed_sigma, method = "L-BFGS-B",
          lower = lower_bounds, upper = rep(Inf, J + 1),
          sigma = sigma_val,
          imp_Y = imp_Y_vec, imp_X = imp_X_mat,
          exp_Y = exp_mom$exp_Y, exp_X = exp_mom$exp_X,
          exp_jmap = exp_mom$jmap,
          exp_sig_V = exp_mom$sig_V, exp_gam_V = exp_mom$gam_V,
          wt_imp = wt_imp_vec, wt_exp = exp_mom$wt_exp,
          ln_gamma_prior = ln_gamma_prior,
          shrinkage_lambda = shrinkage_lambda,
          control = list(maxit = 500)),
    error = function(e) NULL)

  if (is.null(result) || result$convergence != 0) {
    result <- tryCatch(
      optim(par = d_start, fn = het_obj_fixed_sigma, method = "Nelder-Mead",
            sigma = sigma_val,
            imp_Y = imp_Y_vec, imp_X = imp_X_mat,
            exp_Y = exp_mom$exp_Y, exp_X = exp_mom$exp_X,
            exp_jmap = exp_mom$jmap,
            exp_sig_V = exp_mom$sig_V, exp_gam_V = exp_mom$gam_V,
            wt_imp = wt_imp_vec, wt_exp = exp_mom$wt_exp,
            ln_gamma_prior = ln_gamma_prior,
            shrinkage_lambda = shrinkage_lambda,
            control = list(maxit = 1000)),
      error = function(e) NULL)
  }

  if (is.null(result)) return(cell_failure("optimizer_failed"))

  d_hat <- result$par
  gamma_k_hat <- max(d_hat[1], 0)
  gamma_j_hat <- pmax(d_hat[2:(J + 1)], 0)

  # --- Assign tiers to output ---
  est_tier <- ifelse(exporter_order %in% tier1_exp, 1L, 2L)

  # Build estimated portion
  est_dt <- data.table(
    importer     = focal_importer,
    exporter     = c(exporter_order, ref_exporter),
    sigma        = sigma_val,
    gamma        = c(gamma_j_hat, gamma_k_hat),
    ref_exporter = ref_exporter,
    convergence  = result$convergence,
    obj_value    = result$value,
    tier         = c(est_tier, 0L)  # 0 = reference exporter
  )

  # Append Tier 3 with assigned gamma
  if (length(tier3_exp) > 0L) {
    gamma_t3 <- if (!is.na(ln_gamma_prior)) exp(ln_gamma_prior) else median(gamma_j_hat)
    t3_dt <- data.table(
      importer     = focal_importer,
      exporter     = tier3_exp,
      sigma        = sigma_val,
      gamma        = gamma_t3,
      ref_exporter = ref_exporter,
      convergence  = -1L,   # flag: not estimated
      obj_value    = NA_real_,
      tier         = 3L
    )
    est_dt <- rbindlist(list(est_dt, t3_dt))
  }

  est_dt
}


#' Estimate gamma for all importers of one product with fixed sigma + tiers.
estimate_product_fixed_sigma <- function(g, dt_g, cfg) {
  t0 <- proc.time()["elapsed"]
  results_g <- list(); failures_g <- list()
  n_cells <- 0L; n_ok <- 0L; n_skipped <- 0L

  imp_stats <- dt_g[, .(n_exp = uniqueN(exporter),
                         max_pd = max(period_count)), by = importer]
  viable <- imp_stats[n_exp >= cfg$min_exporters &
                      max_pd >= cfg$min_periods, importer]
  n_skipped <- uniqueN(dt_g$importer) - length(viable)

  # Pre-compute destination counts ONCE for this product
  exporter_dests <- compute_exporter_dest_counts(dt_g)

  # Pre-split per-exporter data ONCE for this product (HS6 perf).
  # build_export_moments uses this to avoid O(N_exporters) repeated
  # filtering of dt_g on every cell.
  exp_lookup <- compute_exporter_lookup(dt_g)

  for (imp in viable) {
    n_cells <- n_cells + 1L
    est <- tryCatch(
      estimate_importer_product_fixed_sigma(dt_g, imp, dt_g, cfg,
                                             exporter_dests = exporter_dests,
                                             exp_lookup = exp_lookup),
      error = function(e) cell_failure(paste0("error: ", conditionMessage(e))))
    if (inherits(est, "cell_failure")) {
      failures_g[[length(failures_g) + 1L]] <- list(importer=imp, good=g, reason=est$reason)
      next
    }
    if (!is.null(est)) {
      est[, good := g]
      trade_wt <- dt_g[importer == imp,
                       .(avg_trade = mean(cusval, na.rm = TRUE)), by = exporter]
      est <- trade_wt[est, on = "exporter"]

      # -------------------------------------------------------------------
      #  OPTIMAL TARIFF COMPUTATION
      #
      #  Tier 3 exporters all share the same assigned gamma (the
      #  product-level regional prior). Including them in opt_tariff
      #  biases it toward the prior and masks the heterogeneous
      #  tariff that the identifying data supports.
      #
      #  Primary output (opt_tariff): computed from Tier 0/1/2 only —
      #  i.e., exporters whose gamma was directly estimated.
      #
      #  Secondary output (opt_tariff_all): computed using all exporters
      #  including Tier 3 imputations. Retained for downstream users who
      #  prefer full-coverage at the cost of bias toward the prior.
      # -------------------------------------------------------------------
      est_tiers <- if ("tier" %in% names(est)) est$tier else rep(0L, nrow(est))
      is_estimated <- !is.na(est_tiers) & est_tiers < 3L

      if (any(is_estimated)) {
        ot <- optimal_tariff(est$gamma[is_estimated],
                              est$sigma[is_estimated][1],
                              est$avg_trade[is_estimated])
      } else {
        ot <- NA_real_
      }
      ot_all <- optimal_tariff(est$gamma, est$sigma[1], est$avg_trade)

      est[, `:=`(opt_tariff = ot, opt_tariff_all = ot_all)]
      results_g[[length(results_g) + 1L]] <- est
      n_ok <- n_ok + 1L
    }
  }
  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  if (length(results_g) > 0L) {
    out <- rbindlist(results_g, fill = TRUE)
    attr(out, "timing") <- list(product=g, seconds=elapsed,
                                cells=n_cells, succeeded=n_ok, skipped=n_skipped)
    attr(out, "failures") <- failures_g
    out
  } else NULL
}


#' Run fixed-sigma gamma estimation across all products.
#' @param cfg Config list with sigma_lookup, shrinkage_lambda, shrinkage_priors.
#' @param ncores Number of cores.
#' @param prepared_dt Optional pre-prepared data from prepare_data().
estimate_all_fixed_sigma <- function(cfg, ncores = NULL, prepared_dt = NULL) {

  if (is.null(ncores)) ncores <- max(1L, detectCores() - 2L)
  ncores <- min(ncores, detectCores())
  lam <- if (!is.null(cfg$shrinkage_lambda)) cfg$shrinkage_lambda else 0
  cat(sprintf("FIXED-SIGMA GAMMA ESTIMATION: %d cores, lambda=%.3f\n\n", ncores, lam))

  if (is.null(prepared_dt)) {
    prep <- prepare_data(cfg)
    dt <- prep$dt; qlog <- prep$qlog
  } else {
    dt <- prepared_dt
    qlog <- new_quality_log()
    qlog$add("Data (pre-prepared)", n_obs = nrow(dt))
  }

  products <- unique(dt$good)
  n_products <- length(products)
  dt_by_product <- split(dt, by = "good", keep.by = TRUE)

  # Functions needed by workers
  worker_fns <- c("estimate_product_fixed_sigma",
                   "estimate_importer_product_fixed_sigma",
                   "classify_exporter_tiers",
                   "compute_exporter_dest_counts",
                   "compute_exporter_lookup",
                   "compute_exporter_weights",
                   "choose_reference", "bw_weight", "optimal_tariff",
                   "assign_regions", "build_region_map",
                   "build_export_moments", "cell_failure",
                   "het_obj_fixed_sigma")

  is_windows <- .Platform$OS.type == "windows"
  t_start <- Sys.time()

  if (ncores == 1L) {
    checkpoint_file <- paste0(build_output_prefix(cfg), "_fs_checkpoint.rds")
    results_list <- list()
    start_idx <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_idx <- ckpt$next_idx
      cat(sprintf("  Resuming from checkpoint: %d/%d products done\n",
                  start_idx - 1L, n_products))
    }
    if (start_idx <= n_products) {
      for (idx in start_idx:n_products) {
        g <- products[idx]
        if (idx %% 10 == 0 || idx == start_idx) {
          el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
          done <- idx - start_idx
          rate <- if (done > 0L) el / done else NA
          eta  <- if (done > 0L) rate * (n_products - idx) else NA
          cat(sprintf("  [%d/%d] %.1f min elapsed%s\n",
                      idx, n_products, el,
                      if (!is.na(eta)) sprintf(", ~%.1f min left", eta) else ""))
        }
        results_list[[idx]] <- estimate_product_fixed_sigma(g, dt_by_product[[g]], cfg)
        if (idx %% 50 == 0L) {
          saveRDS(list(results = results_list, next_idx = idx + 1L), checkpoint_file)
        }
      }
    }
    saveRDS(list(results = results_list, next_idx = n_products + 1L), checkpoint_file)
  } else if (is_windows) {
    tmp_dir <- file.path(tempdir(), "het_fixed_sigma")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    cat(sprintf("  Writing %d product slices to temp dir...\n", n_products))
    for (g in products) saveRDS(dt_by_product[[g]], file.path(tmp_dir, paste0(g, ".rds")))
    rm(dt_by_product); gc()

    cl <- makeCluster(ncores)
    on.exit({ tryCatch(stopCluster(cl), error = function(e) NULL)
              unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

    # Export all worker functions + config
    clusterExport(cl, varlist = c(worker_fns, "cfg", "tmp_dir"),
                  envir = environment())

    # Also export het_obj (needed by R fallback wrapper) and the Rcpp originals
    if (exists("het_obj", envir = .GlobalEnv)) {
      clusterExport(cl, varlist = "het_obj", envir = .GlobalEnv)
    }

    # Load data.table + try Rcpp on workers
    clusterEvalQ(cl, {
      library(data.table)
      # Try fixed-sigma Rcpp first
      .loaded_fs <- FALSE
      tryCatch({
        if (requireNamespace("Rcpp", quietly = TRUE)) {
          if (file.exists("het_obj_fixed_sigma_rcpp.cpp")) {
            Rcpp::sourceCpp("het_obj_fixed_sigma_rcpp.cpp")
            het_obj_fixed_sigma <- het_obj_fixed_sigma_rcpp
            .loaded_fs <- TRUE
          }
        }
      }, error = function(e) NULL)
      # If Rcpp failed, het_obj_fixed_sigma (the R wrapper) was already
      # exported via clusterExport above — workers will use it.
      # But the R wrapper needs het_obj, which also needs to be available:
      if (!.loaded_fs && !exists("het_obj")) {
        tryCatch({
          if (requireNamespace("Rcpp", quietly = TRUE) &&
              file.exists("het_obj_rcpp.cpp")) {
            Rcpp::sourceCpp("het_obj_rcpp.cpp")
            het_obj <- het_obj_rcpp
          } else {
            source("het_obj.R")
          }
        }, error = function(e) source("het_obj.R"))
      }
    })

    batch_size <- ncores * 2L
    n_batches <- ceiling(n_products / batch_size)
    checkpoint_file <- paste0(build_output_prefix(cfg), "_fs_checkpoint.rds")

    results_list <- list()
    start_batch <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_batch <- ckpt$next_batch
      cat(sprintf("  Resuming from checkpoint: %d/%d batches done\n",
                  start_batch - 1L, n_batches))
    }

    cat(sprintf("Estimating %d products in %d batches of ~%d...\n\n",
                n_products, n_batches, batch_size))

    for (b in seq(start_batch, n_batches)) {
      if (b > n_batches) break
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      batch_res <- parLapply(cl, batch_products, function(g) {
        dt_g <- readRDS(file.path(tmp_dir, paste0(g, ".rds")))
        estimate_product_fixed_sigma(g, dt_g, cfg)
      })
      results_list <- c(results_list, batch_res)

      # Checkpoint after every batch
      saveRDS(list(results = results_list, next_batch = b + 1L),
              checkpoint_file)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d (%.0f%%) | %.1f min\n",
                  b, n_batches, idx_s, idx_e, 100*idx_e/n_products, el))
    }
  } else {
    # --------------------------------------------------------------
    #  Forked parallel path (Unix/Mac) with batch checkpointing.
    #  Previously ran as a single mclapply with no resume capability —
    #  a crash at hour N of an M-hour run wasted all N hours. The
    #  batch-checkpoint scheme matches the Windows path's resilience.
    #
    #  Batch size on forked is larger than on socket because fork has
    #  effectively zero per-batch overhead (no serialization), so we
    #  want fewer, larger batches to minimize checkpoint I/O.
    # --------------------------------------------------------------
    checkpoint_file <- paste0(build_output_prefix(cfg), "_fs_checkpoint.rds")
    batch_size <- max(ncores * 4L, 50L)
    n_batches <- ceiling(n_products / batch_size)

    results_list <- list()
    start_batch <- 1L
    if (file.exists(checkpoint_file)) {
      ckpt <- readRDS(checkpoint_file)
      results_list <- ckpt$results
      start_batch <- ckpt$next_batch
      cat(sprintf("  Resuming from checkpoint: %d/%d batches done\n",
                  start_batch - 1L, n_batches))
    }

    cat(sprintf("  Forked parallel: %d products in %d batches of ~%d\n\n",
                n_products, n_batches, batch_size))

    for (b in seq(start_batch, n_batches)) {
      if (b > n_batches) break
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      batch_res <- mclapply(batch_products, function(g)
        estimate_product_fixed_sigma(g, dt_by_product[[g]], cfg),
        mc.cores = ncores)
      results_list <- c(results_list, batch_res)

      saveRDS(list(results = results_list, next_batch = b + 1L),
              checkpoint_file)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d (%.0f%%) | %.1f min\n",
                  b, n_batches, idx_s, idx_e, 100*idx_e/n_products, el))
    }
  }

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0L) { cat("\nNo estimates.\n"); return(NULL) }

  timing_info <- lapply(results_list, function(r) attr(r, "timing"))
  timing_info <- timing_info[!sapply(timing_info, is.null)]
  failure_info <- unlist(lapply(results_list, function(r) attr(r, "failures")),
                         recursive = FALSE)
  if (is.null(failure_info)) failure_info <- list()

  output <- rbindlist(results_list, fill = TRUE)
  n_raw <- nrow(output)
  n_succeeded <- length(results_list)

  # --- Post-estimation trimming ---
  # Trim bounds are computed from directly-estimated exporters only
  # (tier 0/1/2). Tier 3 rows share a clustered value at the product-level
  # prior and would bias the quantile cuts toward the prior if included.
  # The bounds are then applied to all rows.
  n_trim_total <- 0L
  trim_pct <- cfg$tail_trim_pct
  if (!is.na(trim_pct) && trim_pct > 0) {
    n_b <- nrow(output)
    output <- output[!is.na(sigma) & !is.na(gamma)]

    if ("tier" %in% names(output)) {
      trim_src <- output[is.na(tier) | tier < 3L]
    } else {
      trim_src <- output
    }
    # Safety fallback: if tier filtering leaves < 100 rows (e.g. thin
    # HS6 product), fall back to full distribution to avoid degenerate
    # quantiles. This is loud — we cat() a warning.
    if (nrow(trim_src) < 100L) {
      cat(sprintf("  [trim] warning: only %d non-Tier-3 rows, using full distribution\n",
                  nrow(trim_src)))
      trim_src <- output
    }

    sig_lo <- quantile(trim_src$sigma, trim_pct, na.rm = TRUE)
    sig_hi <- quantile(trim_src$sigma, 1 - trim_pct, na.rm = TRUE)
    gam_lo <- quantile(trim_src$gamma, trim_pct, na.rm = TRUE)
    gam_hi <- quantile(trim_src$gamma, 1 - trim_pct, na.rm = TRUE)
    output <- output[sigma >= sig_lo & sigma <= sig_hi &
                     gamma >= gam_lo & gamma <= gam_hi]
    n_trim_total <- n_b - nrow(output)
  }

  # Ensure tier column exists
  if (!"tier" %in% names(output)) output[, tier := NA_integer_]

  # Retain avg_trade so downstream recomputations (e.g. plateau fallback
  # in run_estimation.R) can re-weight optimal tariffs without re-reading
  # the full data.
  keep_cols <- intersect(names(output),
    c("importer","exporter","good","sigma","gamma","ref_exporter",
      "opt_tariff","opt_tariff_all","convergence","obj_value","tier",
      "avg_trade"))
  output <- output[, ..keep_cols]

  # --- Summary ---
  cat(sprintf("\nFixed-sigma estimation: %.1f min, %s estimates\n",
              t_elapsed, format(nrow(output), big.mark = ",")))
  cat(sprintf("  sigma median=%.3f, gamma median=%.3f, opt_tariff median=%.3f\n",
              median(output$sigma, na.rm=TRUE), median(output$gamma, na.rm=TRUE),
              median(output$opt_tariff, na.rm=TRUE)))
  if ("tier" %in% names(output)) {
    tier_tab <- output[, .N, by = tier]
    setorder(tier_tab, tier)
    for (i in seq_len(nrow(tier_tab))) {
      cat(sprintf("  Tier %s: %s estimates (%.1f%%)\n",
                  as.character(tier_tab$tier[i]),
                  format(tier_tab$N[i], big.mark = ","),
                  100 * tier_tab$N[i] / nrow(output)))
    }
  }

  attr(output, "run_meta") <- list(
    qlog = qlog, timing_info = timing_info, failure_info = failure_info,
    n_products = n_products, n_succeeded = n_succeeded,
    n_failed = n_products - n_succeeded, t_elapsed = t_elapsed,
    ncores = ncores, rcpp_loaded = .het_obj_fs_rcpp_loaded,
    trim_pct = trim_pct,
    trim_bounds = if (exists("sig_lo")) list(
      sig_lo=sig_lo, sig_hi=sig_hi, gam_lo=gam_lo, gam_hi=gam_hi) else NULL,
    n_pre_trim = n_raw, n_trimmed = n_trim_total,
    timestamp = Sys.time()
  )

  # Clean up checkpoint file after successful completion
  ckpt_cleanup <- paste0(build_output_prefix(cfg), "_fs_checkpoint.rds")
  if (file.exists(ckpt_cleanup)) {
    file.remove(ckpt_cleanup)
    cat(sprintf("  Checkpoint removed: %s\n", ckpt_cleanup))
  }

  output
}


cat("Core library loaded (v2: two-stage + tiered estimator).\n")
cat("  Use run_estimation.R for the three-stage pipeline.\n")
