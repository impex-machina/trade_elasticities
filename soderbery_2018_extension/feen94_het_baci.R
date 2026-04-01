#' ===========================================================================
#' feen94_het_baci.R
#'
#' Core function library for Soderbery (2018) heterogeneous elasticity
#' estimation, adapted for CEPII BACI data.
#'
#' This file contains ALL shared functions. It does NOT contain config
#' or estimation calls. Source this, then use one of:
#'   - run_regional.R   (regional aggregation, Soderbery Table 1)
#'   - run_country.R    (individual country pairs)
#'
#' CITATION:
#'   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
#'   Tariffs," JIE, 114, 2018, pp. 44-62.
#'
#' Last updated: 2026-03-30
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
#' @param cfg Config list.
#' @return Named list with exp_Y, exp_X, exp_jmap, sig_V, gam_V, wt_exp.
build_export_moments <- function(exporter_order, focal_importer, all_dt, cfg) {

  # ---------------------------------------------------------------
  # PURPOSE: Construct the export-side moment conditions (Eq. 11)
  # for each non-reference exporter j that sells to the focal
  # importer I. Identification comes from comparing exporter j's
  # price/share movements in market I against its movements in a
  # *reference destination* V (the largest alternative market for
  # exporter j, excluding I). This cross-market variation pins
  # down the export supply elasticity gamma_j separately from
  # the demand elasticity sigma.
  #
  # RETURNS: A list with export-side Y vector, X matrix, and
  # metadata (jmap, sig_V, gam_V, weights). If no exporters have
  # enough cross-market data, returns zero-length structures so
  # the optimizer falls back to import-side moments only.
  # ---------------------------------------------------------------

  exp_Y_list   <- list()
  exp_X_list   <- list()
  exp_jmap_vec <- integer(0)   # Maps each export obs to its position in d[]
  sig_V_vec    <- numeric(0)   # sigma at reference destination V (fixed)
  gam_V_vec    <- numeric(0)   # gamma at reference destination V (fixed)

  for (j_idx in seq_along(exporter_order)) {
    exp_j <- exporter_order[j_idx]

    # --- Check that exporter j sells to enough destinations ---
    # Eq. 11 requires variation across destinations. If j only
    # exports to 1 market, there's no reference destination to
    # difference against.
    exp_flows <- all_dt[exporter == exp_j]
    n_dest <- uniqueN(exp_flows$importer)
    if (n_dest < cfg$min_destinations) next

    # --- Select reference destination V ---
    # V is the largest market (by total trade value) for exporter j,
    # EXCLUDING the focal importer I. This ensures the reference
    # market is economically significant and avoids using I as its
    # own reference.
    dest_stats <- exp_flows[importer != focal_importer,
                            .(dest_val = sum(cusval, na.rm = TRUE)),
                            by = importer]
    if (nrow(dest_stats) == 0L) next

    ref_dest <- dest_stats$importer[which.max(dest_stats$dest_val)]

    # --- Extract reference destination time series ---
    # These are the first-differenced log prices and log export
    # shares of exporter j in market V. They enter Eq. 11 as the
    # "V" variables (D lp_V, D ls_V).
    ref_dest_vals <- exp_flows[importer == ref_dest,
                               .(t, ref_lp_exp = lp_dif,
                                 ref_ls_exp = ls_exp_dif)]

    # --- Merge with focal importer time series ---
    # Join on time period so each obs has BOTH j's movements in
    # market I (lp_dif, ls_exp_dif) and j's movements in market V
    # (ref_lp_exp, ref_ls_exp) for the same year. The LHS of
    # Eq. 11 is the squared price difference (lp_I - lp_V)^2.
    focal_vals <- exp_flows[importer == focal_importer]
    focal_vals <- ref_dest_vals[focal_vals, on = "t"]
    focal_vals <- focal_vals[!is.na(ref_lp_exp) & !is.na(ref_ls_exp) &
                             !is.na(lp_dif) & !is.na(ls_exp_dif)]

    if (nrow(focal_vals) == 0L) next

    # --- Build moment variables (9 RHS terms for Eq. 11) ---
    # These are the cross-products of differenced log shares and
    # prices across the two markets (I and V). Each term corresponds
    # to a coefficient in the Leamer hyperbola that is a nonlinear
    # function of sigma, gamma_I (the exporter's supply elasticity
    # in market I), and gamma_V / sigma_V (the structural defaults
    # at the reference destination, treated as known).
    focal_vals[, `:=`(
      exp_y  = (lp_dif - ref_lp_exp)^2,        # LHS: squared price gap I vs V
      exp_x1 = ls_exp_dif^2,                    # (D ls_I)^2
      exp_x2 = ls_exp_dif * lp_dif,             # (D ls_I)(D lp_I)
      exp_x3 = ref_ls_exp^2,                    # (D ls_V)^2
      exp_x4 = ref_ls_exp * lp_dif,             # (D ls_V)(D lp_I)
      exp_x5 = ref_ls_exp * ref_lp_exp,         # (D ls_V)(D lp_V)
      exp_x6 = ls_exp_dif * ref_lp_exp,         # (D ls_I)(D lp_V)
      exp_x7 = ls_exp_dif * ref_ls_exp,         # (D ls_I)(D ls_V)
      exp_x8 = ref_lp_exp^2,                    # (D lp_V)^2
      exp_x9 = ref_lp_exp * lp_dif              # (D lp_V)(D lp_I)
    )]
    focal_vals <- focal_vals[!is.na(exp_y)]
    if (nrow(focal_vals) == 0L) next

    # --- BW-weighted time average ---
    # Same Broda-Weinstein weighting as import side: down-weights
    # obs with large period-to-period value swings.
    setorder(focal_vals, t)
    focal_vals[, `:=`(cusval_lag = shift(cusval, 1L), pd_e = .N)]
    focal_vals[, bw_w_e := bw_weight(cusval, cusval_lag, pd_e)]

    exp_cols <- c("exp_y","exp_x1","exp_x2","exp_x3","exp_x4",
                  "exp_x5","exp_x6","exp_x7","exp_x8","exp_x9")
    exp_mom <- focal_vals[,
      lapply(.SD, weighted.mean, w = bw_w_e, na.rm = TRUE),
      .SDcols = exp_cols
    ]

    # --- Collect into output structures ---
    exp_Y_list[[length(exp_Y_list) + 1L]] <- exp_mom$exp_y
    exp_X_list[[length(exp_X_list) + 1L]] <- as.numeric(
      exp_mom[, .(exp_x1,exp_x2,exp_x3,exp_x4,exp_x5,
                  exp_x6,exp_x7,exp_x8,exp_x9)]
    )
    # j_idx + 2 maps this exporter to its position in the parameter
    # vector d[]: d[1]=sigma, d[2]=gamma_k, d[3:]=gamma_j's.
    # So exporter j_idx in exporter_order corresponds to d[j_idx+2].
    exp_jmap_vec <- c(exp_jmap_vec, j_idx + 2L)
    # sigma_V and gamma_V at the reference destination are treated
    # as FIXED KNOWN VALUES in Eq. 11 — they are structural defaults
    # from config, not estimated. This is where the two-step
    # calibration matters: Step 1 uses paper defaults, Step 2
    # replaces them with BACI-specific medians from Step 1.
    sig_V_vec    <- c(sig_V_vec, cfg$sigma_V_default)
    gam_V_vec    <- c(gam_V_vec, cfg$gamma_V_default)
  }

  M <- length(exp_Y_list)

  # Return M export-side observations stacked into a single Y
  # vector and X matrix. If M=0, the objective function receives
  # zero-length inputs and uses import-side moments only (still
  # identified for sigma, but gamma_j estimates rely on import-
  # side variation alone).
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


#' Estimate one (importer, product) cell.
#'
#' Returns a data.table on success, or a cell_failure object with
#' a diagnostic reason on failure. The cell_failure is collected
#' by estimate_product for the failure log.
estimate_importer_product <- function(imp_dt, focal_importer, all_dt, cfg) {

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

  exp_mom <- build_export_moments(exporter_order, focal_importer, all_dt, cfg)
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

  for (imp in viable_importers) {
    n_cells <- n_cells + 1L
    est <- tryCatch(
      estimate_importer_product(dt_g, imp, dt_g, cfg),
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

prepare_data <- function(cfg) {

  validate_config(cfg)

  qlog <- new_quality_log()
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
                   "build_export_moments", "cell_failure")

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
    cat(sprintf("Estimating %d products in %d batches of ~%d...\n\n",
                n_products, n_batches, batch_size))
    results_list <- list()

    for (b in seq_len(n_batches)) {
      idx_s <- (b - 1L) * batch_size + 1L
      idx_e <- min(b * batch_size, n_products)
      batch_products <- products[idx_s:idx_e]

      # Send only product names — workers read their own data from disk
      batch_res <- parLapply(cl, batch_products, function(g) {
        dt_g <- readRDS(file.path(tmp_dir, paste0(g, ".rds")))
        estimate_product(g, dt_g, cfg)
      })
      results_list <- c(results_list, batch_res)

      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("  Batch %d/%d: products %d-%d of %d (%.0f%%) | %.1f min, ~%.1f min left\n",
                  b, n_batches, idx_s, idx_e, n_products,
                  100 * idx_e / n_products, el, el / idx_e * (n_products - idx_e)))
    }

  } else {
    cat("Using forked parallelism (Unix/Mac)...\n")
    results_list <- mclapply(products, function(g)
      estimate_product(g, dt_by_product[[g]], cfg), mc.cores = ncores)
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
    return(tolower(gsub("_", "_", m)))  # e.g., "baci_hs92_v202601"
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


cat("Core library loaded. Use run_regional.R or run_country.R to estimate.\n")
