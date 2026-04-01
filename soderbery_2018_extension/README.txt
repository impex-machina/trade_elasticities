================================================================================
 README: Heterogeneous Elasticity Estimator — BACI Edition (R / data.table)

 R implementation of:
   Soderbery (2018) "Trade Elasticities, Heterogeneity, and Optimal Tariffs"
   Journal of International Economics, 114, pp. 44-62.

 Adapted for CEPII BACI bilateral trade data.
================================================================================

Code last updated: 2026-03-30


OVERVIEW
========
  This code estimates heterogeneous import demand (sigma) and export
  supply (gamma) elasticities at the importer x exporter x product
  level using the structural estimator developed in Soderbery (2018).
  It is adapted to run directly on CEPII BACI reconciled trade data
  rather than raw UN Comtrade extracts.

  The objective function (het_obj.R, or het_obj_rcpp.cpp when Rcpp is
  available) jointly minimizes the distance
  between Leamer (1981) hyperbolae from import-side and export-side
  moment conditions (paper Eqs. 10 and 11), identifying heterogeneous
  export supply elasticities via cross-market variation in prices and
  quantities.


FILES
=====
  README.txt             This file
  het_obj.R              Objective function (Eqs. 10 and 11) — pure R
  het_obj_rcpp.cpp       Objective function — Rcpp/C++ (5-10x faster)
  feen94_het_baci.R      Core function library — no config, no runner
  run_estimation.R       Unified estimation script (regional or country)

  Legacy runner scripts (still functional, but run_estimation.R is preferred):
  run_regional.R         Regional estimation (Soderbery Table 1 aggregation)
  run_country.R          Country-level estimation (individual countries)

  File relationships:
    het_obj.R / het_obj_rcpp.cpp  <--  feen94_het_baci.R  <--  run_estimation.R

  feen94_het_baci.R tries to compile het_obj_rcpp.cpp via Rcpp at load
  time. If Rcpp or Rtools is unavailable, it falls back to het_obj.R
  silently. Both produce identical estimates — Rcpp is a speed
  optimization only. Each parallel worker also tries Rcpp independently.

  run_estimation.R uses a scope parameter ("regional" or "country")
  to control aggregation and parallelization:
    - scope = "regional":  ~20 entities, serial (~4 hours per pass)
    - scope = "country":   ~200 entities, parallelized (hours to days)

  Recommended workflow: run regional first, then country-level.
  Country-level automatically loads regional results (if available)
  to set per-cell starting values, which speeds optimizer convergence.


DEPENDENCIES
============
  data.table       Data manipulation and aggregation
  parallel         Parallel estimation (base R, no install needed)
  Rcpp             (Optional) Compiles C++ objective function for speed.
                   Requires Rtools on Windows. Falls back to pure R
                   if unavailable.
  haven            (Optional) Only needed if loading .dta files

  Base R optim() is used for optimization (L-BFGS-B with Nelder-Mead
  fallback). No other packages required.


INPUT DATA
==========
  CEPII BACI bilateral trade data, any HS revision.
  Download from: https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37

  BACI is distributed as a zip archive per HS revision, containing one
  CSV file per year. The expected column names are:

    t   Year
    i   Exporter (Comtrade numeric country code)
    j   Importer (Comtrade numeric country code)
    k   Product  (HS 6-digit code)
    v   Value    (thousands of USD, FOB)
    q   Quantity (metric tons)

  The code accepts either:
    - A directory path containing per-year CSV files (e.g., "BACI_HS92_V202601/")
    - A single stacked CSV, RDS, or DTA file

  When loading from a directory, the code matches only trade data files
  by their naming pattern (BACI_HS*_Y####_V*.csv). Metadata files such
  as country_codes_V202601.csv and product_codes_HS92_V202601.csv are
  automatically skipped and can safely remain in the same directory.

  IMPORTANT: BACI's product column (k) may be read as numeric by CSV
  parsers, which strips leading zeroes (e.g., 010110 becomes 10110).
  The code detects and corrects this automatically, but for best
  results ensure k is read as character. The load_baci() function
  handles this via fread's colClasses argument.


USAGE
=====
  1. Set scope and BACI data path in run_estimation.R:

       scope <- "regional"   # or "country"
       config$filepath <- "BACI_HS92_V202601/"

  2. Run:

       source("run_estimation.R")

     The script runs two steps automatically with sensitivity comparison.
     Regional runs serially (ncores = 1). Country-level auto-detects
     cores (detectCores() - 2).

  3. Step detection: If a step's RDS output already exists, it is loaded
     instead of re-estimated. Delete the file to force re-estimation.
     The serial path (ncores = 1) also supports checkpoint/resume within
     a step — if the run crashes, re-running picks up where it left off.

  4. Output files are named dynamically based on config parameters:

       {baci_source}_elast_{scope}_{agg_level}_{step}.{ext}

     The BACI source is parsed from the filepath (e.g., "BACI_HS92_V202601/"
     produces "baci_hs92_v202601"). Scope is "regional" or "country".
     Agg level is "hs4" or "hs6".

     Examples:
       baci_hs92_v202601_elast_regional_hs4_step1.rds
       baci_hs92_v202601_elast_regional_hs4_step2.rds
       baci_hs92_v202601_elast_regional_hs4_step2.csv
       baci_hs92_v202601_elast_regional_hs4_summary.rds
       baci_hs92_v202601_elast_regional_hs4_summary.txt
       baci_hs07_v202601_elast_country_hs6_step2.csv

     Country-level runs automatically search for matching regional
     results (same BACI source and agg level) to use as starting
     values.

  5. Output is a data.table with columns:

       importer      Importer code or region label
       exporter      Exporter code or region label
       good          HS4 product code (or HS6 if agg_level = "hs6")
       sigma         Estimated elasticity of substitution (demand)
       gamma         Estimated inverse export supply elasticity
       ref_exporter  Reference exporter used in estimation
       opt_tariff    Optimal non-cooperative tariff (Proposition 1)
       convergence   Optimizer exit code (0 = converged)
       obj_value     Final objective function value

  6. A summary report is generated at the end of each run in both
     machine-readable (_summary.rds) and human-readable (_summary.txt)
     formats. The report includes:
       - Config provenance and performance metrics
       - Data pipeline with trade value retention (cf. Soderbery p. 52)
       - Per-importer distributional statistics (cf. Soderbery Table 2)
       - Within importer-product heterogeneity (cf. Soderbery p. 52)
       - Variance decomposition of gamma (cf. Soderbery p. 52, 72%)
       - Cell failure diagnostics by reason
       - Two-step sensitivity comparison


INTERPRETING THE ESTIMATES
==========================
  Both sigma and gamma are reported as POSITIVE numbers. The economic
  content — downward-sloping demand, upward-sloping supply — is
  built into the structural equations, not expressed through signs
  in the output. This matters for downstream use.

  SIGMA (elasticity of substitution, sigma > 1)
  -----------------------------------------------
  Sigma measures how substitutable imported varieties of a good are
  within a CES demand system. It is NOT the price elasticity of
  demand. The relationships:

    Price elasticity of demand for variety j:  -sigma
    CES price index exponent:                  1 - sigma
    Expenditure share response to own price:   1 - sigma

  Example: sigma = 3.5 means a 1% price increase for German autos
  reduces US quantity demanded of German autos by 3.5%. In formulas
  requiring the price elasticity, use -sigma (negative).

  Common downstream uses and the correct transformation:
    Gravity models (trade flow elasticity):    1 - sigma
    Feenstra exact price index:                1 / (sigma - 1)
    Welfare gains from variety (Broda-Weinstein lambda ratio):
                                               1 / (sigma - 1)
    Markup over marginal cost (CES):           sigma / (sigma - 1)
    Optimal tariff under homogeneity:          1 / (sigma - 1)
      (but use the heterogeneous formula from Proposition 1 instead)

  GAMMA (inverse export supply elasticity, gamma > 0)
  -----------------------------------------------------
  Gamma governs how strongly an exporter's price responds to
  quantity changes. It is the INVERSE of the supply elasticity.
  The supply curve is: p = exp(g) * x^gamma. The relationships:

    Export supply elasticity:                  1 / gamma
    Tariff pass-through to delivered price:    1 / (1 + gamma)
    Terms of trade effect of a tariff:         gamma / (1 + gamma)

  Example: gamma = 1.0 means a supply elasticity of 1.0 (unit
  elastic), 50% tariff pass-through, and 50% terms-of-trade gain.
  Larger gamma = more inelastic supply = more importer market power.

  OPTIMAL TARIFF (opt_tariff, from Proposition 1)
  -------------------------------------------------
  Already computed using the correct structural formula (Eq. 15).
  It is the trade-weighted ratio of terms-of-trade gains to
  efficiency losses across heterogeneous exporters. This column
  can be used directly — no sign or transformation needed.

  IMPORTANT: Do not simply take 1/gamma as the optimal tariff.
  That formula holds only under homogeneous supply elasticities.
  With heterogeneity, the correct formula (Eq. 15) weights each
  exporter's contribution by trade values and accounts for
  interactions with the demand elasticity. The opt_tariff column
  implements this.


CONFIGURATION OPTIONS
=====================
  scope                   "regional" or "country" (in run_estimation.R)
  filepath                Path to BACI data (directory or single file)
  minyear                 First year to include (default: 1995 for HS92)
  maxyear                 Last year to include (default: NULL = all available).
                          Set to e.g. 2023 to exclude incomplete final year,
                          or use with minyear for sub-period estimation.
  agg_level               "hs4" (paper default) or "hs6"
  use_regions             TRUE: aggregate countries to Soderbery Table 1
                          regions; FALSE: keep individual countries.
                          Set automatically from scope in run_estimation.R.
  custom_region_map       Optional data.table(cty_code, region) to
                          override the built-in UN M49 mapping
  min_exporters           Min non-reference exporters per cell (default: 2)
  min_destinations        Min export destinations per exporter (default: 2)
  min_periods             Min years of differenced data (default: 3)
  uv_outlier_threshold    |delta ln(unit value)| cutoff (default: 2.0);
                          set to NA to skip
  sigma_start             Starting value for optimizer (default: 2.88)
  gamma_start             Starting value for optimizer (default: 0.69)
  sigma_V_default         Reference destination sigma (default: 2.88)
  gamma_V_default         Reference destination gamma (default: 0.69)
  tail_trim_pct           Fraction trimmed from each tail of both sigma
                          and gamma distributions (default: 0.005 = 0.5%).
                          Data-driven: uses percentiles of the current
                          estimation's distribution, not hardcoded values.


TWO-STEP STRUCTURAL DEFAULT CALIBRATION
=========================================
  The estimator has two types of parameters that are often conflated:

  STARTING VALUES (sigma_start, gamma_start):
    Where the optimizer begins searching. Affect convergence speed
    and which local minimum is found. Do NOT enter the moment
    equations. Can be updated freely without methodological concern.

  STRUCTURAL DEFAULTS (sigma_V_default, gamma_V_default):
    The reference destination's elasticities, which enter as FIXED
    COEFFICIENTS in the export-side moment equations (Eq. 11). These
    are not estimated — they are treated as known. Incorrect values
    introduce bias in gamma estimates.

  Soderbery sets structural defaults at his global medians from
  Comtrade (sigma=2.88, gamma=0.69). These are endogenous to his
  estimation — he just doesn't iterate, so the circularity is
  hidden. Using his Comtrade values for BACI data is inappropriate
  because the two datasets have different reconciliation procedures,
  sample periods, and unit value construction.

  Both runner scripts use a two-step procedure to set data-appropriate
  structural defaults:

    Step 1 (Discovery): Estimate with paper defaults (2.88/0.69)
      as structural parameters. The purpose is to discover what
      sigma and gamma look like under BACI data.

    Step 2 (Final): Fix structural defaults at Step 1 medians.
      Update starting values for faster convergence. Re-estimate.
      These are the final results. No further updates.

  This is methodologically identical to Soderbery — a single pass
  with fixed structural defaults. The only difference is that the
  defaults are derived from the same dataset rather than borrowed
  from Comtrade. Structural defaults are updated exactly ONCE (from
  Step 1 to Step 2), not iterated. Iteration of structural defaults
  is unstable: each pass produces higher estimates which feed back
  as higher coefficients, creating a positive feedback loop with no
  fixed point.

  The sensitivity comparison at the end of each script reports the
  shift between Step 1 and Step 2 medians and the cross-step
  correlation, documenting how much the structural defaults matter.

  For country-level estimation, run_country.R also loads regional
  results (if available) to set per-cell optimizer starting values.
  This is a speed optimization only — it affects where the optimizer
  begins searching, not the structural parameters.


ADAPTATIONS FROM ORIGINAL (COMTRADE) VERSION
=============================================
  The core estimation logic — objective function, moment construction,
  reference-country selection, BW weighting, and optimization — is
  unchanged. The following adaptations handle differences between
  raw Comtrade and BACI:

  1. COLUMN MAPPING
     BACI uses i=exporter, j=importer (opposite of the paper's
     convention where i=importer, j=exporter). The config block
     maps these correctly.

  2. HS6 TO HS4 AGGREGATION
     BACI provides HS6-level data. The paper estimates at HS4.
     The code aggregates by summing values and quantities within
     each (importer, exporter, HS4, year) group before computing
     unit values. The pad_hs6() function ensures leading zeroes
     are preserved: if k was read as numeric (e.g., 10110 instead
     of 010110), it zero-pads to 6 digits before extracting the
     HS4 prefix.

  3. SINGLE RECONCILED VALUE
     BACI provides one FOB-reconciled value per bilateral flow,
     whereas Soderbery's original uses both importer-reported (CIF)
     and exporter-reported (FOB) values from Comtrade. The code
     computes import shares and export shares from the same
     reconciled value using different denominators, which is the
     source of cross-market variation for identification. The paper
     (p. 49, fn 11) explicitly notes this approach is valid.

  4. REGIONAL AGGREGATION
     The build_region_map() function maps Comtrade numeric country
     codes to the 20 entities in Soderbery Table 1 (13 individual
     countries + 7 UN M49-based regions). Users should verify the
     mapping against their BACI country_codes metadata file.

  5. UNIT VALUE OUTLIER FILTER
     Observations with |delta ln(unit value)| >= 2.0 are dropped
     after first-differencing. This is standard in the Feenstra/
     Broda-Weinstein/Soderbery tradition and removes measurement
     error from quantity misreporting.

  6. POST-ESTIMATION TRIMMING
     A symmetric percentile trim (default: 0.5% each tail) is applied
     to both sigma and gamma distributions. This is data-driven —
     percentiles are computed from the current estimation's results
     rather than hardcoded to values from the paper's Comtrade sample
     (e.g., the paper's sigma cap of 131.05 was its 99.5th percentile).
     Applying the same logic to both parameters ensures consistency.

  7. TWO-STEP STRUCTURAL DEFAULT CALIBRATION
     The sigma_V_default and gamma_V_default parameters enter as fixed
     coefficients in the export-side moment equations (Eq. 11), not just
     as optimizer starting points. Using the paper's Comtrade values
     (2.88/0.69) for BACI data is inappropriate because the datasets
     differ in reconciliation, sample period, and unit values. Both
     runner scripts use a two-step procedure: Step 1 estimates with
     paper defaults to discover BACI-specific medians; Step 2 fixes
     structural defaults at those medians and re-estimates. Structural
     defaults are updated exactly once, not iterated, because iteration
     produces a positive feedback loop with no fixed point.

  8. DATA QUALITY DIAGNOSTICS
     A quality log tracks observations at each stage: raw load,
     positive-value filter, year filter, HS4 aggregation, regional
     aggregation, first-differencing, outlier filter, estimation
     outcomes, and post-estimation trimming. A formatted report
     prints at completion.


REGIONAL MAPPING NOTES
======================
  The built-in region map has been validated against BACI's
  country_codes_V202601.csv. It covers 233 of 238 BACI country
  codes, including historical entities (e.g., Czechoslovakia,
  East/West Germany, Serbia and Montenegro) and overseas
  territories (e.g., Macao, Mayotte, Curacao, New Caledonia).

  The 5 unmapped codes are assigned to "OTHER":
    86   Br. Indian Ocean Terr.
    260  Fr. South Antarctic Terr.
    316  Guam
    612  Pitcairn
    810  USSR (...1990) — excluded by minyear = 1995

  EDGE CASES AND DEPARTURES FROM STRICT M49:
    Soderbery's regions are composites (e.g., "Northern/Western
    Europe" bundles M49 Northern + Western + parts of Southern
    Europe). The following assignments depart from strict M49
    sub-region classification but are consistent with the paper's
    Table 6 statistics and the economic logic of trade patterns:

    - Spain (724), Portugal (620), Greece (300), Malta (470):
      Assigned to NWU (Northern/Western Europe), not M49 Southern
      Europe. These are established EU members with low tariffs
      aligned with the EU common external tariff.

    - Turkey (792): Assigned to SEU (Southern/Eastern Europe),
      not M49 Western Asia. Turkey's EU customs union membership
      and European-oriented trade patterns support this placement.

    - Tajikistan (762): Assigned to ASA (Asian). M49 classifies
      Tajikistan as Central Asia (code 143), under Asia.

  To customize:

    1. Load BACI's country_codes_V202601.csv (or equivalent)
    2. Cross-reference with build_region_map() output
    3. Supply corrections via config$custom_region_map

  Example:
    my_map <- fread("my_region_mapping.csv")  # cols: cty_code, region
    config$custom_region_map <- my_map


COMPUTATIONAL NOTES
===================
  - Regional estimation runs serially (ncores = 1). With only ~20
    importers per product, each product estimates in seconds —
    parallelization overhead exceeds the computation time. Expect
    ~2 hours per step with Rcpp (~4 hours total for both steps).

  - Country-level estimation runs in parallel across products.
    Cores are set automatically to detectCores() - 2. With 10
    cores and BACI's 30-year panel, expect 1-3 days per step.
    On Windows with 32GB RAM, 4-6 cores is typically the maximum
    before memory pressure causes failures.

  - Both steps detect previously completed results via the output
    RDS file. If the file exists, the step is skipped and results
    are loaded from disk. Delete the file to force re-estimation.

  SPEED OPTIMIZATIONS:

  - Rcpp objective function (het_obj_rcpp.cpp): The objective
    function is called hundreds of times per cell by optim().
    The C++ implementation via Rcpp is 5-10x faster per evaluation
    than pure R. feen94_het_baci.R detects and compiles it
    automatically at load time; each parallel worker does the
    same. Requires Rcpp package and Rtools (Windows) or Xcode
    command line tools (Mac). Falls back to pure R silently.

  - Regional-informed starting values: Country-level runs
    automatically load regional results (matching the same BACI
    source and agg level) if present and use them to set per-cell
    optimizer starting values. Instead of starting every cell at
    global defaults, each cell starts at its region-product median.
    This reduces optimizer iterations by 30-50%.

  - Pre-filtering: Before calling the estimation function for each
    importer, the code checks whether the cell meets minimum data
    requirements (min_exporters, min_periods). Cells that cannot
    meet the requirements are skipped without entering the function,
    avoiding unnecessary data subsetting and moment construction.

  - Recommended workflow: run scope="regional" first (~4 hours),
    then scope="country". The country run benefits from Rcpp and
    regional-informed starting values automatically.

  CHECKPOINT AND RESUME:

  - The serial path (ncores = 1) saves incremental checkpoints
    every 50 products to {out_prefix}_checkpoint.rds. If the run
    crashes, re-running the same script resumes from the last
    checkpoint. Only one checkpoint file exists at a time (each
    save overwrites the previous). The checkpoint is deleted after
    successful completion.

  - The checkpoint system combines with step detection: if Step 1
    completes and saves _step1.rds, subsequent runs skip Step 1
    entirely. If Step 2 crashes mid-run, it resumes from its own
    checkpoint.

  WINDOWS PARALLEL NOTES:

  - Windows uses socket clusters (makeCluster/parLapply), which
    require data serialization to worker processes. To avoid
    memory failures with large datasets, the code writes product
    slices to temporary RDS files and has workers read from disk
    rather than receiving data through the socket.

  - On Linux/Mac, forked parallelism (mclapply) shares memory
    and avoids this overhead entirely. If running on Linux with
    64+ cores and 128GB+ RAM, country-level estimation completes
    in 2-4 hours.

  OTHER NOTES:

  - Progress reporting: Serial path reports every 10 products
    with elapsed time and ETA. Parallel path reports per batch.

  - The optimizer uses L-BFGS-B (bounded, gradient-based) with a
    Nelder-Mead fallback (maxit=1000) for cells where L-BFGS-B
    fails to converge.

  - Cell failure diagnostics: When cells fail to produce estimates,
    the reason is logged (insufficient_data, no_nonref_exporters,
    no_valid_moments, optimizer_failed, or error with message).
    Failure counts by reason are printed to the console and
    included in the summary report.

  - The convergence column in the output flags non-converged cells
    (convergence != 0). Consider excluding these from analysis or
    investigating individually.

  - Config validation runs before data loading to catch obvious
    misconfigurations (missing fields, invalid paths, inconsistent
    year ranges) before expensive operations.


CITATION
========
  Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
  Tariffs," Journal of International Economics, 114, 2018, pp. 44-62.
