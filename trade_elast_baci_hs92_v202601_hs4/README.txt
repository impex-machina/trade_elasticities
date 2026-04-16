================================================================================
 README: Heterogeneous Elasticity Estimator — BACI Edition (R / data.table)
         Three-stage pipeline: Feenstra sigma -> regional gamma -> country gamma

 R implementation adapted from:
   Soderbery (2018) "Trade Elasticities, Heterogeneity, and Optimal Tariffs"
   Journal of International Economics, 114, pp. 44-62.

 Applied to CEPII BACI bilateral trade data.
================================================================================

Code last updated: 2026-04-16


OVERVIEW
========
  This code estimates heterogeneous import demand (sigma) and export
  supply (gamma) elasticities at the importer x exporter x product
  level using a three-stage extension of Soderbery's (2018) structural
  estimator, adapted to run directly on CEPII BACI reconciled trade
  data at the country level rather than only Soderbery's regional
  aggregation.

  METHODOLOGICAL DEPARTURE FROM SODERBERY (2018)

  Soderbery's joint nonlinear SUR estimator works at his 20 regional
  entities but fails at the country level (~233 entities) due to two
  distinct identification pathologies:

    1. The sigma-gamma ridge. The joint objective function has a
       degenerate manifold where sigma and gamma inflate together
       along directions that leave moment conditions approximately
       satisfied. At the country level on BACI 1995-2024 HS4, the
       joint estimator produces pathological values (sigma median
       27.3, gamma median 877, with 63% of cells on the plateau).

    2. The gamma plateau. Even with sigma fixed, gamma enters the
       moment conditions through gamma/(1+gamma), which saturates
       above gamma ~= 50. The objective surface becomes flat in
       gamma, so the optimizer drifts into the tail without
       corresponding gains in fit. This pathology is present even
       at Soderbery's 20-entity regional aggregation on BACI data
       (41% of regional cells land above gamma = 20 without
       shrinkage; see empirical validation below).

  This code implements a three-stage pipeline that resolves both
  pathologies:

    STAGE 1: Feenstra (1994) sigma + gamma_common per (importer,
      product). Two-parameter estimation from import-side moments
      under the homogeneity assumption. Well-identified. Breaks
      the ridge by fixing sigma for subsequent stages.

    STAGE 2a: Regional gamma with fixed sigma and moderate
      hierarchical shrinkage (lambda = 0.05 at HS4) anchored to
      Stage 1's gamma_common. The shrinkage prevents the plateau
      at regional level — empirically required on BACI data even
      at 20-entity aggregation.

    STAGE 2b: Country gamma with fixed sigma and shrinkage
      (lambda = 0.1 at HS4) anchored to Stage 2a regional estimates.
      Dense bilateral cells speak through the moments; thin cells
      are regularized toward the regional prior. A Tier 3 safety
      net assigns the regional prior directly to exporters with
      insufficient identification; these rows are separated from
      primary optimal-tariff calculations.

  LAMBDA CALIBRATION AT HS4. Stage 2a lambda was chosen empirically
  via a four-point sweep {0.01, 0.02, 0.05, 0.1}. The R-squared and
  within-pair MAD targets from Soderbery (0.72 and 0.125) trace out
  a trade-off curve along lambda rather than a single point on this
  sample. Selection at lambda = 0.05 minimizes plateau share while
  preserving meaningful within-cell heterogeneity (see DIAGNOSTICS
  below). Stage 2b lambda = 0.1 produces country-level MAD = 0.104,
  close to Soderbery's 0.125 benchmark.

  A fourth design choice — trade-value weighting of exporter moment
  rows — applies across all three stages. Weighting exporters by
  cell trade value rather than uniformly brings the variance
  decomposition R-squared into agreement with Soderbery's benchmark.
  This matches the "large and persistent" intuition in Soderbery
  (p. 50, fn 14), which the paper applies within exporter over time
  via BW weights but leaves implicit across exporters. A period-count
  floor is applied to trade-value weights to protect against lumpy
  single-year shipments (aircraft, LNG cargoes, rough diamonds)
  dominating the objective via raw trade value.

  REFERENCE-DESTINATION ELASTICITIES. Soderbery's Eq. (11) export-
  side moments depend on (sigma_V, gamma_V) at each Tier 1 exporter's
  reference destination V. This code uses per-(destination, product)
  lookups from Stage 1 (sigma_V) and Stage 2a (gamma_V), with
  fallback to global defaults when no match exists. The baseline
  Soderbery treatment uses a single global default for all exporters;
  the bilateral lookup is a structurally correct refinement
  exploiting BACI's granularity.


FILES
=====
  README.txt                          This file
  het_obj.R                           Joint objective function (Eqs. 10+11) — pure R
  het_obj_rcpp.cpp                    Joint objective function — Rcpp/C++ (5-10x faster)
  het_obj_fixed_sigma_rcpp.cpp        Fixed-sigma objective with shrinkage penalty (Rcpp)
  feen94_het_baci.R                   Core function library
  run_est_baci_hs92_v202601_hs4.R                    Three-stage runner

  File relationships:

    het_obj.R / het_obj_rcpp.cpp                 (original joint estimator)
    het_obj_fixed_sigma_rcpp.cpp                 (fixed-sigma gamma estimator)
         |                 |
         +--- feen94_het_baci.R ---+
                    |
              run_est_baci_hs92_v202601_hs4.R

  feen94_het_baci.R compiles both Rcpp files at load time. If Rcpp
  or Rtools is unavailable, it falls back silently: het_obj.R for
  the joint objective, and an R-level wrapper that reuses het_obj
  for the fixed-sigma case. Both paths produce identical estimates.

  run_est_baci_hs92_v202601_hs4.R is the single entry point. It detects previously
  completed stages (raw cache, Stage 1 sigma, Stage 2a regional,
  Stage 2b country) via their output files and skips them when
  present. Delete an output file to force re-estimation of that
  stage. All downstream stages re-run if their input stage changes.


DEPENDENCIES
============
  data.table       Data manipulation and aggregation
  parallel         Parallel estimation (base R, no install needed)
  Rcpp             (Recommended) Compiles C++ objective functions.
                   Requires Rtools on Windows or Xcode CLT on Mac.
                   Falls back to pure R if unavailable.
  haven            (Optional) Only needed if loading .dta files

  Base R optim() is used for optimization (L-BFGS-B with Nelder-Mead
  fallback). No other packages required.


WORKING DIRECTORY
=================
  run_est_baci_hs92_v202601_hs4.R locates feen94_het_baci.R and the Rcpp/het_obj
  files via a three-step precedence:

    1. Environment variable SODERBERY_DIR if set and exists.
       Recommended for production:
         export SODERBERY_DIR=/home/ubuntu/estimation

    2. Directory of the script itself, when invoked via Rscript.

    3. Current working directory (unchanged).

  This replaces the previous hardcoded setwd() call.


INPUT DATA
==========
  CEPII BACI bilateral trade data, any HS revision.
  Download from: https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37

  BACI is distributed as a zip archive per HS revision, containing one
  CSV file per year. Expected column names:

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
  The code detects and corrects this automatically.


USAGE
=====
  1. Set the BACI data path in run_est_baci_hs92_v202601_hs4.R:

       config$filepath <- "BACI_HS92_V202601/"

  2. Set the working directory (once per environment):

       export SODERBERY_DIR=/home/ubuntu/estimation

  3. Run:

       Rscript run_est_baci_hs92_v202601_hs4.R

     The script runs all three stages automatically. Each stage
     detects previously completed output and skips if present.

  4. Stage detection works via output files:

       {prefix}_raw_cache.rds                      Raw cleaned data
       {prefix}_feenstra_sigma.rds                 Stage 1 sigma/gamma_common
       {regional_prefix}_fixed_sigma.rds           Stage 2a regional gamma
       {prefix}_fixed_sigma.rds                    Stage 2b country gamma

     Delete any file to force re-estimation of that stage AND
     all downstream stages. The raw cache should be deleted if
     the BACI data changes (e.g., a new release). Stage 1 should
     be deleted if sigma_start/gamma_start change. Stage 2a
     should be deleted if shrinkage_lambda changes at either stage.

     All stages support per-batch checkpointing on both the forked
     (Linux/Mac) and socket (Windows) parallel paths. Checkpoint
     files are removed on successful completion. A crashed run
     resumes from the last completed batch on re-run.

  5. Output files are named dynamically based on config parameters:

       {baci_source}_elast_{scope}_{agg_level}_{artifact}.{ext}

     The BACI source is parsed from the filepath (e.g., "BACI_HS92_V202601/"
     produces "baci_hs92_v202601"). Scope is "regional" or "country".
     Agg level is "hs4" or "hs6".

     Examples of complete outputs after a successful run:

       baci_hs92_v202601_elast_country_hs4_raw_cache.rds
       baci_hs92_v202601_elast_country_hs4_feenstra_sigma.rds
       baci_hs92_v202601_elast_regional_hs4_fixed_sigma.rds
       baci_hs92_v202601_elast_regional_hs4_fixed_sigma.csv
       baci_hs92_v202601_elast_regional_hs4_summary.rds
       baci_hs92_v202601_elast_regional_hs4_summary.txt
       baci_hs92_v202601_elast_country_hs4_fixed_sigma.rds
       baci_hs92_v202601_elast_country_hs4_fixed_sigma.csv
       baci_hs92_v202601_elast_country_hs4_summary.rds
       baci_hs92_v202601_elast_country_hs4_summary.txt

  6. Output is a data.table with columns:

       importer        Importer code or region label
       exporter        Exporter code or region label
       good            HS4 product code (or HS6 if agg_level = "hs6")
       sigma           Import demand elasticity of substitution (from Stage 1)
       gamma           Inverse export supply elasticity (from Stage 2a/2b)
       ref_exporter    Reference exporter used in import-side estimation
       opt_tariff      Optimal non-cooperative tariff (Proposition 1),
                       computed from directly-estimated exporters only
                       (Tier 0/1/2)
       opt_tariff_all  Optimal tariff including Tier 3 imputed exporters.
                       Retained for full-coverage analyses at the cost
                       of bias toward the regional prior.
       avg_trade       Average trade value within the cell, retained for
                       downstream re-weighting.
       convergence     Optimizer exit code (see below).
       obj_value       Final objective function value (NA for Tier 3
                       and plateau-fallback rows).
       tier            Estimation tier (see below).

     CONVERGENCE CODES:
        0 = Optimizer converged normally
        1 = Max iterations reached (estimates still usable if shrinkage
            keeps gamma near the prior; check distribution before use)
       -1 = Tier 3 assignment: gamma set from regional prior, not estimated
       -2 = Plateau fallback at Stage 2a: gamma estimate exceeded
            plateau threshold (>20) and was replaced with the Feenstra
            product-level anchor (Stage 2a only)

     TIER COLUMN:
       0 = Reference exporter (import-side reference within the cell)
       1 = Full import + export side estimation
       2 = Import-side only estimation (moderate density)
       3 = Assigned from product-level regional prior (thin cells)

     Tier thresholds are configurable (see CONFIGURATION OPTIONS below).
     With default thresholds set to match Soderbery's minimum data
     requirements, most estimated cells are Tier 1.

  7. Summary reports are generated for BOTH regional (Stage 2a)
     and country (Stage 2b) outputs, in both machine-readable
     (_summary.rds) and human-readable (_summary.txt) formats.
     Each report includes:
       - Config provenance and performance metrics
       - Data pipeline with trade value retention
       - Per-importer distributional statistics (cf. Soderbery Table 2)
       - Within importer-product heterogeneity (cf. Soderbery p. 52)
       - Variance decomposition of gamma (cf. Soderbery p. 52)
       - Cell failure diagnostics by reason


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
    Welfare gains from variety:                1 / (sigma - 1)
    Markup over marginal cost (CES):           sigma / (sigma - 1)

  NOTE ON SIGMA MAGNITUDE. Because sigma comes from Feenstra's
  import-side-only estimator under the homogeneity assumption, it
  is structurally higher than Soderbery's joint estimates (HS4 BACI
  1995-2024 median 3.88 vs Soderbery's 2.88). This is the known
  upward bias from homogeneity when gamma is in fact heterogeneous.
  The bias is bounded (gamma/(1+gamma) is bounded above by 1, which
  caps the bias), predictable in direction (always upward), and
  empirically small relative to the catastrophic failure of joint
  estimation at disaggregate resolution.

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
  efficiency losses across heterogeneous exporters.

  PRIMARY opt_tariff is computed from directly-estimated exporters
  only (Tier 0/1/2). A secondary opt_tariff_all includes Tier 3
  imputations for downstream users who want full coverage at the
  cost of bias toward the regional prior. On HS4 these medians
  differ by less than 1%; the divergence is larger at HS6 where
  Tier 3 shares are expected to be higher.

  IMPORTANT: Do not simply take 1/gamma as the optimal tariff.
  That formula holds only under homogeneous supply elasticities.
  With heterogeneity, the correct formula (Eq. 15) weights each
  exporter's contribution by trade values and accounts for
  interactions with the demand elasticity. The opt_tariff column
  implements this.


DIAGNOSTICS FROM HS4 BACI 1995-2024
====================================
  HEADLINE STATISTICS:

    Stage 1 (Feenstra sigma):
      Cells:                          243,035 (of 243,036 attempts)
      Sigma median:                   3.884
      Sigma IQR:                      [2.769, 6.084]
      Gamma_common median:            1.000

    Stage 2a (regional, lambda=0.05):
      Estimates:                      439,913
      Gamma median:                   1.075
      Plateau fallback invoked:       2,081 cells (0.47%)
      Tier 1 share:                   88.5%

    Stage 2b (country, lambda=0.10):
      Estimates:                      8,082,613
      Gamma median:                   1.223
      Gamma IQR:                      [0.710, 1.986]
      Gamma > 20 share:               0.9% (Tier 3 values at high-gamma
                                      products; not plateau escape)
      Opt_tariff median:              1.045
      Opt_tariff_all median:          1.048 (0.3% divergence)
      Within-pair MAD:                0.104 (Soderbery target: 0.125)
      Convergence = 0 share:          99% of attempted cells
      Tier 1 share:                   70.9%
      Tier 3 share:                   26.1%

  STRUCTURAL RATIOS:

                             Country        Regional       Soderbery
    gamma/(1+gamma):          0.550          0.518           0.408
    1/(sigma-1):              0.333          0.338           0.532
    gamma/((1+g)(s-1)):       0.168          0.167           0.217

  Deviations from Soderbery's ratios are attributable to:
    - Bounded upward bias in Feenstra sigma (drives 1/(sigma-1) low)
    - Sample composition differences (1995-2024 BACI vs 1991-2007
      Comtrade; different country coverage and sector weighting)
    Both deviations are documented features, not bugs.

  BEFORE-AFTER COMPARISON (country-level, HS4):

                             Pre-patch           Current
    Gamma median:              877.014             1.223
    Gamma IQR:                 [0.96, 11873]       [0.71, 1.99]
    Gamma > 20 share:           62.8%               0.9%
    Sigma median:               27.320              4.007

  The pre-patch pipeline's country-level output was a failure mode of
  the joint estimator on BACI country-level data, not a valid
  estimation. The three-stage pipeline with fixed sigma and moderate
  shrinkage brings all metrics into Soderbery's neighborhood.


INTERPRETING VARIANCE DECOMPOSITION R-SQUARED
==============================================
  Under the headline specification (exporter_weight = "trade_value",
  Stage 2a lambda = 0.05, Stage 2b lambda = 0.1), the full-sample
  R-squared (importer x product FEs on log gamma) varies with the
  shrinkage level. The lambda sweep at Stage 2a produced:

    lambda   plateau>20   MAD       R-squared
    0.01      1.73%       0.188      0.320
    0.02      0.96%       0.107      0.428
    0.05      0.48%       0.047      0.593    (headline)
    0.10      0.33%       0.024      0.716

  The Soderbery targets (MAD = 0.125, R-squared = 0.72) trace out a
  trade-off curve rather than a single point on BACI HS4 data. Higher
  lambda compresses within-cell variation (lower MAD) while letting
  cross-cell means be driven by the product-level prior (higher R²).
  Selection at lambda = 0.05 balances the two, with plateau share
  well under 1%.

  Under the uniform weighting robustness specification, R-squared
  values collapse substantially (approximately 0.13 regional, 0.06
  country). The gap between weighting schemes reflects how heavily
  the objective weights noise-dominated moment contributions from
  small bilateral relationships; trade-value weighting lets the
  large, persistent relationships dominate identification as in
  Soderbery's intuition (p. 50, fn 14).

  For methodological appendices:
    - Report the headline (lambda=0.05 regional / lambda=0.1 country,
      trade_value weighting)
    - Report the uniform-weighted run as a robustness comparison
    - Report the Stage 2a lambda sweep table


CONFIGURATION OPTIONS
=====================
  filepath                Path to BACI data (directory or single file)
  minyear                 First year to include (default: 1995 for HS92)
  maxyear                 Last year to include (default: NULL = all available)
  agg_level               "hs4" (paper default) or "hs6"
  custom_region_map       Optional data.table(cty_code, region) to
                          override the built-in UN M49 mapping
  min_exporters           Min non-reference exporters per cell (default: 2)
  min_destinations        Min export destinations per exporter (default: 2)
  min_periods             Min years of differenced data (default: 3)
  uv_outlier_threshold    |delta ln(unit value)| cutoff (default: 2.0)
  sigma_start             Starting value for Stage 1 optimizer (default: 2.88)
  gamma_start             Starting value for Stage 1 optimizer (default: 0.69)
  sigma_V_default         Default reference destination sigma.
                          Updated to Stage 1 median for Stages 2a/2b.
  gamma_V_default         Default reference destination gamma.
                          Updated to Stage 1 median for Stages 2a/2b.
  shrinkage_lambda        Penalty strength in Stages 2a and 2b.
                          Stage 2a default: 0.05 (HS4-calibrated).
                          Stage 2b default: 0.1.
                          Higher lambda = stronger pull toward prior.
                          See LAMBDA CALIBRATION above for sweep results.
                          At HS6, re-run the sweep to pick new lambdas —
                          moment magnitudes and plateau behavior shift
                          with sample composition.
  exporter_weight         Across-exporter WLS weighting mode.
                          "trade_value" (default): weight each
                          exporter's moment row by cell trade value,
                          normalized to sum to J. Large bilateral
                          relationships contribute more to the
                          objective. Applied symmetrically to import
                          and export side moments.
                          "uniform": rep(1, J). Every exporter gets
                          equal weight regardless of size. Restores
                          the baseline Soderbery specification;
                          useful as a robustness check.
  weight_period_floor     Under trade_value weighting, exporters with
                          fewer than this many periods have their
                          weight scaled by n_periods / weight_period_floor.
                          Protects against single-year mega-shipments
                          dominating the objective via raw trade value.
                          Default: 10L. At HS6 consider 5-7 for thinner
                          bilateral panels.
  tier1_min_periods       Min periods for Tier 1 estimation (default: 3L,
                          matching Soderbery minimum; effectively all
                          qualifying exporters go to Tier 1)
  tier1_min_dests         Min destinations for Tier 1 (default: 2L)
  tier2_min_periods       Min periods for Tier 2 (default: 3L)
  tail_trim_pct           Fraction trimmed from each tail of both sigma
                          and gamma distributions (default: 0.005 = 0.5%).
                          Trim bounds computed from Tier 0/1/2 cells only
                          (excludes Tier 3 imputations which cluster at
                          the regional prior).

  INTERNAL / POPULATED BY THE PIPELINE (users typically don't touch):
  sigma_lookup            (importer, good) -> sigma from Stage 1
  sigma_V_lookup          (ref_dest, good) -> sigma for export-side
                          moments (bilateral reference elasticities)
  gamma_V_lookup          (ref_dest, good) -> gamma for export-side
                          moments at Stage 2b (maps country ref_dest
                          to its region's Stage 2a gamma)
  shrinkage_priors        (good) -> ln_gamma_prior for the shrinkage
                          penalty. Feenstra-derived at Stage 2a,
                          Stage-2a-derived at Stage 2b.
  regional_starts         (region, good) -> (sigma, gamma) for
                          per-cell optimizer starting values at Stage 2b


METHODOLOGICAL NOTES
====================
  1. STAGE 1 — FEENSTRA (1994) SIGMA
     Estimates sigma per (importer, product) from import-side moments
     under the homogeneity assumption gamma_j = gamma_k = gamma_common.
     Two-parameter optimization over (sigma, gamma_common). Well-
     identified because multiple exporters' hyperbolae must intersect
     at a single point. Sigma is used as a fixed input for Stages 2a
     and 2b. Gamma_common is used as the Stage 2a shrinkage anchor
     and as plateau-fallback value for cells that escape bounds.

     Known issue: Under heterogeneous gamma, the Feenstra sigma has
     a bounded upward bias (Soderbery 2015, Section 4). The bias is
     accepted as the cost of breaking the identification ridge.
     Empirically, HS4 BACI 1995-2024 sigma median is 3.88 vs
     Soderbery's 2.88, consistent with the predicted direction.

  2. STAGE 2a — REGIONAL GAMMA WITH MODERATE SHRINKAGE
     Estimates heterogeneous gamma for Soderbery's 20 regional
     entities with sigma fixed from Stage 1. A shrinkage penalty
     (lambda = 0.05 at HS4) pulls gamma toward the Stage 1
     gamma_common aggregated to product-level median. This prevents
     the gamma plateau while letting well-identified regional cells
     deviate from the prior.

     Important empirical finding: at lambda = 0 (fully unshrunk),
     51% of Stage 2a regional cells land above gamma = 20 on BACI
     HS4 data. The 20-entity regional panel is NOT sufficiently
     well-identified to go fully unshrunk on 1995-2024 BACI data,
     contrary to the README's original premise. Moderate shrinkage
     at lambda = 0.05 reduces plateau share to ~0.5% while
     preserving within-cell heterogeneity.

     Regional results serve two purposes: they are publication
     outputs in their own right (replicating Soderbery's framework
     on BACI data), and they provide the product-level shrinkage
     priors for Stage 2b.

     A plateau-fallback mechanism replaces any residual outlier cells
     (gamma > 20) with the Feenstra anchor. Affected cells are
     flagged with convergence code -2; opt_tariff is recomputed to
     exclude their contributions. At HS4 with lambda = 0.05, only
     2,081 cells (0.47%) trigger this mechanism.

  3. STAGE 2b — COUNTRY GAMMA WITH SHRINKAGE
     Estimates heterogeneous gamma for individual countries with
     sigma fixed from Stage 1 and shrinkage (lambda = 0.1 at HS4)
     toward the Stage 2a regional priors. Dense bilateral cells have
     strong moments that dominate the penalty; thin cells are pulled
     toward the regional mean for their product.

     Thin cells with insufficient identification to produce reliable
     estimates are assigned the regional prior directly as Tier 3
     (convergence code -1). Tier 3 rows are excluded from the
     primary opt_tariff calculation; opt_tariff_all retains them for
     full-coverage analyses.

  4. SHRINKAGE PENALTY
     lambda * sum_j (ln(gamma_j) - ln_gamma_prior)^2, applied to
     each exporter's gamma in the cell. The log-space formulation
     penalizes multiplicative deviations symmetrically (gamma = 10
     is as far from gamma = 1 as gamma = 0.1 is). With
     ln_gamma_prior varying by product, the penalty respects
     product-level structure while allowing within-product
     heterogeneity.

     LAMBDA CALIBRATION. Stage 2a lambda = 0.05 and Stage 2b
     lambda = 0.1 were chosen on HS4 empirical diagnostics.
     Stage 2a lambda was swept across {0.01, 0.02, 0.05, 0.1} and
     selected to minimize plateau share while preserving within-cell
     heterogeneity. Stage 2b lambda was not swept; its default of
     0.1 produces within-pair MAD = 0.104 at country level, close
     to Soderbery's 0.125 benchmark.

     A Stage 2b lambda sweep is recommended for appendix-level
     robustness and is enabled by the lambda_calibration_diagnostic
     and lambda_calibration_sweep helper functions.

     AT HS6: Lambda calibration should be re-run. Product count
     rises from ~1,240 to ~5,300; within-product panels become
     thinner; plateau behavior and MAD will shift. The existing
     lambda values are HS4-specific.

  5. TRADE-VALUE WEIGHTING OF EXPORTER MOMENTS
     BW weights (Soderbery Eq. 14, Broda-Weinstein 2006) operate
     WITHIN an exporter's time series to handle measurement error
     in per-year unit values. They do not weight across exporters.
     The baseline Soderbery specification gives every exporter row
     equal weight in the across-exporter sum of squared residuals,
     treating a tiny bilateral flow with 3 observations identically
     to a dense, persistent relationship with 30 observations.

     Under exporter_weight = "trade_value", each exporter's moment
     row is weighted by its total trade value within the cell,
     normalized to sum to J. Large and persistent relationships
     thereby dominate identification, consistent with Soderbery's
     intuition (p. 50, fn 14) applied across rather than within
     exporters. The same weighting scheme is applied to import-side
     and export-side moments to preserve coherence of the joint
     objective.

     WEIGHT PERIOD FLOOR. To protect against lumpy single-year
     shipments (aircraft, LNG cargoes, rough diamonds) receiving
     outsized weights via raw trade value despite providing little
     identifying information, exporters with fewer than
     weight_period_floor periods (default 10) have their trade-
     value weight scaled by n_periods / weight_period_floor.

     Within-pair heterogeneity statistics (gamma MAD, SD) decline
     under trade-value weighting because fewer exporters are
     effectively "speaking" in each cell. This is a deliberate
     tradeoff: the weighting concentrates identification on
     relationships with sufficient data to support it. Report both
     specifications in the appendix to make the tradeoff transparent.

  6. REFERENCE-DESTINATION ELASTICITY LOOKUPS
     Soderbery's export-side Eq. (11) depends on (sigma_V, gamma_V)
     at each Tier 1 exporter's reference destination V. The baseline
     treatment uses a single global default for all exporters; this
     code uses per-(ref_dest, good) lookups with fallback to defaults.

     At Stage 2a (regional):
       sigma_V_lookup: regional-aggregated Stage 1 sigma indexed by
                       (region, good). ref_dest is a region.
       gamma_V_lookup: NULL (regional gamma hasn't been estimated
                       yet at this stage). Falls back to default.

     At Stage 2b (country):
       sigma_V_lookup: country-level Stage 1 sigma indexed by
                       (country, good). ref_dest is a country.
       gamma_V_lookup: Stage 2a regional gamma mapped to country
                       level via the build_region_map(). Each
                       country inherits its region's median gamma
                       for that product. Coarse but structurally
                       consistent with the regional priors already
                       in use for shrinkage.

     The bilateral lookup is a structural refinement over the
     Soderbery baseline, exploiting BACI's country-level granularity
     while maintaining compatibility when no match exists.

  7. TIER CLASSIFICATION
     With default thresholds (3 periods, 2 destinations, matching
     Soderbery's minimum data requirements), nearly all qualifying
     exporters enter Tier 1 — full import + export side estimation.
     The tier framework was originally designed to manage compute
     costs on small-core machines; on 62+ core infrastructure it
     is largely vestigial at HS4 (Tier 1 ~89% of regional cells,
     71% of country cells; Tier 3 26% of country cells driven by
     thin bilateral relationships rather than thresholds).

     The tier column is retained in the output for downstream
     analysis (e.g., running robustness checks on Tier 1 estimates
     only), for computing opt_tariff from directly-estimated cells
     only, and for backward compatibility with tighter thresholds
     if compute is constrained.

  8. REGIONAL AGGREGATION
     The build_region_map() function maps Comtrade numeric country
     codes to the 20 entities in Soderbery Table 1 (13 individual
     countries + 7 UN M49-based regions). Users should verify the
     mapping against their BACI country_codes metadata file.

  9. UNIT VALUE OUTLIER FILTER
     Observations with |delta ln(unit value)| >= 2.0 are dropped
     after first-differencing. Standard in the Feenstra/Broda-
     Weinstein/Soderbery tradition; removes measurement error from
     quantity misreporting.

 10. POST-ESTIMATION TRIMMING
     A symmetric percentile trim (default: 0.5% each tail) is
     applied to both sigma and gamma distributions. Trim bounds
     are computed from Tier 0/1/2 cells only (excludes Tier 3
     imputations which cluster at the regional prior and would
     bias the quantile cuts). Bounds are then applied to all rows.

 11. DATA CACHING
     Raw BACI data (through HS4 aggregation) is loaded once and
     cached to {prefix}_raw_cache.rds. All three stages reuse
     this cache rather than re-loading and re-aggregating 110M+
     rows. The country-level and regional finalizations (shares,
     first-differences, regional grouping) are computed from the
     cache and held in memory across stages, not re-computed per
     stage.

 12. DATA QUALITY DIAGNOSTICS
     A quality log tracks observations at each stage: raw load,
     positive-value filter, year filter, HS4 aggregation, regional
     aggregation, first-differencing, outlier filter, estimation
     outcomes, and post-estimation trimming. A formatted report
     prints at completion and is embedded in the summary outputs.


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
    Soderbery's regions are composites. The following assignments
    depart from strict M49 sub-region classification but are
    consistent with the paper's Table 6 statistics and economic
    logic of trade patterns:

    - Spain (724), Portugal (620), Greece (300), Malta (470):
      Assigned to NWU (Northern/Western Europe), not M49 Southern
      Europe. Established EU members with the EU common external
      tariff.

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
  STAGE TIMINGS (measured on AWS c7a.16xlarge, 62 cores, 1995-2024
  BACI HS92, HS4):

    Raw data load and cache:         5-10 min (once; reused across stages)
    Stage 1 (Feenstra sigma):        ~3 min on 62 cores
    Stage 2a (regional gamma):       ~3 min on 62 cores
    Stage 2b (country gamma):        ~35 min on 62 cores
    Summary writing:                 ~1 min

  TOTAL PIPELINE (from fresh raw data):
    Large machine (62 cores, 128GB+): ~45 min
    Small machine (4 cores, 32GB):    5-8 hours (estimated)

  AT HS6 (EXTRAPOLATED):
    Product count rises from ~1,240 to ~5,300 (4.3x).
    Raw data load:                   similar (same underlying rows)
    Stage 1:                         ~10-15 min on 62 cores
    Stage 2a:                        ~10 min on 62 cores
    Stage 2b:                        ~2-3 hours on 62 cores
    Total:                           ~3-4 hours end-to-end

  STAGE DETECTION:
    Every stage detects previously completed output via its RDS file.
    If present, the stage is skipped and results are loaded from disk.
    Delete the file to force re-estimation. Deleting a stage's output
    also requires re-running any downstream stages that depend on it.

  CHECKPOINT AND RESUME:
    All stages support per-batch checkpointing on both the forked
    (Linux/Mac mclapply) and socket (Windows parLapply) parallel
    paths, as well as per-product checkpointing on the serial path.
    Checkpoints save to {prefix}_feenstra_checkpoint.rds (Stage 1) or
    {prefix}_fs_checkpoint.rds (Stage 2a/2b) and are automatically
    removed after successful completion. A crashed run resumes from
    the last completed batch without repeating completed stages.

  SPEED OPTIMIZATIONS:

    - Rcpp objective functions: Both het_obj_rcpp.cpp (joint) and
      het_obj_fixed_sigma_rcpp.cpp (Stage 2) are compiled at load
      time. Each parallel worker compiles independently. Requires
      Rcpp package and Rtools (Windows) or Xcode command line tools
      (Mac). Falls back to pure R silently if unavailable.

    - Data caching: Raw data is loaded once through HS4 aggregation
      and reused across stages. Country-level and regional finalizations
      (shares, first-differences) are held in memory and passed to
      stages directly, avoiding re-preparation.

    - Per-exporter lookup: Product-level data is pre-split by exporter
      once per product (compute_exporter_lookup), converting the
      O(N_exporters) repeated filtering in build_export_moments into
      O(1) list index lookups. Material speedup at HS4; larger at HS6.

    - Regional-informed starting values: Stage 2b automatically uses
      Stage 2a regional estimates as per-cell starting values. Instead
      of starting every cell at global defaults, each cell starts at
      its region-product median. Reduces optimizer iterations 30-50%.

    - Pre-filtering: Before calling the estimation function for each
      importer, the code checks minimum data requirements. Cells that
      cannot meet them are skipped without entering the function.

  WINDOWS PARALLEL NOTES:
    Windows uses socket clusters (makeCluster/parLapply), which
    require data serialization to worker processes. To avoid memory
    failures with large datasets, the code writes product slices to
    temporary RDS files and has workers read from disk rather than
    receiving data through the socket.

    On Linux/Mac, forked parallelism (mclapply) shares memory and
    avoids this overhead entirely.

  OTHER NOTES:
    - Progress reporting: Serial path reports every 10 products
      with elapsed time and ETA. Parallel path reports per batch.
    - The optimizer uses L-BFGS-B (bounded, gradient-based) with a
      Nelder-Mead fallback (maxit=1000) for cells where L-BFGS-B
      fails to converge.
    - Cell failure diagnostics: When cells fail to produce estimates,
      the reason is logged (insufficient_data, no_nonref_exporters,
      no_valid_moments, no_sigma_estimate, all_tier3_no_prior,
      optimizer_failed, or error with message). Counts by reason
      are in the summary.
    - The convergence column flags non-converged cells
      (convergence != 0). Convergence = -1 specifically indicates
      Tier 3 assignment (gamma set from regional prior, not
      estimated); -2 indicates Stage 2a plateau-fallback replacement.
      Consider separating these from analyses of directly estimated
      gamma, or use opt_tariff (Tier 0/1/2 only) rather than
      opt_tariff_all.
    - Config validation runs before data loading to catch obvious
      misconfigurations before expensive operations.


CALIBRATION AND ROBUSTNESS TOOLING
===================================
  The library exposes three helper functions for empirical
  calibration of shrinkage lambda on new samples (e.g., HS6):

    lambda_calibration_diagnostic(results, lambda)
      Returns within-pair MAD, variance-decomposition R², plateau
      share, and Tier 1/3 distinctness for a single run.

    print_lambda_diagnostic(diag)
      Formatted console output.

    lambda_calibration_sweep(cfg_base, prepared_dt, lambda_grid, ncores)
      One-shot driver that runs Stage 2 at each lambda in the grid
      and returns a combined diagnostic table.

  Recommended before committing to a full HS6 run:
    1. Complete Stage 1 + raw cache at HS6
    2. Run lambda_calibration_sweep at regional level across
       {0.01, 0.02, 0.05, 0.1}
    3. Pick Stage 2a lambda minimizing plateau share while
       preserving MAD
    4. Optionally sweep Stage 2b lambda on a product subset
    5. Run full Stage 2b with chosen lambdas


KNOWN LIMITATIONS / FUTURE WORK
================================
  1. STAGE 2b LAMBDA NOT SWEPT AT HS4. The current Stage 2b lambda
     (0.1) produces within-pair MAD of 0.104 (close to Soderbery's
     0.125) but was not formally calibrated against a sweep.
     Appendix-level robustness should report a Stage 2b lambda
     sweep alongside the Stage 2a sweep.

  2. TWO SODERBERY TARGETS ARE JOINTLY INCONSISTENT ON BACI.
     MAD = 0.125 and R-squared = 0.72 trace a trade-off curve
     along lambda rather than a single point. This is a sample
     property, not a methodology failure; document transparently
     in any publication.

  3. REFERENCE DESTINATION TIE-BREAKING. build_export_moments picks
     the reference destination as the largest non-focal destination
     by trade value. At HS6 ties will be more common and this choice
     is less robust. A low-CV-based tie-breaking alternative is
     available via the robustness script.

  4. AT HS6, LAMBDA MUST BE RECALIBRATED. The HS4-selected lambdas
     (0.05 regional, 0.1 country) are sample-specific. Rerun the
     lambda sweep at HS6 before committing to full estimation.

  5. STAGE 2a PLATEAU AT LAMBDA = 0. Empirically, 51% of regional
     cells hit the plateau without shrinkage on BACI HS4 data,
     contradicting the original assumption that the 20-entity
     regional panel would be well-identified enough to go unshrunk.
     A version of the pipeline with plateau fallback as the primary
     regularization mechanism (rather than shrinkage plus fallback)
     was tested and produced a bimodal output distribution where
     the "fallback" half was the Feenstra anchor; the current
     design (moderate shrinkage plus residual fallback) is cleaner.


CITATION
========
  Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
  Tariffs," Journal of International Economics, 114, 2018, pp. 44-62.

  Feenstra, Robert C., "New Product Varieties and the Measurement
  of International Prices," American Economic Review, 84(1), 1994,
  pp. 157-177.

  Soderbery, Anson, "Estimating import supply and demand elasticities:
  Analysis and implications," Journal of International Economics, 96(1),
  2015, pp. 1-17.

  Broda, Christian, and David E. Weinstein, "Globalization and the
  Gains from Variety," Quarterly Journal of Economics, 121(2), 2006,
  pp. 541-585.
