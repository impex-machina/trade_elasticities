// =========================================================================
// het_obj_fixed_sigma_rcpp.cpp
//
// Rcpp objective function for gamma-only estimation with sigma fixed.
// Used in Stage 2 of the two-stage estimator where sigma is pre-estimated
// via Feenstra (1994) and held constant.
//
// The parameter vector d contains ONLY gamma values:
//   d[0] = gamma_k (reference exporter)
//   d[1:] = gamma_1, ..., gamma_J (non-reference exporters)
//
// Sigma is passed as a separate fixed scalar.
//
// Optionally adds a hierarchical shrinkage penalty:
//   lambda * sum_j (ln(gamma_j) - ln_gamma_prior)^2
// which pulls estimates toward a product-level prior (e.g., from
// regional estimation).
//
// CITATION:
//   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
//   Tariffs," JIE, 114, 2018, pp. 44-62.
//
// Last updated: 2026-04-10
// =========================================================================

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double het_obj_fixed_sigma_rcpp(NumericVector d,
                                double sigma,
                                NumericVector imp_Y, NumericMatrix imp_X,
                                NumericVector exp_Y, NumericMatrix exp_X,
                                IntegerVector exp_jmap,
                                NumericVector exp_sig_V, NumericVector exp_gam_V,
                                NumericVector wt_imp, NumericVector wt_exp,
                                double ln_gamma_prior,
                                double shrinkage_lambda) {

  // d[0] = gamma_k, d[1:J] = gamma_j
  // sigma is FIXED (not part of d)
  double gam_k = d[0];
  int J = d.size() - 1;

  double sig = sigma;

  // --- Enforce constraints ---
  if (sig <= 1.0 || gam_k <= 0.0) return 1e12;
  for (int j = 0; j < J; j++) {
    if (d[j + 1] <= 0.0) return 1e12;
  }

  double sm1 = sig - 1.0;

  // ===========================================================
  //  IMPORT-SIDE RESIDUALS — Eq (10)
  //
  //  Identical to het_obj_rcpp but sigma is fixed externally.
  //  Index shift: gamma_j = d[j+1] instead of d[j+2]
  // ===========================================================

  double SSR_imp = 0.0;

  for (int j = 0; j < J; j++) {
    double gam_j = d[j + 1];
    double inv_1pgj = 1.0 / (1.0 + gam_j);

    double pred = (gam_j * inv_1pgj / sm1)                             * imp_X(j, 0) +
                  (gam_j * inv_1pgj)                                    * imp_X(j, 1) +
                  (-1.0 / sm1)                                          * imp_X(j, 2) +
                  (gam_j * gam_k * inv_1pgj / ((1.0 + gam_k) * sm1))   * imp_X(j, 3) +
                  ((gam_j - gam_k) * inv_1pgj / gam_k)                 * imp_X(j, 4);

    double resid = imp_Y[j] - pred;
    SSR_imp += wt_imp[j] * resid * resid;
  }

  // ===========================================================
  //  EXPORT-SIDE RESIDUALS — Eq (11)
  //
  //  exp_jmap uses original 1-based indexing into the FULL d
  //  vector (where d[1]=sigma, d[2]=gam_k, d[3:]=gam_j).
  //  Since sigma is removed, we subtract 1 from the index:
  //  d_fixed[exp_jmap[m] - 2] = gamma_j for that exporter
  //  (C++ 0-based: exp_jmap[m] - 1 was original, now - 2)
  // ===========================================================

  int N_exp = exp_Y.size();
  double SSR_exp = 0.0;

  if (N_exp > 0) {
    for (int m = 0; m < N_exp; m++) {
      // Original: d[exp_jmap[m] - 1] with d = (sigma, gam_k, gam_j...)
      // Fixed:    d[exp_jmap[m] - 2] with d = (gam_k, gam_j...)
      double gam_I = d[exp_jmap[m] - 2];
      double gam_V = exp_gam_V[m];
      double sig_V = exp_sig_V[m];

      double inv_1pgI = 1.0 / (1.0 + gam_I);
      double inv_1pgV = 1.0 / (1.0 + gam_V);

      double pred = (gam_I * inv_1pgI / sm1)                                         * exp_X(m, 0) +
                    ((gam_I * (sig - 2.0) - 1.0) * inv_1pgI / sm1)                   * exp_X(m, 1) +
                    (gam_V * inv_1pgV / sm1)                                          * exp_X(m, 2) +
                    ((1.0 - gam_V * (sig - 2.0)) * inv_1pgV / sm1)                   * exp_X(m, 3) +
                    ((1.0 - gam_V * (sig_V - 2.0)) * inv_1pgV / sm1)                 * exp_X(m, 4) +
                    ((gam_I * (sig_V - 2.0) - 1.0) * inv_1pgI / sm1)                 * exp_X(m, 5) +
                    (-(gam_V * (1.0 + gam_I) + gam_I * (1.0 + gam_V)) *
                       inv_1pgI * inv_1pgV / sm1)                                     * exp_X(m, 6) +
                    ((sig - sig_V) / sm1)                                             * exp_X(m, 7) +
                    ((sig_V - sig) / sm1)                                             * exp_X(m, 8);

      double resid = exp_Y[m] - pred;
      SSR_exp += wt_exp[m] * resid * resid;
    }
  }

  // ===========================================================
  //  SHRINKAGE PENALTY
  //
  //  lambda * sum_j (ln(gamma_j) - ln_gamma_prior)^2
  //
  //  Applied to all gamma parameters (both ref and non-ref)
  //  that are above the boundary (1e-5). Pulls estimates
  //  toward the product-level prior from regional estimation.
  // ===========================================================

  double penalty = 0.0;

  if (shrinkage_lambda > 0.0 && !std::isnan(ln_gamma_prior)) {
    for (int i = 0; i < d.size(); i++) {
      if (d[i] > 1e-5) {
        double dev = std::log(d[i]) - ln_gamma_prior;
        penalty += dev * dev;
      }
    }
    penalty *= shrinkage_lambda;
  }

  return SSR_imp + SSR_exp + penalty;
}
