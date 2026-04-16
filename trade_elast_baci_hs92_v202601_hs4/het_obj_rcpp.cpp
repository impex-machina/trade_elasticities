// =========================================================================
// het_obj_rcpp.cpp
//
// Rcpp implementation of the joint nonlinear SUR objective function
// for Soderbery (2018) Eqs. (10) and (11).
//
// Drop-in replacement for het_obj() in het_obj.R.
// Compile via Rcpp::sourceCpp("het_obj_rcpp.cpp") or automatically
// from feen94_het_baci.R.
//
// CITATION:
//   Soderbery, Anson, "Trade Elasticities, Heterogeneity, and Optimal
//   Tariffs," JIE, 114, 2018, pp. 44-62.
//
// Last updated: 2026-03-28
// =========================================================================

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double het_obj_rcpp(NumericVector d,
                    NumericVector imp_Y, NumericMatrix imp_X,
                    NumericVector exp_Y, NumericMatrix exp_X,
                    IntegerVector exp_jmap,
                    NumericVector exp_sig_V, NumericVector exp_gam_V,
                    NumericVector wt_imp, NumericVector wt_exp) {

  double sig   = d[0];
  double gam_k = d[1];
  int J = d.size() - 2;

  // --- Enforce constraints ---
  if (sig <= 1.0 || gam_k <= 0.0) return 1e12;
  for (int j = 0; j < J; j++) {
    if (d[j + 2] <= 0.0) return 1e12;
  }

  double sm1 = sig - 1.0;

  // ===========================================================
  //  IMPORT-SIDE RESIDUALS — Eq (10)
  // ===========================================================

  double SSR_imp = 0.0;

  for (int j = 0; j < J; j++) {
    double gam_j = d[j + 2];
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
  // ===========================================================

  int N_exp = exp_Y.size();

  if (N_exp == 0) return SSR_imp;

  double SSR_exp = 0.0;

  for (int m = 0; m < N_exp; m++) {
    // exp_jmap uses 1-based R indexing into d
    double gam_I = d[exp_jmap[m] - 1];
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

  return SSR_imp + SSR_exp;
}
