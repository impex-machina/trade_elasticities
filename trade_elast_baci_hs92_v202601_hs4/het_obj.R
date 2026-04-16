#' ===========================================================================
#' het_obj.R
#'
#' Joint nonlinear SUR objective function for Soderbery (2018)
#' "Trade Elasticities, Heterogeneity, and Optimal Tariffs"
#' JIE 114, pp. 44-62.
#'
#' Vectorized: all residuals computed via element-wise vector operations
#' rather than row-by-row for-loops. Mathematically identical output.
#'
#' Last updated: 2026-03-27
#' ===========================================================================


#' Joint nonlinear SUR objective for eqs (10) and (11)
#'
#' @param d Numeric vector of parameters:
#'   d[1]     = sigma (common within import market)
#'   d[2]     = gamma_k (reference exporter)
#'   d[3:end] = gamma_1, ..., gamma_J (non-reference exporters)
#'   Constraints: sigma > 1, all gamma > 0.
#'
#' @param imp_Y Numeric vector (J): import-side LHS, time-averaged per exporter.
#' @param imp_X Numeric matrix (J x 5): import-side RHS.
#'   Cols: (Dk ls)^2, (Dk ls)(Dk lp), (Dk ls)(D lp_j), (Dk ls)(D lp_k), (Dk lp)(D lp_k)
#'
#' @param exp_Y Numeric vector (M): export-side LHS. Length 0 if unavailable.
#' @param exp_X Numeric matrix (M x 9): export-side RHS. 0 rows if unavailable.
#'   Cols: (D ls_I)^2, (D ls_I)(D lp_I), (D ls_V)^2, (D ls_V)(D lp_I),
#'         (D ls_V)(D lp_V), (D ls_I)(D lp_V), (D ls_I)(D ls_V),
#'         (D lp_V)^2, (D lp_V)(D lp_I)
#'
#' @param exp_jmap Integer vector (M): maps each export obs to its index in d.
#' @param exp_sig_V Numeric vector (M): sigma at each reference destination.
#' @param exp_gam_V Numeric vector (M): gamma at each reference destination.
#' @param wt_imp Numeric vector (J): import-side weights.
#' @param wt_exp Numeric vector (M): export-side weights.
#'
#' @return Scalar: weighted sum of squared residuals.
#'
het_obj <- function(d, imp_Y, imp_X, exp_Y, exp_X,
                    exp_jmap, exp_sig_V, exp_gam_V,
                    wt_imp, wt_exp) {

  sig   <- d[1]
  gam_k <- d[2]
  J     <- length(d) - 2L

  # ---- Enforce constraints ----
  if (sig <= 1 || gam_k <= 0) return(1e12)
  if (J > 0L && any(d[3:(J + 2L)] <= 0)) return(1e12)

  sm1 <- sig - 1

  # ===========================================================
  #  IMPORT-SIDE RESIDUALS — Eq (10), vectorized
  #
  #  gam_j is a vector (length J), one per non-ref exporter.
  #  Coefficient vectors are built once, then multiplied
  #  element-wise against the corresponding imp_X column.
  # ===========================================================

  gam_j    <- d[3:(J + 2L)]
  inv_1pgj <- 1 / (1 + gam_j)

  pred_imp <- (gam_j * inv_1pgj / sm1)                         * imp_X[, 1] +
              (gam_j * inv_1pgj)                                * imp_X[, 2] +
              (-1 / sm1)                                        * imp_X[, 3] +
              (gam_j * gam_k * inv_1pgj / ((1 + gam_k) * sm1)) * imp_X[, 4] +
              ((gam_j - gam_k) * inv_1pgj / gam_k)             * imp_X[, 5]

  SSR_imp <- sum(wt_imp * (imp_Y - pred_imp)^2)

  # ===========================================================
  #  EXPORT-SIDE RESIDUALS — Eq (11), vectorized
  #
  #  gam_I is looked up per obs via exp_jmap.
  #  gam_V, sig_V vary per obs (from reference destinations).
  #  All coefficient computations are vector operations.
  # ===========================================================

  N_exp <- length(exp_Y)

  if (N_exp > 0L) {
    gam_I <- d[exp_jmap]
    gam_V <- exp_gam_V
    sig_V <- exp_sig_V

    inv_1pgI <- 1 / (1 + gam_I)
    inv_1pgV <- 1 / (1 + gam_V)

    pred_exp <- (gam_I * inv_1pgI / sm1)                                     * exp_X[, 1] +
                ((gam_I * (sig - 2) - 1) * inv_1pgI / sm1)                   * exp_X[, 2] +
                (gam_V * inv_1pgV / sm1)                                      * exp_X[, 3] +
                ((1 - gam_V * (sig - 2)) * inv_1pgV / sm1)                   * exp_X[, 4] +
                ((1 - gam_V * (sig_V - 2)) * inv_1pgV / sm1)                * exp_X[, 5] +
                ((gam_I * (sig_V - 2) - 1) * inv_1pgI / sm1)                * exp_X[, 6] +
                (-(gam_V * (1 + gam_I) + gam_I * (1 + gam_V)) *
                   inv_1pgI * inv_1pgV / sm1)                                 * exp_X[, 7] +
                ((sig - sig_V) / sm1)                                         * exp_X[, 8] +
                ((sig_V - sig) / sm1)                                         * exp_X[, 9]

    SSR_exp <- sum(wt_exp * (exp_Y - pred_exp)^2)
    return(SSR_imp + SSR_exp)
  } else {
    return(SSR_imp)
  }
}
