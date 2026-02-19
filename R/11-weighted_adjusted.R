#' Weighted Adjusted Fractions
#'
#' Calculates the weighted adjusted potential impact fractions (or population
#' attributable fractions). Each individual fraction \eqn{\widehat{\text{PIF}}_i}
#' is rescaled proportionally so that together they are consistent with
#' the ensemble fraction.
#'
#' @param pif1 A potential impact fraction (class `pif_class`). This is
#'   the first of the individual fractions to be adjusted.
#'
#' @param ... The remaining potential impact fractions (class `pif_class`).
#'   All fractions must be of the same type (all `PIF` or all `PAF`).
#'
#' @param pif_total_link Link to pass to `pif_total` in the denominator
#' of the adjustment
#'
#' @param pif_total_link_inv Inverse of the link to pass to `pif_total` in the denominator
#' of the adjustment
#'
#' @param pif_total_link_deriv Derivative of the link to pass to `pif_total` in the denominator
#' of the adjustment
#'
#' @param pif_ensemble_link Link to pass to `pif_ensemble` in the numerator
#' of the adjustment
#'
#' @param pif_ensemble_link_inv Inverse of the link to pass to `pif_ensemble` in the numerator
#' of the adjustment
#'
#' @param pif_ensemble_link_deriv Derivative of the link to pass to `pif_ensemble` in the numerator
#' of the adjustment
#'
#' @param weights Weights for the ensemble (`pif_ensemble`).
#'   Passed directly to [pif_ensemble()] / [paf_ensemble()].
#'   Defaults to `NULL` (equal weights of 1 for each fraction).
#'
#' @param var_weights Covariance structure for `weights`.
#'   Passed directly to [pif_ensemble()] / [paf_ensemble()].
#'   Defaults to `0`.
#'
#' @param var_pif_weights Covariance matrix between individual
#'   fractions and `weights`. Passed to [pif_ensemble()] /
#'   [paf_ensemble()]. Defaults to `NULL`.
#'
#' @param var_p Covariance matrix for the prevalence parameters `p` across
#'   all fractions. Passed to [cov_total_pif()] when computing cross-covariances.
#'   Defaults to `NULL`.
#'
#' @param var_beta Covariance matrix for the relative-risk parameters `beta`
#'   across all fractions. Passed to [cov_total_pif()] when computing
#'   cross-covariances. Defaults to `NULL`.
#'
#' @param label_ensemble Character label for the internally constructed
#'   ensemble fraction. Defaults to `NULL` (auto-generated).
#'
#' @param label_sum Character label for the internally constructed sum
#'   (total) fraction. Defaults to `NULL` (auto-generated).
#'
#' @inheritParams pifpaf
#' @inheritParams classes
#'
#' @return A named list of `pif_class` objects, one per input fraction,
#'   where each element is \eqn{\widehat{\text{PIF}}_i^{\text{adj}}}.
#'   Names are taken from the `label` slots of the input fractions.
#'
#' @section Formula:
#'
#' The weighted adjusted fraction for the \eqn{i}-th exposure is:
#' \deqn{
#'   \widehat{\text{PIF}}_i^{\text{adj}} =
#'     \frac{\widehat{\text{PIF}}_i}{\sum_{j=1}^{n} \widehat{\text{PIF}}_j}
#'     \cdot \widehat{\text{PIF}}_{\text{Ensemble}}
#' }
#'
#' Using the log-transform, the variance is:
#' \deqn{
#' \begin{aligned}
#'   \operatorname{Var}\!\Big[\ln \widehat{\text{PIF}}_i^{\text{adj}}\Big]
#'     &= \operatorname{Var}\!\left[\ln \widehat{\text{PIF}}_i\right]
#'      + \operatorname{Var}\!\left[\ln \textstyle\sum_j \widehat{\text{PIF}}_j\right]
#'      + \operatorname{Var}\!\left[\ln \widehat{\text{PIF}}_{\text{Ensemble}}\right] \\
#'     &\quad + 2\Bigg[
#'         \operatorname{Cov}\!\Big(\ln \widehat{\text{PIF}}_i,
#'           \ln \widehat{\text{PIF}}_{\text{Ensemble}}\Big)
#'       - \operatorname{Cov}\!\Big(\ln \widehat{\text{PIF}}_i,
#'           \ln \textstyle\sum_j \widehat{\text{PIF}}_j\Big)
#'       - \operatorname{Cov}\!\Big(\ln \widehat{\text{PIF}}_{\text{Ensemble}},
#'           \ln \textstyle\sum_j \widehat{\text{PIF}}_j\Big)
#'       \Bigg]
#' \end{aligned}
#' }
#'
#' where each log-covariance is approximated via the delta method:
#' \deqn{
#'   \operatorname{Cov}\!\big(\ln X, \ln Y\big) \approx
#'     \frac{\operatorname{Cov}(X, Y)}{X \cdot Y}
#' }
#'
#' The covariances \eqn{\operatorname{Cov}(\widehat{\text{PIF}}_i, \cdot)} are
#' computed by [cov_total_pif()] applied to the individual fraction and the
#' internally constructed sum/ensemble objects, which automatically propagates
#' the full uncertainty structure already embedded in those objects.
#'
#' @section Internal construction:
#'
#' Internally the function builds:
#' \itemize{
#'   \item A **sum** fraction via [pif_total_class()] with weights all equal to
#'     `1`, so its `pif` slot equals \eqn{\sum_i \widehat{\text{PIF}}_i}.
#'   \item An **ensemble** fraction via [pif_ensemble()] / [paf_ensemble()].
#' }
#' All cross-covariances are then computed automatically through the recursive
#' [cov_total_pif()] machinery that already handles atomic, total, and ensemble
#' fractions.
#'
#' @examples
#' \dontrun{
#' paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001, label = "Lead")
#' paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001, label = "Radiation")
#'
#' # Adjusted fractions (covariances computed automatically)
#' adj <- weighted_adjusted_fractions(paf_lead, paf_rad)
#' adj$Lead
#' adj$Radiation
#'
#' # With correlated prevalences
#' adj2 <- weighted_adjusted_fractions(paf_lead, paf_rad, var_p = 0.0001)
#' adj2$Lead
#' }
#' @seealso [pif_ensemble()], [paf_ensemble()], [cov_total_pif()],
#'   [weighted_adjusted_paf()], [weighted_adjusted_pif()]
#'
#' @keywords internal
weighted_adjusted_fractions <- function(
    pif1,
    ...,
    weights                  = NULL,
    pif_total_link           = "log-complement",
    pif_total_link_inv       = NULL,
    pif_total_link_deriv     = NULL,
    pif_ensemble_link        = "identity",
    pif_ensemble_link_inv    = NULL,
    pif_ensemble_link_deriv  = NULL,
    var_weights              = 0,
    var_pif_weights          = NULL,
    var_p                    = NULL,
    var_beta                 = NULL,
    conf_level               = 0.95,
    label_ensemble           = NULL,
    label_sum                = NULL,
    quiet                    = FALSE
) {

  # ── 1. Collect and validate inputs ──────────────────────────────────────────

  pif_list <- c(list(pif1), list(...))
  n        <- length(pif_list)

  if (is.null(weights)){
    weights <- rep(1, n)
  }

  if (length(weights) != n){
    cli::cli_abort(
      "{.val {n}} fractions where given but weights has length {length(weights)}. There should be a weight per fraction.")
  }

  if (!all(sapply(pif_list, function(x) S7::S7_inherits(x, pif_class)))) {
    cli::cli_abort(
      "All fractions passed to {.fn weighted_adjusted_fractions} must be {.cls pif_class}."
    )
  }

  types  <- unique(sapply(pif_list, function(x) x@type))
  is_paf <- length(types) == 1L && types == "PAF"


  # ── 3. Label management ─────────────────────────────────────────────────────

  labels <- sapply(pif_list, function(x) x@label)

  # Check for duplicate labels before we embed them into pif_list
  if (any(duplicated(labels))) {
    cli::cli_abort(
      "Duplicate labels found: {labels[duplicated(labels)]}. Use the `label` argument to assign unique labels."
    )
  }

  names(pif_list) <- labels

  sum_label <- if (is.null(label_sum)) {
    paste0("deltapif-sum-", sub("\\.", "", as.character(abs(stats::rnorm(1)))))
  } else {
    label_sum
  }

  ens_label <- if (is.null(label_ensemble)) {
    paste0("deltapif-ensemble-", sub("\\.", "", as.character(abs(stats::rnorm(1)))))
  } else {
    label_ensemble
  }

  # ── 4. Build the sum as pif_total_class with weights = rep(1, n) ───────────
  pif_sum <- pif_total(
    pif1, ...,
    weights         = rep(1, n),
    var_weights     = 0,
    var_pif_weights = NULL,
    link            = pif_total_link,
    link_inv        = pif_total_link_inv,
    link_deriv      = pif_total_link_deriv,
    conf_level      = conf_level,
    label           = sum_label,
    weights_sum_to_1 = FALSE
  )

  # ── 5. Build the ensemble ────────────────────────────────────────────────────
  pif_ens_obj <-
    pif_ensemble(
      pif1, ...,
      link            = pif_ensemble_link,
      link_inv        = pif_ensemble_link_inv,
      link_deriv      = pif_ensemble_link_deriv,
      weights         = weights,
      var_weights     = var_weights,
      var_pif_weights = var_pif_weights,
      conf_level      = conf_level,
      label           = ens_label,
      quiet           = quiet,
      is_paf          = is_paf
    )


  # ── 6. Extract point estimates ───────────────────────────────────────────────
  pif_vals <- sapply(pif_list, function(x) coef(x))   # length n
  sum_pif  <- coef(pif_sum)     # = Σ PIF_i
  ens_pif  <- coef(pif_ens_obj)

  # ── 7. Variances on log scale for sum and ensemble ───────────────────────────
  #
  #   variance() dispatches to cov_total_pif(x, x), propagating the full
  #   recursive uncertainty already embedded in each object.
  #
  #   Var[ln X] ≈ Var[X] / X^2  (first-order delta method)

  var_sum     <- variance(pif_sum)
  var_ens     <- variance(pif_ens_obj)
  var_log_sum <- var_sum / sum_pif^2
  var_log_ens <- var_ens / ens_pif^2

  # ── 8. Cov(sum, ensemble) ────────────────────────────────────────────────────
  cov_sum_ens     <- cov_total_pif(pif_sum, pif_ens_obj, var_p = var_p, var_beta = var_beta, warning = FALSE)
  cov_log_sum_ens <- cov_sum_ens / (sum_pif * ens_pif)

  # ── 9. Build each adjusted fraction ─────────────────────────────────────────

  adjusted_list <- vector("list", n)

  for (i in seq_len(n)) {

    pif_i   <- pif_vals[i]
    var_i   <- variance(pif_list[[i]])   # cov_total_pif(pif_i, pif_i)

    # ── Point estimate ──────────────────────────────────────────────────────
    pif_adj <- (pif_i / sum_pif) * ens_pif

    # ── Var[ln PIF_i] ───────────────────────────────────────────────────────
    var_log_i <- var_i / pif_i^2

    # ── Cov(PIF_i, sum)  via cov_total_pif ──────────────────────────────────
    #   pif_list[[i]] is a child of pif_sum, so cov_total_pif correctly
    #   returns the covariance of PIF_i with Σ_j PIF_j.
    cov_i_sum     <- cov_total_pif(pif_list[[i]], pif_sum, var_p = var_p, var_beta = var_beta, warning = FALSE)
    cov_log_i_sum <- cov_i_sum / (pif_i * sum_pif)
    #cov_ensemble_atomic(pif_ensemble = pif_sum, pif_atomic = pif_list[[i]], var_p = var_p, var_beta = var_beta)

    # ── Cov(PIF_i, ensemble) via cov_total_pif ──────────────────────────────
    #   pif_list[[i]] is also a child of pif_ens_obj, so cov_total_pif
    #   propagates uncertainty through the ensemble transform correctly.
    cov_i_ens     <- cov_total_pif(pif_list[[i]], pif_ens_obj, var_p = var_p, var_beta = var_beta, warning = FALSE)
    cov_log_i_ens <- cov_i_ens / (pif_i * ens_pif)

    # ── Var[ln PIF_adj] (delta method on log scale) ─────────────────────────
    #
    #   ln PIF_adj^i = ln PIF_i - ln S + ln Ensemble
    #
    #   Var[ln PIF_adj^i]
    #     = Var[ln PIF_i] + Var[ln S] + Var[ln Ens]
    #     + 2*(Cov(ln PIF_i, ln Ens) - Cov(ln PIF_i, ln S) - Cov(ln Ens, ln S))

    var_log_adj <- var_log_i + var_log_sum + var_log_ens +
      2 * (cov_log_i_ens - cov_log_i_sum - cov_log_sum_ens)

    # ── Back-transform: Var[X] ≈ X^2 * Var[ln X] ────────────────────────────
    var_adj <- pif_adj^2 * var_log_adj
    var_adj <- max(var_adj, 0, na.rm = TRUE)

    if (!quiet && var_log_adj < 0) {
      cli::cli_warn(
        c(
          "Negative log-scale variance for adjusted fraction {.val {labels[i]}}.",
          "i" = "Covariance terms dominate; variance set to 0.",
          "i" = "Consider providing {.arg var_p} or {.arg var_beta} for improved accuracy."
        )
      )
    }

    adjusted_list[[i]] <- pif_class(
      pif        = pif_adj,
      variance   = var_adj,
      conf_level = conf_level,
      label      = paste0(labels[i], "_adj"),
      type       = types[1L],
      link       = identity,
      link_inv   = identity,
      link_deriv = function(x) 1
    )
  }

  names(adjusted_list) <- labels
  adjusted_list
}


#' Weighted Adjusted PAF
#'
#' Convenience wrapper around [weighted_adjusted_fractions()] for
#' population attributable fractions (PAF).
#'
#' @inheritParams weighted_adjusted_fractions

#' @param paf1 A population attributable fraction
#' @param pif1 A potential impact fraction
#'
#' @return A named list of `pif_class` objects, one per input
#'   fraction, each being the weighted adjusted PAF (or PIF, respectively).
#'
#' @examples
#' paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001, label = "Lead")
#' paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001, label = "Radiation")
#'
#' weighted_adjusted_paf(paf_lead, paf_rad)
#'
#' pif_lead <- pif(0.2, p_cft = 0.1, beta = log(2.2), quiet = TRUE,
#'                 var_p = 0.001, label = "Lead")
#' pif_rad  <- pif(0.1, p_cft = 0.05, beta = log(1.2), quiet = TRUE,
#'                 var_p = 0.0001, label = "Radiation")
#'
#' weighted_adjusted_pif(pif_lead, pif_rad)
#'
#' @seealso [weighted_adjusted_fractions()], [paf_ensemble()] [pif()]
#' @name weighted_adjusted
#' @export
weighted_adjusted_paf <- function(
    paf1,
    ...,
    weights         = NULL,
    var_weights     = 0,
    var_pif_weights = NULL,
    var_p                    = NULL,
    var_beta                 = NULL,
    conf_level               = 0.95,
    label_ensemble           = NULL,
    label_sum                = NULL,
    quiet                    = FALSE
) {

  if (paf1@type != "PAF"){
    cli::cli_abort(
      "All fractions should be population attributable fractions (PAF). Otherwise use `weighted_adjusted_pif`."
    )
  }

  weighted_adjusted_fractions(
    paf1,
    ...,
    weights         = weights,
    var_weights     = var_weights,
    var_pif_weights = var_pif_weights,
    var_p                    = var_p,
    var_beta                 = var_beta,
    conf_level               = conf_level,
    label_ensemble           = label_ensemble,
    label_sum                = label_sum,
    quiet                    = quiet
  )
}


#' @seealso [weighted_adjusted_fractions()], [pif_ensemble()]
#' @rdname weighted_adjusted
#' @export
weighted_adjusted_pif <- function(
    pif1,
    ...,
    weights         = NULL,
    var_weights     = 0,
    var_pif_weights = NULL,
    var_p                    = NULL,
    var_beta                 = NULL,
    pif_total_link           = "log-complement",
    pif_total_link_inv       = NULL,
    pif_total_link_deriv     = NULL,
    pif_ensemble_link        = "identity",
    pif_ensemble_link_inv    = NULL,
    pif_ensemble_link_deriv  = NULL,
    conf_level               = 0.95,
    label_ensemble           = NULL,
    label_sum                = NULL,
    quiet                    = FALSE
) {
  weighted_adjusted_fractions(
    pif1,
    ...,
    weights         = weights,
    var_weights     = var_weights,
    var_pif_weights = var_pif_weights,
    var_p                    = var_p,
    var_beta                 = var_beta,
    pif_total_link           = "log-complement",
    pif_total_link_inv       = NULL,
    pif_total_link_deriv     = NULL,
    pif_ensemble_link        = "identity",
    pif_ensemble_link_inv    = NULL,
    pif_ensemble_link_deriv  = NULL,
    conf_level               = conf_level,
    label_ensemble           = label_ensemble,
    label_sum                = label_sum,
    quiet                    = quiet
  )
}
