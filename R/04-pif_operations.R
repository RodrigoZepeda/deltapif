#' Combine Potential Impact Fractions and Population Attributable Fractions
#'
#' Combine potential impact fractions or the population attributable fractions
#' to either generate the total fraction from the fractions of
#' subpopulations (`pif_total/paf_total`) or the ensemble fraction of a population
#' from different (independent) exposures.
#'
#' @param pif1 A potential impact fraction (class `pif_class`)
#'
#' @param paf1 A population attributable fraction (class `pif_class`)
#'
#' @param ... The remaining potential impact fractions or (respectively)
#' population attributable fractions.
#'
#' @param weights A vector containing the proportion of the population
#' for each of the categories (for each of the pifs given).
#'
#' @param sigma_weights Colink_variance structure for the `weights`. Can be `0` (default) if
#' the weights are not random, a vector if only the link_variances of the weights
#' are available or a colink_variance matrix.
#'
#' @inheritParams pifpaf
#'
#' @section Total potential impact fraction:
#'
#' Assuming the overall population can be subdivided into \eqn{N} distinct
#' subpopulations each of them with a different potential impact fraction
#' (or population attributable fraction) we can estimate
#' the total population attributable fraction or potential impact
#' fraction of the whole population as:
#'
#' \deqn{
#'  \text{PIF}_{\text{Total}} = \sum\limits_{i = 1}^{N} \pi_i \cdot \text{PIF}_i
#' }
#'
#' where each \eqn{\text{PIF}_i} corresponds to the potential impact
#' fraction of the i-th subpopulation and \eqn{\pi_i} correspond to
#' the proportion of the total population occupied by \eqn{\text{PIF}_i}.
#' The weights are such that \eqn{\sum_{i=1}^{N} \pi_i = 1}.
#'
#' @section Ensemble potential impact fraction:
#'
#' If a population is exposed to \eqn{K} different independent risk factors
#' then the ensemble impact fraction of the combination of those factors
#' can be written as:
#'
#' \deqn{
#'  \text{PIF}_{\text{Ensemble}} = 1 - \prod\limits_{\ell = 1}^{K} \Big(1 - \cdot \text{PIF}_{\ell}\Big)
#' }
#'
#' where each \eqn{\text{PIF}_{\ell}} corresponds to the potential impact
#' fraction of the \eqn{\ell}-th risk factor for the same population.
#'
#' @examples
#' #Potential impact fraction for women
#' pif_women <- pif(0.32, 0.1, 1.2, quiet = TRUE, var_p = 0.1)
#'
#' #Potential impact fraction for men
#' pif_men <- pif(0.27, 0.1, 1.3, quiet = TRUE, var_p = 0.1)
#'
#' #Population potential impact fraction with 49% men and 51% women
#' pif_total(pif_men, pif_women, weights = c(0.49, 0.51), link = "logit")
#'
#' #Population attributable  fraction for women
#' paf_women <- paf(0.32, 1.3, quiet = TRUE, var_p = 0.1)
#'
#' #Population attributable  fraction for men
#' paf_men <- paf(0.27, 1.3, quiet = TRUE, var_p = 0.1)
#' paf_total(paf_men, paf_women, weights = c(0.49, 0.51), link = "logit")
#'
#' # Calculate the ensemble from lead and radiation exposure
#' paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001)
#' paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001)
#' pif_ensemble(paf_lead, paf_rad)
#'
#' # Totals and ensembles can be combined
#' pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001)
#' pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001)
#' pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women)
#' pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001)
#' pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001)
#' pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men)
#'
#' pif_total(pif_men, pif_women, weights = c(0.49, 0.51))
#' @name totalpifpaf
NULL

#' @rdname totalpifpaf
#' @export
paf_total <- function(paf1, ..., weights, sigma_weights = 0, conf_level = 0.95,
                      link = "log-complement",
                      link_inv = NULL,
                      link_deriv = NULL,
                      quiet = FALSE){

  paf_list <- append(list(paf1), list(...))
  for (i in seq_along(length(paf_list))){
    if (fraction_type(paf_list[[i]]) == "PIF"){
      cli::cli_abort(
        paste0(
          "Element {i} is not a Population Attributable Fraction. Did you mean to use `pif_total`?"
        )
      )
    }
  }

  pif_total(paf1, ..., weights = weights, sigma_weights = sigma_weights,
            conf_level = conf_level, link = link, link_inv = link_inv,
            link_deriv = link_deriv, quiet = quiet)

}

#' @rdname totalpifpaf
#' @export
pif_total <- function(pif1, ..., weights, sigma_weights = 0, conf_level = 0.95,
                      link = "log-complement",
                      link_inv = NULL,
                      link_deriv = NULL,
                      quiet = FALSE){

  #Get the fractions in list form
  pif_list <- append(list(pif1), list(...))
  npifs    <- length(pif_list)

  # Get the links
  link_name <- link #Save for later evaluation
  if (is.null(link_inv) & is.character(link)) {
    link_inv <- parse_inv_link(link)
  }

  # Get the functions
  link    <- parse_link(link)

  # Get the derivative
  if (is.null(link_deriv)) {
    link_deriv <- Deriv::Deriv(link)
  }

  if (abs(sum(weights) - 1) > sqrt(.Machine$double.eps)){
    cli::cli_abort(
      "`weights` should sum to 1 but `sum(weights)` = {sum(weights)}"
    )
  }

  #Check sigma weights
  if (length(sigma_weights) == 1){
    sigma_weights <- matrix(sigma_weights, nrow = npifs, ncol = npifs)
  } else if (is.vector(sigma_weights)) {
    sigma_weights <- sigma_weights %*% t(sigma_weights)
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "Assuming parameters `weights` are correlated but correlation is unknown. ",
          "If they are uncorrelated redefine `sigma_weights = diag(sigma_weights)` to ",
          "transform them into a matrix with no correlations."
        )
      )
    }
  } else if (is.matrix(sigma_weights)){

    if (ncol(sigma_weights) != npifs){
      cli::cli_abort(
        paste0(
          "Matrix `sigma_weights` has incorrect dimensions. Should be an ",
          "{npifs} x {npifs} matrix"
        )
      )
    }

    if (!isSymmetric(sigma_weights, trans = "T")){
      cli::cli_abort("Matrix `sigma_weights` is not symmetric.")
    }
  } else {
    cli::cli_abort(
      "`sigma_weights` should be a number, vector or matrix."
    )
  }

  pif <- pif_total_class(pif_list = pif_list,
                  weights = weights,
                  sigma_weights = sigma_weights,
                  conf_level = conf_level,
                  link = link,
                  link_inv = link_inv,
                  link_deriv = link_deriv)

  if (is.character(link_name) && link_name == "logit" && coef(pif) <= 0){
    cli::cli_alert_danger(
      paste0(
        "Value for {fraction_type(pif)} = {round(coef(pif),2)} <= 0. ",
        "Change link to a different value as `logit` is only ",
        "valid for strictly positive {fraction_type(pif)}s."
      )
    )
  }
  return(pif)
}

#' @rdname totalpifpaf
#' @export
pif_ensemble <- function(pif1, ..., conf_level = 0.95, quiet = FALSE){

  #Get the fractions in list form
  pif_list <- append(list(pif1), list(...))
  pif_ensemble_class(pif_list = pif_list, conf_level = conf_level)

}
