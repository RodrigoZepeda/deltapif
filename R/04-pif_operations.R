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
#' @param var_weights link_covariance structure for the `weights`. Can be `0` (default) if
#' the weights are not random, a vector if only the link_variances of the weights
#' are available or a link_covariance matrix.
#'
#' @inheritParams pifpaf
#' @inheritParams classes
#' @param is_paf Whether the computed quantity is a population attributable fraction or not
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
#'  \text{PIF}_{\text{Ensemble}} = 1 - \prod\limits_{\ell = 1}^{K} \Big(1 - \pi_{\ell} \cdot  \text{PIF}_{\ell}\Big)
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


#' Verify the variance between pif weights
#' @keywords internal
verify_var_pif_weights <- function(sigma_list, pif_list) {
  return(TRUE)
  n <- length(pif_list)
#FIXME: This function doesn't verify anything right now
  # Case 1: NULL is valid
  if (is.null(sigma_list)) {
    return(TRUE)
  }

  # Case 2: Check if it's a list
  if (!is.list(sigma_list)) {
    cli::cli_abort(
      "`sigma_list` must be a list or NULL"
    )
  }

  # Case 3: Check length
  if (length(sigma_list) != n) {
    cli::cli_abort(
      "`var_pif_weights` must have length {n} but has length {length(sigma_list)}"
    )
  }

  # Check each element
  for (i in 1:n) {
    # Check if each row element is a list
    if (!is.list(sigma_list[[i]])) {
      cli::cli_abort(
        "Element {i} var_pif_weights must be a list"
      )
    }

    # Check length of each row
    if (length(sigma_list[[i]]) != n) {
      cli::cli_abort(
        "Element {i} must have length {n} but has length {length(sigma_list[[i]])}"
      )
    }

    for (j in 1:n) {
      # Check if element is a matrix
      if (!is.matrix(sigma_list[[i]][[j]]) && !is.null(sigma_list[[i]][[j]])) {
        cli::cli_abort(
          "Element [{i}][{j}] must be a matrix or `NULL`"
        )
      }

      if (S7::S7_inherits(pif_list[[i]], pif_atomic_class) || S7::S7_inherits(pif_list[[j]], pif_atomic_class)){
        if (is.matrix(sigma_list[[i]][[j]])){
          cli::cli_abort(
            "Element [{i}][{j}] must be `NULL` as either the i-th or j-th fraction has no weights (is atomic)"
          )
        }
      }

      # Check matrix dimensions
      if (S7::S7_inherits(pif_list[[i]], pif_ensemble_class) && S7::S7_inherits(pif_list[[j]], pif_ensemble_class)){
        if (is.matrix(sigma_list[[i]][[j]])){
          nrows <- length(weights(pif_list[[i]]))
          ncols <- length(weights(pif_list[[j]]))
          if (nrow(sigma_list[[i]][[j]]) != nrows || ncol(sigma_list[[i]][[j]]) != ncols) {
            cli::cli_abort(
              paste0(
                "Matrix at [{i}][{j}] must be {nrows} x {ncols}",
                "but is, {nrow(sigma_list[[i]][[j]])} x {ncol(sigma_list[[i]][[j]])}"
              )
            )
          }
        }
      }
    }
  }
  return(TRUE)
}

#' Prepare a global_ensemble
#'
#' Helper function that helps validate the inputs on any global ensemble
#'
#' @inheritParams totalpifpaf
#' @param weights_sum_to_1 Whether to check if the weights sum to 1
#'
#' @return A list of validated values for use in `pif_total` or `pif_ensemble`.
#' @keywords internal
pif_validate_ensemble <- function(pif1, ..., weights, var_weights,
                                  var_pif_weights, conf_level,
                                  link, link_inv, link_deriv, quiet,
                                  is_paf = FALSE, weights_sum_to_1 = FALSE,
                                  label){

  #Get the fractions in list form
  pif_list <- append(list(pif1), list(...))
  npifs    <- length(pif_list)

  #Check they are all pifs
  class_val <- sapply(pif_list, function(x) S7::S7_inherits(x, pif_atomic_class) || S7::S7_inherits(x, pif_global_ensemble_class))
  if (!all(class_val)){
    cli::cli_abort(
      "Element {which(!class_val)} is not of `pif_atomic_class` or `pif_global_ensemble_class`"
    )
  }

  #Check that labels are different
  label_names <- sapply(pif_list, function(x) x@label)
  dups        <- duplicated(label_names)
  if (any(dups)){
    cli::cli_abort(
      "Element {label_names[dups]} is duplicated. Use the `label` argument to set up a different label. Cannot combine fractions with duplicated labels."
    )
  }

  # Get the inverse function
  if (is.null(link_inv) & is.character(link)) {
    link_inv <- parse_inv_link(link)
  }

  if (is.null(link_deriv) & is.character(link)){
    link_deriv <- parse_deriv_link(link)
  }

  # Get the functions
  link <- parse_link(link)

  if (!is.function(link_deriv) && is.null(link_deriv)) {
    link_deriv <- Deriv::Deriv(link)
  }

  #Check weights and set default to 1
  if (is.null(weights)){
    weights <- rep(1, length(pif_list))
  }


  #Check the weights
  if (weights_sum_to_1 && abs(sum(weights) - 1) > sqrt(.Machine$double.eps)){
    cli::cli_abort(
      "`weights` should sum to 1 but `sum(weights)` = {sum(weights)}"
    )
  }

  if (!is.numeric(var_weights)){
    cli::cli_abort(
      "`var_weights` should be a number, vector or matrix."
    )
  }

  if (is.matrix(var_weights)){
    if (length(weights) != ncol(var_weights) || length(weights) != nrow(var_weights)){
      cli::cli_abort(
        "Invalid dimensions for `var_weights`"
      )
    }
  }

  if (is.vector(var_weights) && length(var_weights) > 1){
    if (length(weights) != length(var_weights)){
      cli::cli_abort(
        "Invalid dimensions for `var_weights`"
      )
    }
  }

  #Validate the matrix
  #verify_var_pif_weights(var_pif_weights, pif_list)

  #Check in case is paf that each element is a paf
  if (is_paf){
    for (i in seq_along(1:length(pif_list))){
      if (fraction_type(pif_list[[i]]) == "PIF"){
        cli::cli_abort(
          paste0(
            "Element {i} is not a Population Attributable Fraction. Did you mean to use `pif_total`?"
          )
        )
      }
    }
  }

  #Check sigma weights
  if (length(var_weights) == 1 && is.numeric(var_weights)){
    var_weights <- matrix(var_weights, nrow = npifs, ncol = npifs)
  } else if (is.vector(var_weights)) {
    var_weights <- var_weights %*% t(var_weights)
    if (!quiet){
      cli::cli_warn(
        paste0(
          "Assuming parameters `weights` are correlated but correlation is unknown. ",
          "If they are uncorrelated redefine `var_weights = diag(var_weights)` to ",
          "transform them into a matrix with no correlations."
        )
      )
    }
  }

  #Check sigma weights
  if (length(var_pif_weights) == 1 && is.numeric(var_pif_weights)){
    var_pif_weights <- matrix(var_pif_weights, nrow = npifs, ncol = npifs)
  }

  if (is.matrix(var_pif_weights)){

    if (ncol(var_pif_weights) != npifs || nrow(var_pif_weights) != npifs){
      cli::cli_abort(
        paste0(
          "Matrix `var_pif_weights` has incorrect dimensions. Should be an ",
          "{npifs} x {npifs} matrix"
        )
      )
    }

    if (!isSymmetric(var_pif_weights, trans = "T")){
      cli::cli_abort("Matrix `var_pif_weights` is not symmetric.")
    }
  } else if (is.vector(var_pif_weights)) {
    cli::cli_abort(
      "Matrix `var_pif_weights` should be a matrix or `0`."
    )
  }

  if (is.null(label)){
    label <- paste0("deltapif-", sub("\\.", "", as.character(abs(stats::rnorm(1)))))
  }

  #Get the names into pif_list
  labels <- sapply(pif_list, function(x) x@label)
  names(pif_list) <- labels


  return(
    list(
      pif_list = pif_list,
      weights = weights,
      var_weights = var_weights,
      var_pif_weights = var_pif_weights,
      conf_level = conf_level,
      link = link,
      link_inv = link_inv,
      link_deriv = link_inv,
      label = label
    )
  )

}

#' @rdname totalpifpaf
#' @export
paf_total <- function(paf1, ..., weights, var_weights = 0,
                      var_pif_weights = NULL,
                      conf_level = 0.95,
                      link = "log-complement",
                      link_inv = NULL,
                      link_deriv = NULL,
                      quiet = FALSE,
                      label = NULL){


  pif_total(paf1, ..., weights = weights, var_weights = var_weights,
            var_pif_weights = var_pif_weights,
            conf_level = conf_level, link = link, link_inv = link_inv,
            link_deriv = link_deriv, quiet = quiet, is_paf = TRUE, label = label)

}

#' @rdname totalpifpaf
#' @export
pif_total <- function(pif1, ..., weights,
                      var_weights = 0,
                      var_pif_weights = NULL,
                      conf_level = 0.95,
                      link = "log-complement",
                      link_inv = NULL,
                      link_deriv = NULL,
                      quiet = FALSE,
                      label = NULL,
                      is_paf = FALSE){



  pif_params <- pif_validate_ensemble(pif1 = pif1, ..., weights = weights,
                                      var_weights = var_weights,
                                      var_pif_weights = var_pif_weights,
                                      conf_level = conf_level,
                                      link = link, link_inv = link_inv,
                                      link_deriv = link_deriv,
                                      quiet = quiet,
                                      label = label,
                                      is_paf = is_paf, weights_sum_to_1 = TRUE)


  pif <- do.call(pif_total_class, pif_params)

  if (is.character(link) && link == "logit" && coef(pif) <= 0){
    cli::cli_warn(
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
paf_ensemble <- function(paf1, ..., weights = NULL, var_weights = 0,
                         var_pif_weights = NULL,
                      link = "identity",
                      link_inv = NULL,
                      link_deriv = NULL,
                      conf_level = 0.95,
                      label = NULL,
                      quiet = FALSE){


  pif_ensemble(paf1, ..., weights = weights, var_weights = var_weights,
            var_pif_weights = var_pif_weights,
            conf_level = conf_level, link = link, link_inv = link_inv,
            link_deriv = link_deriv, quiet = quiet, is_paf = TRUE,
            label = label)

}

#' @rdname totalpifpaf
#' @export
pif_ensemble <- function(pif1, ...,
                         weights = NULL,
                         var_weights = 0,
                         var_pif_weights = NULL,
                         link = "identity",
                         link_inv = NULL,
                         link_deriv = NULL,
                         conf_level = 0.95,
                         quiet = FALSE,
                         label = NULL,
                         is_paf = FALSE){


  pif_params <- pif_validate_ensemble(pif1 = pif1, ..., weights = weights,
                                      var_weights = var_weights,
                                      var_pif_weights = var_pif_weights,
                                      conf_level = conf_level,
                                      link = link, link_inv = link_inv,
                                      link_deriv = link_deriv,
                                      quiet = quiet,
                                      label = label,
                                      is_paf = is_paf, weights_sum_to_1 = FALSE)

  pif <- do.call(pif_ensemble_class, pif_params)

  if (is.character(link) && link == "logit" && coef(pif) <= 0){
    cli::cli_warn(
      paste0(
        "Value for {fraction_type(pif)} = {round(coef(pif),2)} <= 0. ",
        "Change link to a different value as `logit` is only ",
        "valid for strictly positive {fraction_type(pif)}s."
      )
    )
  }
  return(pif)

}
