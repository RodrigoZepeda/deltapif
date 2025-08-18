#' Covariance between two atomic pifs
#'
#' Calculates the approximate covariance between two atomic
#' pif class expressions using the delta method.
#'
#' @param pif1 A `pif_atomic_class` object.
#' @param pif2 A second `pif_atomic_class` object.
#' @param var_p covariance matrix for the prevalences in both `pif1` and `pif2`.
#' Default `NULL` (no correlation).
#' @param var_beta covariance matrix for the parameter `beta` in both `pif1`
#' and `pif2`. Default `NULL` (no correlation).
#' @param uncorrelated_p If `pif1` and `pif2` were constructed with independent
#' (uncorrelated) data for the prevalence `p`. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param uncorrelated_beta If `pif1` and `pif2` were constructed with independent
#'  (uncorrelated) data for the relative risk parameter `beta`. Either
#' `TRUE`, `FALSE` or `guess` (default).
#'
#' @param quiet Whether to throw warnings or other messages.
#'
#' @section Covariance matrices:
#'
#' By default if `var_p = NULL` is specified this assumes the parameters `p`
#' of `pif1` and `pif2` are uncorrelated. However, if `pif1` and `pif2` share
#' the same prevalence estimates
#' (i.e. share the same `p`s, then the user should set `uncorrelated_p = FALSE`
#' to account for that correlation).
#'
#' The same thing happens with `var_beta`. If no `var_beta` is specified (set to `NULL`)
#' then the assumption is that the `beta` parameters from `pif1` and
#' from `pif2` are uncorrelated (unless `uncorrelated_beta` is set to `FALSE`).
#'
#' If the user provides a covariance matrix for `var_p` then `uncorrelated_p`
#' is disregarded. Similarly, if the user provides `var_beta` then
#' `uncorrelated_beta` is ignored.
#'
#' @section Computation:
#'
#' The  `cov_atomic_pif` computes the approximate covariance between
#' two potential impact fractions \eqn{\text{PIF}_1} and \eqn{\text{PIF}_2}
#' as follows:
#' \deqn{
#' \text{Cov}\Big(\text{PIF}_1, \text{PIF}_2 \Big) \approx \Bigg(
#'  \dfrac{\partial \text{PIF}_1}{\partial p_1}^{\top}
#'  \Sigma_p
#'  \dfrac{\partial \text{PIF}_2}{\partial p_2} +
#'  \dfrac{\partial \text{PIF}_1}{\partial \beta_1}^{\top}
#'  \Sigma_{\beta}
#'  \dfrac{\partial \text{PIF}_2}{\partial \beta_2}
#' \Bigg)
#' }
#' where \eqn{\Sigma_p} is `var_p` and  \eqn{\Sigma_{\beta}} is `var_beta`.
#' The parameters \eqn{p_1} and \eqn{\beta_1} refer to the prevalence (`p`) and
#' relative risk exposure parameter `beta` of the first potential
#' impact fraction `pif1` and  \eqn{p_2} and \eqn{\beta_2} to the respective
#' parameters of `pif2`.
#'
#'
#' @seealso [from_parameters_covariance_p_component()]
#' @keywords internal
cov_atomic_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                           uncorrelated_p = "guess", uncorrelated_beta = "guess",
                           quiet = FALSE) {

  # Check they are pif objects
  if (!S7::S7_inherits(pif1, pif_atomic_class) ||
      !S7::S7_inherits(pif2, pif_atomic_class)) {
    cli::cli_abort(
      paste0(
        "One of the arguments passed to link_covariance estimation is not a ",
        "`deltapif::pif_atomic_class` object."
      )
    )
  }


  # Check that the ps don't appear similar
  if (uncorrelated_p == "guess" && is.null(var_p) &&
      length(pif1@p) == length(pif2@p) &&
      all(pif1@p == pif2@p) && all(pif1@var_p == pif2@var_p)) {
    uncorrelated_p <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The prevalence parameters `p` for the potential impact fractions appear ",
          "to be the same. If they are, set `uncorrelated_p = FALSE`. Otherwise ",
          "set `uncorrelated_p = TRUE"
        )
      )
    }
  }

  if (uncorrelated_beta == "guess" && is.null(var_beta) &&
      length(pif1@beta) == length(pif2@beta) &&
      all(pif1@beta == pif2@beta) &&
      all(pif1@var_beta == pif2@var_beta)) {
    uncorrelated_beta <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The relative risk `beta` parameters for the potential impact fractions appear ",
          "to be the same. If they are, set `uncorrelated_beta = FALSE`. Otherwise ",
          "set `uncorrelated_beta = TRUE"
        )
      )
    }
  }

  # If p's are the same set var_p as the link_covariance, otherwise assume independence
  if (is.null(var_p) && (uncorrelated_p == "guess" || uncorrelated_p)) {
    var_p <- matrix(0, nrow = length(pif1@p), ncol = length(pif2@p))
  } else if (is.null(var_p)) {
    var_p <- as.matrix(pif1@var_p)
  } else {
    var_p <- as.matrix(var_p)
  }

  # If beta's are the same set var_beta as the link_covariance, otherwise assume independence
  if (is.null(var_beta) && (uncorrelated_beta == "guess" || uncorrelated_beta)) {
    var_beta <- matrix(0, nrow = length(pif1@beta), ncol = length(pif2@beta))
  } else if (is.null(var_beta)) {
    var_beta <- as.matrix(pif1@var_beta)
  } else {
    var_beta <- as.matrix(var_beta)
  }


  if (nrow(var_p) != length(pif1@p) || ncol(var_p) != length(pif2@p)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `var_p`. It should be an ",
        "{length(pif1@p)} x {length(pif2@p)} matrix."
      )
    )
  }

  if (nrow(var_beta) != length(pif1@beta) || ncol(var_beta) != length(pif2@beta)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `var_beta`. It should be an ",
        "{length(pif1@beta)} x {length(pif2@beta)} matrix."
      )
    )
  }


  from_parameters_pif_covariance(
    p1       = pif1@p,
    p2       = pif2@p,
    p1_cft   = pif1@p_cft,
    p2_cft   = pif2@p_cft,
    rr1      = pif1@rr,
    rr2      = pif2@rr,
    mu_obs1  = pif1@mu_obs,
    mu_obs2  = pif2@mu_obs,
    mu_cft1  = pif1@mu_cft,
    mu_cft2  = pif2@mu_cft,
    var_p    = var_p,
    var_beta = var_beta,
    rr_link_deriv_vals1 = pif1@rr_link_deriv_vals,
    rr_link_deriv_vals2 = pif2@rr_link_deriv_vals,
  )

}


#' Covariance between the j-th weight of a `pif_global_ensemble_class` and
#' another `pif_global_ensemble_class`
#'
#' Calculates the covariance of the jth-weight
#' of a potential impact fraction of `pif_global_ensemble_class`
#' with a second `pif_global_ensemble_class` or `pif_atomic_class`.
#'
#' @param pif1 A `pif_global_ensemble_class` from which the weight is taken.
#'
#' @param pif2 A `pif_global_ensemble_class` or `pif_atomic_class` to compute the covariance
#'
#' @param j Weight indicator (covariance is for `weights(pif1)[j]`)
#'
#' @param sigma_weights Covariance vector between the weights of
#' `pif1` and the j-th weight of `pif2`.
#'
#' @param sigma_pif_weights Covariance vector between the potential
#' impact fractions in `pif1` and the j-th weight in `pif2`.
#'
#' @param uncorrelated_weights Whether the weights from `pif1` and `pif2`
#' are uncorrelated. By default it tries to `guess` by analyzing whether
#' they have the same weights.
#'
#' @param uncorrelated_pif_weights
#'
#' @section Formula:
#' Given a `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{1} =
#' g^{-1}\Bigg(\sum\limits_{i=1}^{M_1} g(\hat{q}_i \cdot \widehat{\textrm{PIF}}_{1,i})\Bigg)
#' }
#' and a second `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{2} = h^{-1}\Bigg(\sum\limits_{j=1}^{M_2} g(\hat{w}_j \cdot
#' \widehat{\textrm{PIF}}_{2,j})\Bigg)
#' }
#' This function computes the covariance between the first impact fraction and the
#' jth weight \eqn{w_j} of the second fraction by computing:
#' \deqn{
#' \operatorname{Cov}(\widehat{\textrm{PIF}}_A, \hat{w}_j) \approx
#' \frac{1}{g'\big(\widehat{\textrm{PIF}}_A\big)}\sum\limits_{i=1}^{M_1}g'(\hat{q}_i
#' \cdot \widehat{\textrm{PIF}}_{A,i}) \Bigg[ \hat{q} \operatorname{Cov}\Big(
#' \widehat{\textrm{PIF}}_{A,i}, \hat{w}_j\Big) +  \widehat{\textrm{PIF}}_{A,i} \operatorname{Cov}\Big( \hat{q}_j,
#' \hat{w}_j\Big) \Bigg]
#' }
#' where \eqn{\widehat{\textrm{PIF}}_{1,:} = (\widehat{\textrm{PIF}}_{1,1},
#' \widehat{\textrm{PIF}}_{1,2}, \dots, \widehat{\textrm{PIF}}_{1,M_1})^{\top}}
#'
#' @note The model currently works under the assumption that all `pif_atomic_class`
#' are independent from weights. Hence will return `0` if `pif2` is an atomic pif.
#'
#' @seealso [cov_atomic_pif()], [cov_ensemble_atomic()]
#' @keywords internal
cov_ensemble_weight <- function(pif1, pif2, j = 1, sigma_weights = NULL, sigma_pif_weights = NULL,
                                 uncorrelated_weights = "guess", uncorrelated_pif_weights = "guess",
                                 quiet = FALSE){

  #Otherwise
  if (!S7::S7_inherits(pif1, pif_global_ensemble_class)){
    cli::cli_abort(
      "Variable `pif1` must be a `pif_global_ensemble_class` object."
    )
  }

  if (!S7::S7_inherits(pif2, pif_global_ensemble_class) && !S7::S7_inherits(pif2, pif_atomic_class)){
    cli::cli_abort(
      "Variable `pif2` must be a `pif_global_ensemble_class` or a `pif_atomic_class` object."
    )
  }

  #By assumption weights are uncorrelated with atomic pifs
  if (S7::S7_inherits(pif2, pif_atomic_class)){
    return(0)
  }

  #If they have the same weights they are correlated and the correlation is in sigma_weights
  if (is.null(sigma_weights) &&
      uncorrelated_weights == "guess" &&
      length(pif1@weights) == length(pif2@weights) &&
      all(pif1@weights == pif2@weights) &&
      all(pif1@sigma_weights == pif2@sigma_weights)
      ) {

    uncorrelated_weights <- FALSE

    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The weights the potential impact fractions appear ",
          "to be the same. If they are, set `uncorrelated_weights = FALSE`. Otherwise ",
          "set `uncorrelated_weights = TRUE"
        )
      )
    }

    sigma_weights <- pif1@sigma_weights[,j]

  }


  if (is.null(sigma_weights) && (uncorrelated_weights == "guess" || uncorrelated_beta)) {
    sigma_weights <- rep(0, length(pif2@weights))
  }

  #In this case we compute the new weights
  if (is.null(sigma_pif_weights)){
    sigma_pif_weights <- rep(NA, length(pif2@coefs))
    for (k in 1:length(pif2@coefs)){
      sigma_pif_weights[k] <- cov_ensemble_weight(pif1 = pif1,
                                                  pif2 = pif2@pif_list[[k]],
                                                  j = j,
                                                  sigma_weights = sigma_weights,
                                                  sigma_pif_weights = NULL) #FIXME: Check this weight inheritance
    }
  }

  #Check the dimensions
  if (length(sigma_weights) != length(pif2@weights)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `sigma_weights` should be a vector of ",
        "length = {length(pif2@weights)}"
      )
    )
  }

  #Check the dimensions
  if (length(sigma_pif_weights) != length(pif2@coefs)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `sigma_pif_weights` should be a vector of ",
        "length = {length(pif2@coefs)}"
      )
    )
  }

  #Get the factors involved in covariance
  cov_val <- 0
  for (k in 1:length(pif2@coefs)){
    cov_val <- cov_val +
      pif2@pif_deriv_transform(pif2@weights[k]*pif2@coefs[k])*(
        pif2@weights[k]*sigma_pif_weights[k] +
          pif2@coefs[k]*sigma_weights[k]
    )
  }

  cov_val <- cov_val / pif2@pif_deriv_transform(pif2@pif)

  return(cov_val)

}

#' Covariance between `pif_global_ensemble_class` and a `pif_atomic`
#'
#' Calculates the covariance of a potential impact fraction of a
#' `pif_global_ensemble_class`  with a `pif_atomic`
#'
#' @param pif_ensemble A `pif_global_ensemble_class`
#'
#' @param pif_atomic A `pif_atomic_class`
#'
#' @param sigma_pifs Covariance vector between the potential
#' impact fractions in `pif_ensemble` and `pif_atomic`.
#'
#' @param sigma_weights_pif Covariance vector between the weights in  `pif_ensemble` and
#' the `pif_atomic`. For most applications this should be left as `NULL` which assumes
#' independence between the `pif_atomic` and the weights used in the `pif_ensemble`.
#'
#' @param uncorrelated_pifs Whether to "guess" if there is a correlation
#' between the atomic impact fraction and the impact fraction ensemble.
#'
#' @param uncorrelated_pif_weights Whether to "guess" if there is a correlation
#' between the atomic impact fraction and the weights. Currently ignored.
#'
#' @section Formula:
#' Given a `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{1} =
#' g^{-1}\Bigg(\sum\limits_{i=1}^{M_1} g(\hat{q}_i \cdot \widehat{\textrm{PIF}}_{1,i})\Bigg)
#' }
#' and a second `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{2} = h^{-1}\Bigg(\sum\limits_{j=1}^{M_2} g(\hat{w}_j \cdot
#' \widehat{\textrm{PIF}}_{2,j})\Bigg)
#' }
#' This function computes the covariance between the first impact fraction and the
#' coefficients of the second fraction by computing:
#' \deqn{
#' \operatorname{Cov}(\widehat{\textrm{PIF}}_A, \widehat{\textrm{PIF}}_{B,j})
#' \approx \frac{1}{g'\big(\widehat{\textrm{PIF}}_A\big)}\sum\limits_{i=1}^{M_1}g'(\hat{q}_i
#' \cdot \widehat{\textrm{PIF}}_{A,i}) \Bigg[ \mathbb{E}[\hat{q}_i] \operatorname{Cov}\Big( \textrm{PIF}_{A,i},
#' \widehat{\textrm{PIF}}_{B,j}\Big) +  \mathbb{E}[\textrm{PIF}_{A,i}] \operatorname{Cov}\Big( \hat{q}_i,
#' \widehat{\textrm{PIF}}_{B,j}\Big) \Bigg]
#' }
#' where \eqn{\widehat{\textrm{PIF}}_{1,:} = (\widehat{\textrm{PIF}}_{1,1},
#' \widehat{\textrm{PIF}}_{1,2}, \dots, \widehat{\textrm{PIF}}_{1,M_1})^{\top}}
#'
#' @seealso [cov_ensemble_weight()], [cov_atomic_pif()]
#' @keywords internal
cov_ensemble_atomic <- function(pif_ensemble, pif_atomic, sigma_pifs = NULL, sigma_weights_pif = NULL,
                                uncorrelated_pifs = "guess", uncorrelated_pif_weights = "guess"){

  if (!S7::S7_inherits(pif_ensemble, pif_global_ensemble_class) && !S7::S7_inherits(pif_ensemble, pif_atomic_class)){
    cli::cli_abort(
      "Entry `pif_ensemble` should be a `pif_global_ensemble_class` or a `pif_atomic` object."
    )
  }

  if (!S7::S7_inherits(pif_atomic, pif_atomic_class)){
    cli::cli_abort(
      "Entry `pif_atomic` should be a `pif_atomic_class` object."
    )
  }

  if (S7::S7_inherits(pif_ensemble, pif_atomic_class)){
    return(
      cov_atomic_pif(pif_ensemble, pif_atomic) #FIXME: Fix inheritance here
    )
  }

  #Guess the uncorrelated weights
  if (uncorrelated_pifs == "guess" && is.null(sigma_pifs)){
    sigma_pifs <- rep(0, length(pif_ensemble@coefs))
    for (k in 1:length(pif_ensemble@pif_list)){
      sigma_pifs[k] <- cov_ensemble_atomic(pif_ensemble = pif_ensemble@pif_list[[k]], pif_atomic = pif_atomic) #FIXME: Fix inheritance here
    }
  } else if (is.null(sigma_pifs)){
    sigma_pifs <- rep(0, ncol = length(pif_ensemble@coefs))
  }

  if (is.null(sigma_weights_pif)){
    sigma_weights_pif <- rep(0, length(pif_ensemble@coefs))
  }

  #Check the dimensions
  if (length(sigma_pifs) != length(pif_ensemble@coefs)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `sigma_pifs` should be a vector of ",
        "length {length(pif_ensemble@coefs)}"
      )
    )
  }

  #Check the dimensions
  if (length(sigma_weights_pif) != length(pif_ensemble@coefs)){
    cli::cli_abort(
      paste0(
        "Invalid dimensions for `sigma_weights_pif` should be a vector of ",
        "length {length(pif_ensemble@weights)}"
      )
    )
  }

  covariance_value <- 0
  for (k in 1:length(pif_ensemble@coefs)){
    covariance_value <- covariance_value +
      pif_ensemble@pif_deriv_transform(pif_ensemble@weights[k]*pif_ensemble@coefs[k])*(
        pif_ensemble@weights[k]*sigma_pifs[k] +
          pif_ensemble@coefs[k]*sigma_weights_pif[k]
      )
  }
  covariance_value <- covariance_value / pif_ensemble@pif_deriv_transform(pif_ensemble@pif)

  return(covariance_value)
}



#' Covariance function for a `pif_global_ensemble_class`
#'
#' Recursively obtains the covariance matrix of a `pif_global_ensemble_class`
#'
#' @param pif1 Either a `pif_atomic_class` or a `pif_global_ensemble_class`
#' @param pif2 Either a `pif_atomic_class` or a `pif_global_ensemble_class`
#'
#' @inheritParams cov_atomic_pif
#'
#' @inheritSection cov_atomic_pif Covariance matrices
#'
#' @section Computation:
#'
#' @seealso [from_parameters_covariance_p_component()], [cov_atomic_pif()],
#' [cov_ensemble_atomic()], [cov_ensemble_weight()]
#'
#' @keywords internal
cov_total_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                          uncorrelated_p = "guess", uncorrelated_beta = "guess",
                          quiet = FALSE) {

  # Base case: both are pif_atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_atomic_class)) {
    return(
      cov_atomic_pif(pif1 = pif1, pif2 = pif2, var_p = var_p,
                     var_beta = var_beta, uncorrelated_p = uncorrelated_p,
                     uncorrelated_beta = uncorrelated_beta, quiet = quiet) #FIXME
    )
  }

  # If pif1 is an ensemble and pif2 is atomic
  if (S7::S7_inherits(pif1, pif_global_ensemble_class) && S7::S7_inherits(pif2, pif_atomic_class)) {
    return(
      cov_ensemble_atomic(pif_ensemble = pif1, pif_atomic = pif2,
                          sigma_pifs = NULL, sigma_weights_pif = NULL) #FIXME
    )
  }

  # If pif2 is an ensemble and pif1 is atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_global_ensemble_class)) {
    return(
      cov_ensemble_atomic(pif_ensemble = pif2, pif_atomic = pif1,
                          sigma_pifs = NULL, sigma_weights_pif = NULL) #FIXME
    )
  }

  # If both are ensemble
  if (S7::S7_inherits(pif1, pif_global_ensemble_class) && S7::S7_inherits(pif2, pif_global_ensemble_class)) {
    total <- 0
    for (i in seq_along(pif1@pif_list)) {

      multiplier <- pif1@pif_deriv_transform(pif1@pif)

      total <- total +
        pif1@pif_deriv_transform(pif1@weights[i]*pif1@coefs[i])*(
          #q*cov(pif_a, pif_b)
          pif1@weights[i]*cov_total_pif(
            pif1 = pif1@pif_list[[i]],
            pif2 = pif2,
            var_p = var_p,
            var_beta = var_beta,
            uncorrelated_p = uncorrelated_p,
            uncorrelated_beta = uncorrelated_beta, #FIXME:
            quiet = quiet
          ) +
          #pifa*cov(q, pif_b)
          pif1@coefs[i]*cov_ensemble_weight(
            pif1 = pif1,
            pif2 = pif2,
            j = i #FIXME:
          )
        )

      total <- total / multiplier
    }
    return(total)
  }

  # If we get here, unsupported types
  cli::cli_abort(
    "Unsupported types for covariance calculation"
  )
}


#' Covariance matrix, correlation matrix, variance and standard deviation
#' for potential impact fractions
#'
#' Computes the covariance (`covariance`) or correlation (`correlation`) for multiple
#' potential impact fractions and the variance `variance` and standard deviation
#' `standard_deviation`for a potential impact fractions.
#'
#' @param x A potential impact fraction
#'
#' @param ... Multiple additional potential impact fraction objects
#' separated by commas.
#'
#' @param var_p covariance matrix for the prevalences in all `pif1` and
#' the ones included in `...`.
#'
#' @param var_beta covariance matrix for the parameter `beta` in all `pif1`
#' and the ones included in `...`.
#'
#' @param uncorrelated_p If all the pifs share the same prevalence data. Either
#' `TRUE`, `FALSE`, `guess` (default) or a matrix. If a matrix is given then
#' `uncorrelated_p[i,j] = 1` if the i-th and j-th pifs share the same prevalence data
#' `uncorrelated_p[i,j] = 0` if the i-th and j-th pifs don't share the same prevalence data.
#'
#' @param uncorrelated_beta If all the pifs share the same `beta` parameter. Either
#' `TRUE`, `FALSE` or `guess` (default) or a matrix. If a matrix is given then
#' `uncorrelated_beta[i,j] = 1` if the i-th and j-th pifs share the same relative risk parameters
#' `uncorrelated_beta[i,j] = 0` if the i-th and j-th pifs don't share the same relative risk parameters.
#'
#' @param quiet Whether to throw warnings and other messages
#'
#' @examples
#' # Get the approximate link_variance of a pif object
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' variance(my_pif)
#'
#' # This is the same as link_covariance with just 1 pIF
#' covariance(my_pif)
#'
#' # Calculate the link_covariance between 3 fractions with shared relative risk
#' beta <- 0.3
#' var_beta <- 0.1
#' pif1 <- pif(0.5, 0.2, beta, var_p = 0.5 * (1 - 0.5) / 100, var_beta = var_beta)
#' pif2 <- pif(0.3, 0.1, beta, var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta)
#' pif3 <- pif(0.7, 0.3, beta, var_p = 0.7 * (1 - 0.7) / 100, var_beta = var_beta)
#' covariance(pif1, pif2, pif3, uncorrelated_beta = FALSE)
#'
#' # The link_covariance between a pif and itself only has the link_variance as entries
#' covariance(pif1, pif1, uncorrelated_beta = FALSE, uncorrelated_p = FALSE)
#'
#' # Or if there is a link_covariance structure between different betas you can specify with
#' # var_beta in the link_covariance
#' betas <- c(1.3, 1.2, 1.27)
#'
#' # link_covariance among all betas
#' var_beta <- matrix(c(
#'   1.0000000, -0.12123053, 0.35429369,
#'   -0.1212305, 1.00000000, -0.04266409,
#'   0.3542937, -0.04266409, 1.00000000
#' ), byrow = TRUE, ncol = 3)
#' pif1 <- pif(0.5, 0.2, betas[1], var_p = 0.5 * (1 - 0.5) / 100, var_beta = var_beta[1, 1])
#' pif2 <- pif(0.3, 0.1, betas[2], var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta[2, 2])
#' pif3 <- pif(0.3, 0.1, betas[3], var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta[3, 3])
#' covariance(pif1, pif2, pif3, var_beta = var_beta)
#'
#' # Compute the correlation
#' correlation(pif1, pif2, pif3, var_beta = var_beta, quiet = TRUE)
#' @name covcor
NULL

#' @rdname covcor
#' @export
covariance <- S7::new_generic(
  "covariance", "x",
  function(x, ..., var_p = NULL, var_beta = NULL,
           uncorrelated_p = "guess", uncorrelated_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(covariance, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL,
                                                                              uncorrelated_p = "guess", uncorrelated_beta = "guess",
                                                                              quiet = FALSE) {

  # Get the list of fractions
  pif_list <<- append(list(x), list(...))
  npifs    <- length(pif_list)

  # If the independent ps are not a matrix nor "guess" then assign the same value
  if (!is.matrix(uncorrelated_p)) {
    uncorrelated_p <- matrix(uncorrelated_p, ncol = npifs, nrow = npifs)
  }

  if (!is.matrix(uncorrelated_beta)) {
    uncorrelated_beta <- matrix(uncorrelated_beta, ncol = npifs, nrow = npifs)
  }

  if (is.matrix(var_beta) && (ncol(var_beta) != npifs || nrow(var_beta) != npifs)){
    cli::cli_abort(
      paste0(
        "`var_beta` has dimensions {nrow(var_beta)} x {ncol(var_beta)} ",
        "but should be {npifs} x {npifs}"
      )
    )
  }

  if (is.matrix(var_p) && (ncol(var_p) != npifs || nrow(var_p) != npifs)){
    cli::cli_abort(
      paste0(
        "`var_p` has dimensions {nrow(var_p)} x {ncol(var_p)} ",
        "but should be {npifs} x {npifs}"
      )
    )
  }

  # link_covariance matrix
  if (npifs > 1) {
    cov_mat <- matrix(0, ncol = npifs, nrow = npifs)
    for (i in 1:(npifs - 1)) {
      for (j in (i + 1):npifs) {

        if (!is.null(var_p)) {
          nps <-  pif_class_apply_1st(x, length, "p")
          sub_var_p <- var_p[(nps * (i - 1) + 1):(nps * i), (nps * (j - 1) + 1):(nps * j)]
        } else {
          sub_var_p <- NULL
        }

        if (!is.null(var_beta)) {
          nbetas <- nps <-  pif_class_apply_1st(x, length, "beta")
          sub_var_beta <- var_beta[(nbetas * (i - 1) + 1):(nbetas * i), (nbetas * (j - 1) + 1):(nbetas * j)]
        } else {
          sub_var_beta <- NULL
        }


        cov_mat[i, j] <- cov_total_pif(pif_list[[i]], pif_list[[j]],
                                       var_p = sub_var_p, var_beta = sub_var_beta,
                                       uncorrelated_p = uncorrelated_p[i, j],
                                       uncorrelated_beta = uncorrelated_beta[i, j],
                                       quiet = quiet)
      }
    }

    # Add lower triangle
    cov_mat <- cov_mat + t(cov_mat) + diag(sapply(pif_list, variance))

  } else {
    cov_mat <- matrix(variance(x), ncol = 1, nrow = 1)
  }

  return(cov_mat)
}

#' @rdname covcor
#' @export
variance <- S7::new_generic("variance", "x")
S7::method(variance, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ...) {
  if (length(list(...)) > 0) {
    cli::cli_warn(
      "Currently this function does not support more than 1 argument. Ignoring the rest."
    )
  }
  #Get the covariance with uncorrelated = FALSE as they are correlated being the same pif
  cov_total_pif(x, x, uncorrelated_p = FALSE, uncorrelated_beta = FALSE)
}

#' @rdname covcor
#' @export
standard_deviation <- S7::new_generic("standard_deviation", "x")
S7::method(standard_deviation, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ...) {
  sqrt(variance(x, ...))
}

#' @rdname covcor
#' @export
correlation <- S7::new_generic(
  "correlation", "x",
  function(x, ..., var_p = NULL, var_beta = NULL,
           uncorrelated_p = "guess", uncorrelated_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(correlation, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL,
                                                                              uncorrelated_p = "guess", uncorrelated_beta = "guess",
                                                                              quiet = FALSE) {
  cov2cor(
    covariance(x, ...,
        var_p = var_p, var_beta = var_beta,
        uncorrelated_p = uncorrelated_p, uncorrelated_beta = uncorrelated_beta,
        quiet = quiet
    )
  )
}
