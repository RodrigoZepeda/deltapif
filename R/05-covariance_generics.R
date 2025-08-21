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
cov_atomic_pif <- function(pif1, pif2, var_p, var_beta) {

  # Check they are pif objects
  if (!S7::S7_inherits(pif1, pif_atomic_class) ||
      !S7::S7_inherits(pif2, pif_atomic_class)) {
    cli::cli_abort(
      paste0(
        "One of the arguments `pif1` or `pif2` passed to cov_atomic_pif estimation is not a ",
        "`deltapif::pif_atomic_class` object."
      )
    )
  }

  #Change to matrices just in case numbers were given
  if (S7::S7_inherits(var_p, covariance_structure_class)){
    var_p    <- as.matrix(var_p[[pif1@label]][[pif2@label]])
  } else {
    var_p    <- as.matrix(var_p)
  }

  if (S7::S7_inherits(var_beta, covariance_structure_class)){
    var_beta <- as.matrix(var_beta[[pif1@label]][[pif2@label]])
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
#' @param var_weights Covariance vector between the weights of
#' `pif1` and the j-th weight of `pif2`.
#'
#' @param sigma_pif_weights Covariance vector between the potential
#' impact fractions in `pif1` and the weights in `pif2`.
#'
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
cov_ensemble_weights <- function(pif1, pif2, var_weights, var_pif_weights, recursive = !is.null(var_pif_weights)){

  #Return 0 if pif1 or pif2 are atomic as there is no covariance between weights
  if (S7::S7_inherits(pif1, pif_atomic_class) || S7::S7_inherits(pif2, pif_atomic_class)){
    return(0)
  }

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


  if (recursive){
    for (k in 1:length(pif1@coefs)){
      name_pif_k <- names(pif1@pif_list)[k]
      var_pif_weights[[name_pif_k]][[pif2@label]] <-
        cov_ensemble_weights(pif1 = pif1@pif_list[[k]],
                             pif2 = pif2,
                             var_weights = var_weights,
                             var_pif_weights = var_pif_weights,
                             recursive = TRUE)
    }
  }

  # #Check the dimensions
  # if (!identical(pif1, pif2) && length(var_weights) != length(flatten_names(pif1)) +  length(flatten_names(pif2))){
  #   cli::cli_abort(
  #     paste0(
  #       "Invalid dimensions for `var_weights` should be a `covariance_structure` of ",
  #       "length {length(flatten_names(pif1)) +  length(flatten_names(pif2))}"
  #     )
  #   )
  # } else if (identical(pif1, pif2) && length(var_weights) != length(flatten_names(pif1))){
  #   cli::cli_abort(
  #     paste0(
  #       "Invalid dimensions for `var_weights` should be a `covariance_structure` of ",
  #       "length {length(flatten_names(pif1))}"
  #     )
  #   )
  # }


  #Get the factors involved in covariance
  alabel <- pif1@label
  klabel <- pif2@label


  #Loop applying the transform
  for (k in 1:length(pif1@coefs)){
    pif_deriv_transform_vals <- pif1@pif_deriv_transform(pif1@weights[k]*pif1@coefs[k])
  }

  #Expand the variances to matrix form if needed
  var_pif_w <- var_pif_weights[[klabel]][[alabel]]
  if (!is.matrix(var_pif_w)){
    var_pif_w <- matrix(var_pif_w, ncol = length(pif2@weights), nrow = length(pif1@weights))
  }

  var_w <- var_weights[[klabel]][[alabel]]
  if (!is.matrix(var_w)){
    var_w <- matrix(var_w, ncol = length(pif2@weights), nrow = length(pif1@weights))
  }

  #Transform so that it makes sense
  as.numeric(
    pif_deriv_transform_vals * (
      t(pif1@weights) %*% var_pif_w + t(pif1@coefs) %*% var_w
    ) / pif2@pif_deriv_transform(pif2@pif)
  )


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
#' @param var_pifs Covariance vector between the potential
#' impact fractions in `pif_ensemble` and `pif_atomic`.
#'
#' @param var_pif_weights Covariance vector between the weights in  `pif_ensemble` and
#' the `pif_atomic`.
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
cov_ensemble_atomic <- function(pif_ensemble, pif_atomic, var_p, var_beta,
                                var_pifs, var_pif_weights,
                                recursive = TRUE){


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
      cov_atomic_pif(pif1     = pif_ensemble,
                     pif2     = pif_atomic,
                     var_p    = var_p[[pif_ensemble@label]][[pif_atomic@label]],
                     var_beta = var_beta[[pif_ensemble@label]][[pif_atomic@label]])
    )
  }

  #Calculate the inner weights of the pifs
  if (recursive){
    #Loop through each fraction in the ensemble and repeat
    for (k in 1:length(pif_ensemble@pif_list)){
      name_pif_k <- names(pif_ensemble@pif_list)[k]
      var_pifs[[name_pif_k]][[pif_atomic@label]] <-
        cov_ensemble_atomic(
          pif_ensemble = pif_ensemble@pif_list[[k]],
          pif_atomic   = pif_atomic,
          var_p        = var_p,
          var_beta     = var_beta,
          var_pifs   = var_pifs,
          var_pif_weights = var_pif_weights,
          recursive    = TRUE)

      #Symmetry assignment
      var_pifs[[pif_atomic@label]][[name_pif_k]]  <- var_pifs[[name_pif_k]][[pif_atomic@label]]

    }
  }

  #Loop through the sum of 1/g'*g'(qpif)*(q cov(pif1, pifj) +pifi*cov(q*pifj)
  alabel <- pif_atomic@label
  gprime <- rep(NA, length(pif_ensemble@coefs))
  for (k in 1:length(gprime)){
    gprime[k] <- pif_ensemble@pif_deriv_transform(pif_ensemble@weights[k]*pif_ensemble@coefs[k])
  }
  var_pif         <- subset(var_pifs, cols = alabel, rows = children(pif_ensemble))
  var_pif_weight  <- subset(var_pif_weights, cols = alabel, rows = children(pif_ensemble))


  1/pif_ensemble@pif_deriv_transform(pif_ensemble@pif)*(
    t(gprime * pif_ensemble@weights) %*% as.matrix(var_pif) +
      t(gprime * pif_ensemble@coefs) %*% as.matrix(var_pif_weight)
  )

}



#' Covariance function for a `pif_global_ensemble_class`
#'
#' Recursively obtains the covariance matrix of a `pif_global_ensemble_class`
#'
#' @param pif1 Either a `pif_atomic_class` or a `pif_global_ensemble_class`
#' @param pif2 Either a `pif_atomic_class` or a `pif_global_ensemble_class`
#' @param var_weights Covariance between the weights of all the fractions
#' contained in `pif1` and all the fractions contained in `pif2`.
#'
#' @inheritParams cov_atomic_pif
#'
#' @inheritSection cov_atomic_pif Covariance matrices
#'
#' @section Computation:
#' This computes:
#' \deqn{
#' \operatorname{Cov}(\widehat{\textrm{PIF}}_{A}, \widehat{\textrm{PIF}}_{B})  \approx
#' \frac{1}{g'\big(\widehat{\textrm{PIF}}_A\big)}\sum\limits_{i=1}^{M_1}g'(\hat{q}_i \cdot
#' \widehat{\textrm{PIF}}_{A,i}) \Bigg[ \mathbb{E}[\hat{q}_i] \operatorname{Cov}\Big( \textrm{PIF}_{A,i},
#' \widehat{\textrm{PIF}}_{B}\Big) +  \mathbb{E}[\textrm{PIF}_{A,i}] \operatorname{Cov}\Big( \hat{q}_i,
#' \widehat{\textrm{PIF}}_{B}\Big) \Bigg]
#' }
#' where \eqn{\operatorname{Cov}\Big( \textrm{PIF}_{A,i},
#' \widehat{\textrm{PIF}}_{B}\Big)} is computed as in [cov_ensemble_atomic()]
#' in the case \eqn{\widehat{\textrm{PIF}}_{A,i}} is atomic or with the same formula
#' in the case it is an ensemble. The expression \eqn{\operatorname{Cov}\Big( \hat{q}_i,
#' \widehat{\textrm{PIF}}_{B}\Big) } is computed as in  [cov_ensemble_weight()]
#' for each weight \eqn{q_j}.
#'
#' @seealso [from_parameters_covariance_p_component()], [cov_atomic_pif()],
#' [cov_ensemble_atomic()], [cov_ensemble_weight()]
#'
#' @keywords internal
cov_total_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                          var_weights = NULL, var_pif_weights = NULL,
                          var_pifs = NULL,
                          quiet = FALSE) {

  #If varp and varbeta are null set defaults
  if (is.null(var_p)){
    var_p <- default_parameter_covariance_structure2(pif1, pif2, parameter = "p")
  } else if (!S7::S7_inherits(var_p, covariance_structure_class)){
    cli::cli_abort(
      paste0(
        "`var_p` should be a `covariance_structure_class`. Use ",
        "`default_parameter_covariance_structure2(pif1, pif2, parameter = 'p')`",
        "to create a prototype and work from there."
      )
    )
  }

  if (is.null(var_beta)){
    var_beta <- default_parameter_covariance_structure2(pif1, pif2, parameter = "beta")
  } else if (!S7::S7_inherits(var_beta, covariance_structure_class)){
    cli::cli_abort(
      paste0(
        "`var_beta` should be a `covariance_structure_class`. Use ",
        "`default_parameter_covariance_structure2(pif1, pif2, parameter = 'beta')`",
        "to create a prototype and work from there."
      )
    )
  }

  if (is.null(var_weights)){
    var_weights <- default_weight_covariance_structure2(pif1, pif2)
  } else if (!S7::S7_inherits(var_weights, covariance_structure_class)){
    cli::cli_abort(
      paste0(
        "`var_weights` should be a `covariance_structure_class`. Use ",
        "`default_weight_covariance_structure2(pif1, pif2, parameter = 'beta')`",
        "to create a prototype and work from there."
      )
    )
  }

  #If sigma_*_weights are null set defaults
  if (is.null(var_pif_weights)){
    var_pif_weights <- default_weight_pif_covariance_structure2(pif1, pif2)
  } else if (!S7::S7_inherits(var_pif_weights, covariance_structure_class)){
    cli::cli_abort(
      paste0(
        "`var_pif_weights` should be a `covariance_structure_class`. Use ",
        "`covariance_structure2(pif1, pif2)`",
        "to create a prototype and work from there."
      )
    )
  }

  if (is.null(var_pifs)){
    var_pifs <- default_pif_covariance_structure2(pif1, pif2)
  } else if (!S7::S7_inherits(var_pifs, covariance_structure_class)){
    cli::cli_abort(
      paste0(
        "`var_pifs` should be a `covariance_structure_class`. Use ",
        "`covariance_structure2(pif1, pif2)`",
        "to create a prototype and work from there."
      )
    )
  }

  # Base case: both are pif_atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_atomic_class)) {

    return(
      cov_atomic_pif(pif1     = pif1,
                     pif2     = pif2,
                     var_p    = var_p,
                     var_beta = var_beta)
    )
  }

  # If pif1 is an ensemble and pif2 is atomic
  if (S7::S7_inherits(pif1, pif_global_ensemble_class) && S7::S7_inherits(pif2, pif_atomic_class)) {
    return(
      cov_ensemble_atomic(pif_ensemble = pif1, pif_atomic = pif2,
                          var_p = var_p, var_beta = var_beta,
                          var_pifs = var_pifs,
                          var_pif_weights = var_pif_weights,
                          recursive = TRUE)
    )
  }

  # If pif2 is an ensemble and pif1 is atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_global_ensemble_class)) {
    return(
      cov_ensemble_atomic(pif_ensemble = pif2, pif_atomic = pif1,
                          var_p = var_p, var_beta = var_beta,
                          var_pifs = var_pifs,
                          var_pif_weights = var_pif_weights,
                          recursive = TRUE)
    )
  }

  # If both are ensemble
  if (S7::S7_inherits(pif1, pif_global_ensemble_class) && S7::S7_inherits(pif2, pif_global_ensemble_class)) {

    #Get the denominator
    multiplier <- pif1@pif_deriv_transform(pif1@pif)

    #Get covariance betweem the weights
    cov_weights <- cov_ensemble_weights(pif1 = pif1, pif2 = pif2, var_weights = var_weights,
                                        var_pif_weights = var_pif_weights, recursive = TRUE)

    #Get the covariance vector for the sub-pifs
    total_cov <- rep(NA, length(pif1@coefs))


    for (i in 1:length(total_cov)) {
      total_cov[i] <- cov_total_pif(pif1 = pif1@pif_list[[i]],
                                    pif2 = pif2, var_p = var_p,
                                    var_beta = var_beta,
                                    var_weights = var_weights,
                                    var_pifs = var_pifs,
                                    var_pif_weights = var_pif_weights)
    }


    return(
      as.numeric(
        pif1@pif_deriv_transform(pif1@weights[i]*pif1@coefs[i])*(
          pif1@weights %*% total_cov +
          pif1@coefs %*% cov_weights
        ) / multiplier
      )
    )

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
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(covariance, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL,
                                                                              quiet = FALSE) {

  # Get the list of fractions
  pif_list <- append(list(x), list(...))
  npifs    <- length(pif_list)



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
          sub_var_p <- var_p[(nps * (i - 1) + 1):(nps * i), (nps * (j - 1) + 1):(nps * j)] #FIXME: This should be a list
        } else {
          sub_var_p <- NULL
        }

        if (!is.null(var_beta)) {
          nbetas <- nps <-  pif_class_apply_1st(x, length, "beta")
          sub_var_beta <- var_beta[(nbetas * (i - 1) + 1):(nbetas * i), (nbetas * (j - 1) + 1):(nbetas * j)] #FIXME: This should be a list
        } else {
          sub_var_beta <- NULL
        }


        cov_mat[i, j] <- cov_total_pif(pif_list[[i]], pif_list[[j]],
                                       var_p = sub_var_p, var_beta = sub_var_beta)
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
  cov_total_pif(x, x)
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
