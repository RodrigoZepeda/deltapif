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
#'
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
#' @examples
#' \dontrun{
#' p1 <- pif(p = 0.5, p_cft = 0.1, beta = 0.2, 0.003, 0.04, label = "Population 1")
#' p2 <- pif(p = 0.4, p_cft = 0.1, beta = 0.2, 0.003, 0.04, label = "Population 2")
#'
#' #Positive covariance as program automatically detects same beta and p
#' cov_atomic_pif(p1, p2)
#'
#' #Zero covariance as they have no covariates in common
#' p3 <- pif(p = 0.6, p_cft = 0.1, beta = 0.12, 0.003, 0.04, label = "Population 3")
#' cov_atomic_pif(p1, p3)
#'
#' #Covariance has to be given to specify the covariance between p and/or beta parameters
#' cov_atomic_pif(p1, p3, var_p = 0.12, var_beta = 0.11)
#'
#' #Works with parameters of different dimensions with entry i-j being the
#' #covariance between entry p[i] of pif1 and entry p[j] of pif2
#' #(same with the betas)
#' p_dim_3 <- paf(p = c(0.5, 0.2, 0.1), beta = c(0.2, 0.1, 0.4),
#'   label = "Population 1", quiet = TRUE)
#' p_dim_2 <- pif(p = c(0.4, 0.1), p_cft = c(0.1, 0.05), beta = c(0.3, 0.8),
#'   label = "Population 2", quiet = TRUE)
#' cov_atomic_pif(p_dim_2, p_dim_3,
#'   var_p = matrix(c(0.1, 0.2, 0.4, 0.5, 0.6, 0.08), ncol = 3))
#' }
#'
#' @seealso [from_parameters_covariance_p_component()] to calculate the covariance
#' from first principles, [cov_ensemble_atomic] and [cov_ensemble_weights()]
#' for the covariance between an ensemble fraction and an atomic fraction
#' and the covariance between the weights of an ensemble fraction
#' and another fraction.
#'
#' @keywords internal
cov_atomic_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL) {
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

  #If varp and varbeta are null set defaults
  if (is.numeric(var_p)){
    var_p <- as.matrix(var_p)
  }

  #Check if is zero matrix
  if (identical(var_p, matrix(0, ncol = 1, nrow = 1))){
    var_p <- matrix(0, ncol = length(pif2@p), nrow = length(pif1@p))
  }

  if (!is.matrix(var_p)){
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
  } else {
    if (ncol(var_p) != length(pif2@p) || nrow(var_p) != length(pif1@p)){
      cli::cli_abort(
        paste0(
          "`var_p` should be a matrix of dimensions {length(pif1@p)} x {length(pif2@p)}"
        )
      )
    }
  }

  if (is.numeric(var_beta)){
    var_beta <- as.matrix(var_beta)
  }

  if (identical(var_beta, matrix(0, ncol = 1, nrow = 1))){
    var_beta <- matrix(0, ncol = length(pif2@beta), nrow = length(pif1@beta))
  }

  if (!is.matrix(var_beta)){
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
  } else {
    if (ncol(var_beta) != length(pif2@beta) || nrow(var_beta) != length(pif1@beta)){
      cli::cli_abort(
        paste0(
          "`var_beta` should be a matrix of dimensions {length(pif1@beta)} x {length(pif2@beta)}"
        )
      )
    }
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

  #Double check if they are collapsed to 0 make them bigger
  if (length(var_p) == 1 && (var_p == 0 || identical(var_p, matrix(0)))){
    var_p <- matrix(0, nrow = length(pif1@p), ncol = length(pif2@p))
  }

  if (length(var_beta) == 1 && (var_beta == 0 || identical(var_beta, matrix(0)))){
    var_beta <- matrix(0, nrow = length(pif1@beta), ncol = length(pif2@beta))
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


#' Covariance between a pif and a weight
#'
#' Calculates the covariance of the weight
#' of a potential impact fraction of `pif_global_ensemble_class`
#' with a second `pif_global_ensemble_class` or `pif_atomic_class`.
#'
#' @param pif1 A `pif_global_ensemble_class` from which the weight is taken.
#'
#' @param pif2 A `pif_global_ensemble_class` or `pif_atomic_class` to compute the covariance
#'
#' @param var_weights Covariance matrix between the weights of
#' `pif1` and the weights of `pif2`. Entry `var_weights[i,j]` is the
#' covariance between the `i`th weight of `pif1` and the `j`th weight of `pif2`.
#' This refers to the term \eqn{\operatorname{Cov}\Big( \hat{q}_i,\hat{w}_j\Big)}
#' in the equation before.
#'
#' @param var_pif_weights Covariance structure between the potential
#' impact fractions in `pif1` and the weights in `pif2`. Entry `var_pif_weights[i,j]` is the
#' covariance between the `i`th fraction of `pif1` and the `j`th weight of `pif2`.
#' This refers to the term \eqn{\operatorname{Cov}\Big(\widehat{\textrm{PIF}}_{A,i}, \hat{w}_j\Big)}
#' in the equation before.
#'
#' @section Formula:
#' Given a `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{1} =
#' g^{-1}\Bigg(\sum\limits_{i=1}^{M_1} g(\hat{q}_i \cdot
#' \widehat{\textrm{PIF}}_{1,i})\Bigg)
#' }
#' and a second `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{2} = h^{-1}\Bigg(\sum\limits_{j=1}^{M_2} h(\hat{w}_j \cdot
#' \widehat{\textrm{PIF}}_{2,j})\Bigg)
#' }
#' This function computes the covariance between the first impact fraction and the
#' weights of the second fraction. This returns a vector where the `j`th entry
#' is given by:
#' \deqn{
#' \operatorname{Cov}(\widehat{\textrm{PIF}}_A, \hat{w}_j) \approx
#' \frac{1}{g'\big(\widehat{\textrm{PIF}}_A\big)}\sum\limits_{i=1}^{M_1}g'(\hat{q}_i
#' \cdot \widehat{\textrm{PIF}}_{A,i})
#' \Bigg[ \hat{q}_i \operatorname{Cov}\Big(
#' \widehat{\textrm{PIF}}_{A,i}, \hat{w}_j\Big) +
#' \widehat{\textrm{PIF}}_{A,i} \operatorname{Cov}\Big( \hat{q}_i,
#' \hat{w}_j\Big) \Bigg]
#' }
#' where \eqn{\widehat{\textrm{PIF}}_{1,:} = (\widehat{\textrm{PIF}}_{1,1},
#' \widehat{\textrm{PIF}}_{1,2}, \dots, \widehat{\textrm{PIF}}_{1,M_1})^{\top}}
#'
#' @note The model currently works under the assumption that all `pif_atomic_class`
#' are independent from weights. Hence will return `0` if `pif2` is an atomic pif.
#'
#' @examples
#' \dontrun{
#' p1 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
#'   var_beta = 0, label ="1")
#' p2 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
#'   var_beta = 0, label = "2")
#' var_weights <- matrix(c(0.1, 0.01, -0.03, 0.01, 0.4, 0.02, -0.03, 0.02, 0.01), ncol = 3)
#' pt1 <- paf_total(p1, p2, weights = c(0.1, 0.9), var_weights = var_weights[1:2,1:2],
#'   label ="total")
#'
#' #This gives 0 as the assumption is that atomic pifs don't correlate with anything
#' cov_ensemble_weights(p1, p2, var_weights = NULL, var_pif_weights = NULL)
#'
#' #This returns the correlation of pt1 with its weights
#' cov_ensemble_weights(pt1, pt1, var_weights = NULL, var_pif_weights = NULL)
#'
#' #Changes if it changes to ensemble
#' p3 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
#'   var_beta = 0, label = "3")
#' p4 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
#'   var_beta = 0, label = "4")
#' p5 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
#'   var_beta = 0, label = "5")
#' pe1 <- paf_ensemble(p3, p4, p5, weights = c(0.1, 0.5, 0.4),
#'   var_weights = var_weights, label = "ensemble")
#'
#' #This returns the correlation of pt1 with its weights
#' cov_ensemble_weights(pe1, pe1, var_weights = NULL, var_pif_weights = NULL)
#'
#' #You can also specify var_weights as either a matrix
#' var_weights     <- matrix(c(0.1, 0.2, 0.1, 0.3, 0.4, -0.2), ncol = 2)
#' var_pif_weights <- matrix(c(0.01, 0.03, 0.021, -0.02, 0.004, -0.5), ncol = 2)
#' cov_ensemble_weights(pe1, pt1, var_weights = var_weights,
#'   var_pif_weights = var_pif_weights)
#'
#' #or a covariance structure:
#' var_weights2     <- default_weight_covariance_structure2(pe1, pt1)
#' var_weights2[[pe1@label]][[pt1@label]] <- var_weights
#' var_weights2[[pt1@label]][[pe1@label]] <- var_weights
#'
#' var_pif_weights2 <- default_weight_pif_covariance_structure2(pe1, pt1)
#' var_pif_weights2[[pe1@label]][[pt1@label]] <- var_pif_weights
#' var_pif_weights2[[pt1@label]][[pe1@label]] <- var_pif_weights
#' cov_ensemble_weights(pe1, pt1, var_weights = var_weights2,
#'   var_pif_weights = var_pif_weights2)
#'
#' }
#'
#'
#' @seealso [cov_atomic_pif()], [cov_ensemble_atomic()] [cov_total_pif()]
#' @keywords internal
cov_ensemble_weights <- function(pif1, pif2, var_weights = NULL, var_pif_weights = NULL,
                                 recursive = !is.null(var_pif_weights)){

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

  #Assign default weights if missing
  if (is.numeric(var_weights)){
    var_weights <- as.matrix(var_weights)
  }

  if (identical(var_weights, matrix(0, ncol = 1, nrow = 1))){
    var_weights <- matrix(0, ncol = length(pif2@weights), nrow = length(pif1@coefs))
  }

  if (!is.matrix(var_weights)){
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
  } else {
    if (ncol(var_weights) != length(pif2@weights) || nrow(var_weights) != length(pif1@coefs)){
      cli::cli_abort(
        paste0(
          "Invalid dimensions for `var_weights` matrix. Should be a matrix of ",
          "{length(pif1@coefs)} x {length(pif2@weights)}"
        )
      )
    }
  }

  if (is.numeric(var_pif_weights)){
    var_pif_weights <- as.matrix(var_pif_weights)
  }

  if (identical(var_pif_weights, matrix(0, ncol = 1, nrow = 1))){
    var_pif_weights <- matrix(0, ncol = length(pif2@weights), nrow = length(pif1@coefs))
  }

  if (!is.matrix(var_pif_weights)){
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
  } else {
    if (ncol(var_pif_weights) != length(pif2@weights) || nrow(var_pif_weights) != length(pif1@coefs)){
      cli::cli_abort(
        paste0(
          "Invalid dimensions for `var_pif_weights` matrix. Should be a matrix of ",
          "{length(pif1@coefs)} x {length(pif2@weights)}"
        )
      )
    }
  }

  #Loop deeper to get the covariance Cov(pif_1, w)
  if (recursive && !is.matrix(var_pif_weights)){ #TODO: This is the part that I have not yet checked
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

  #Get the factors involved in covariance
  alabel <- pif1@label
  klabel <- pif2@label

  #Loop applying the transform to each q*PIF
  pif_deriv_transform_vals <- sapply(pif1@weights*pif1@coefs, pif1@pif_deriv_transform)

  #Expand the variances to matrix form if needed
  if (S7::S7_inherits(var_pif_weights, covariance_structure_class)){
    var_pif_w <- var_pif_weights[[klabel]][[alabel]]
  } else {
    var_pif_w <- var_pif_weights
  }

  if (!is.matrix(var_pif_w)){
    var_pif_w <- matrix(var_pif_w, ncol = length(pif2@weights), nrow = length(pif1@coefs))
  }

  #Expand the variances to matrix form if needed
  if (S7::S7_inherits(var_weights, covariance_structure_class)){
    var_w <- var_weights[[klabel]][[alabel]]
  } else {
    var_w <- var_weights
  }

  if (!is.matrix(var_w)){
    var_w <- matrix(var_w, ncol = length(pif2@weights), nrow = length(pif1@weights))
  }

  as.numeric(
    (1 / pif1@pif_deriv_transform(pif1@pif)) * (
      (pif_deriv_transform_vals * pif1@weights) %*% var_pif_w +
        (pif_deriv_transform_vals * pif1@coefs) %*% var_w
    )
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
#' impact fractions in `pif_ensemble` and `pif_atomic`. This refers to
#' the term \eqn{\operatorname{Cov}\Big( \textrm{PIF}_{A,i}, \widehat{\textrm{PIF}}_{B,j}\Big)}
#' in the equation below. If set to `NULL` its automatically calculated.
#'
#' @param var_pif_weights Covariance vector between the weights in  `pif_ensemble` and
#' the `pif_atomic`. This refers to
#' the term \eqn{\operatorname{Cov}\Big( \hat{q}_i,\widehat{\textrm{PIF}}_{B,j}\Big)}
#' in the equation below. If set to `NULL` its automatically calculated.
#'
#' @section Formula:
#' Given a `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{1} =
#' g^{-1}\Bigg(\sum\limits_{i=1}^{M_1} g(\hat{q}_i \cdot \widehat{\textrm{PIF}}_{1,i})\Bigg)
#' }
#' and a second `pif_global_ensemble`:
#' \deqn{
#' \widehat{\textrm{PIF}}_{2} = h^{-1}\Bigg(\sum\limits_{j=1}^{M_2} h(\hat{w}_j \cdot
#' \widehat{\textrm{PIF}}_{2,j})\Bigg)
#' }
#' This function computes the covariance between the first impact fraction and the
#' coefficients of the second fraction by computing:
#' \deqn{
#' \operatorname{Cov}(\widehat{\textrm{PIF}}_A, \widehat{\textrm{PIF}}_{B,j})
#' \approx \frac{1}{g'\big(\widehat{\textrm{PIF}}_A\big)}\sum\limits_{i=1}^{M_1}g'(\hat{q}_i
#' \cdot \widehat{\textrm{PIF}}_{A,i}) \Bigg[ \hat{q}_i \operatorname{Cov}\Big( \textrm{PIF}_{A,i},
#' \widehat{\textrm{PIF}}_{B,j}\Big) +  \widehat{\textrm{PIF}}_{A,i} \operatorname{Cov}\Big( \hat{q}_i,
#' \widehat{\textrm{PIF}}_{B,j}\Big) \Bigg]
#' }
#' where \eqn{\widehat{\textrm{PIF}}_{1,:} = (\widehat{\textrm{PIF}}_{1,1},
#' \widehat{\textrm{PIF}}_{1,2}, \dots, \widehat{\textrm{PIF}}_{1,M_1})^{\top}}
#'
#' @examples
#' \dontrun{
#' #Works for two atomic fractions
#' p1 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
#'   var_beta = diag(c(0.51, 0.3)), label ="1")
#' p2 <- paf(c(0.2, 0.3), beta = c(0.01, 0.14), var_p = diag(c(0.06, 0.052)),
#'   var_beta = diag(c(0.51, 0.3)), label = "2")
#' cov_ensemble_atomic(p1, p2)
#'
#' #Works for fractions with ensembles
#' pt1 <- paf_total(p1, p2, weights = c(0.1, 0.9), var_weights = diag(c(0.6, 0.2)),
#'   label ="total")
#' p3 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
#'   var_beta = 0, label = "3")
#' cov_ensemble_atomic(pt1, p3)
#'
#' #Changes if it changes to ensemble
#' p4 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
#'   var_beta = 0, label = "4")
#' p5 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
#'   var_beta = 0, label = "5")
#' pe1 <- paf_ensemble(p3, p4, p5, weights = c(0.1, 0.5, 0.4),
#'   var_weights = diag(c(0.6, 0.2, 0.1)), label = "ensemble")
#'
#' #This returns the correlation of pt1 with its weights
#' cov_ensemble_atomic(pe1, p1)
#'
#' #You can also specify var_weights as either a matrix
#' var_pif_weights <- matrix(c(0.01, 0.03, 0.021), ncol = 1)
#' cov_ensemble_atomic(pe1, p1, var_pif_weights = var_pif_weights)
#'
#' #or a covariance structure:
#' var_pif_weights2     <- default_weight_covariance_structure2(pe1, p1)
#' var_pif_weights2[[p3@label]][[p1@label]]  <- var_pif_weights[1]
#' var_pif_weights2[[p1@label]][[p3@label]]  <- var_pif_weights[1]
#' var_pif_weights2[[p4@label]][[p1@label]]  <- var_pif_weights[2]
#' var_pif_weights2[[p1@label]][[p4@label]]  <- var_pif_weights[2]
#' var_pif_weights2[[p5@label]][[p1@label]]  <- var_pif_weights[3]
#' var_pif_weights2[[p1@label]][[p5@label]]  <- var_pif_weights[3]
#'
#' cov_ensemble_atomic(pe1, p1, var_pif_weights = var_pif_weights2)
#'
#' }
#'
#'
#'
#' @seealso [cov_ensemble_weights()], [cov_atomic_pif()], [cov_total_pif()].
#' @keywords internal
cov_ensemble_atomic <- function(pif_ensemble, pif_atomic,
                                var_p = NULL, var_beta = NULL,
                                var_pifs = NULL, var_pif_weights = NULL,
                                recursive = !is.null(var_pifs)){

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
                     var_p    = var_p,
                     var_beta = var_beta)
    )
  }

  if (is.numeric(var_pif_weights)){
    var_pif_weights <- as.matrix(var_pif_weights)
  }

  if (identical(var_pif_weights, matrix(0, ncol = 1, nrow = 1))){
    var_pif_weights <- matrix(0, ncol = 1, nrow = length(pif_ensemble@pif_list))
  }


  if (!is.matrix(var_pif_weights)){
    if (is.null(var_pif_weights)){
      var_pif_weights <- default_weight_pif_covariance_structure2(pif_ensemble, pif_atomic)
    } else if (!S7::S7_inherits(var_pif_weights, covariance_structure_class)){
      cli::cli_abort(
        paste0(
          "`var_pif_weights` should be a `covariance_structure_class`. Use ",
          "`covariance_structure2(pif_ensemble, pif_atomic)`",
          "to create a prototype and work from there."
        )
      )
    }
  } else {
    if (ncol(var_pif_weights) != 1 || nrow(var_pif_weights) != length(pif_ensemble@pif_list)){
      cli::cli_abort(
        paste0(
          "`var_pif_weights` should be a `matrix` of {length(pif_ensemble@pif_list)} x 1"
        )
      )
    }
  }

  if (is.numeric(var_pifs)){
    var_pifs <- as.matrix(var_pifs)
  }

  if (identical(var_pifs, matrix(0, ncol = 1, nrow = 1))){
    var_pifs <- matrix(0, ncol = 1, nrow = length(pif_ensemble@pif_list))
  }


  if (!is.matrix(var_pifs)){
    if (is.null(var_pifs)){
      var_pifs <- default_pif_covariance_structure2(pif_ensemble, pif_atomic)
    } else if (!S7::S7_inherits(var_pifs, covariance_structure_class)){
      cli::cli_abort(
        paste0(
          "`var_pifs` should be a `covariance_structure_class`. Use ",
          "`covariance_structure2(pif_ensemble, pif_atomic)`",
          "to create a prototype and work from there."
        )
      )
    }
  } else {
    if (ncol(var_pifs) != 1 || nrow(var_pifs) != length(pif_ensemble@pif_list)){
      cli::cli_abort(
        paste0(
          "`var_pifs` should be a `matrix` of {length(pif_ensemble@pif_list)} x 1"
        )
      )
    }
  }

  #Calculate the inner weights of the pifs
  if (recursive && !is.matrix(var_pifs)){
    #Loop through each fraction in the ensemble and repeat
    for (k in 1:length(pif_ensemble@pif_list)){
      name_pif_k <- names(pif_ensemble@pif_list)[k]
      var_pifs[[name_pif_k]][[pif_atomic@label]] <-
        cov_ensemble_atomic(
          pif_ensemble = pif_ensemble@pif_list[[k]],
          pif_atomic   = pif_atomic,
          var_p        = var_p,
          var_beta     = var_beta,
          var_pifs     = var_pifs,
          var_pif_weights = var_pif_weights,
          recursive    = TRUE)

      #Symmetry assignment
      var_pifs[[pif_atomic@label]][[name_pif_k]]  <- var_pifs[[name_pif_k]][[pif_atomic@label]]

    }
  }

  #Loop through the sum of 1/g'*g'(qpif)*(q cov(pif1, pifj) +pifi*cov(q*pifj)
  alabel          <- pif_atomic@label
  gprime          <- sapply(pif_ensemble@weights*pif_ensemble@coefs, pif_ensemble@pif_deriv_transform)

  if (!is.matrix(var_pifs)){
    var_pif <- subset(var_pifs, cols = alabel, rows = children(pif_ensemble)) |> as.matrix()
  } else {
    var_pif <- var_pifs
  }

  if (!is.matrix(var_pif_weights)){
    var_pif_weight  <- subset(var_pif_weights, cols = alabel, rows = children(pif_ensemble)) |> as.matrix()
  } else {
    var_pif_weight <- var_pif_weights
  }


  as.numeric(
    1/pif_ensemble@pif_deriv_transform(pif_ensemble@pif)*(
      t(gprime * pif_ensemble@weights) %*% as.matrix(var_pif) +
        t(gprime * pif_ensemble@coefs) %*% as.matrix(var_pif_weight)
    )
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
#' \widehat{\textrm{PIF}}_{B}\Big) } is computed as in  [cov_ensemble_weights()]
#' for each weight \eqn{q_j}.
#'
#' @seealso [from_parameters_covariance_p_component()], [cov_atomic_pif()],
#' [cov_ensemble_atomic()], [cov_ensemble_weights()]
#'
#' @keywords internal
cov_total_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                          var_weights = NULL, var_pif_weights = NULL,
                          var_pifs = NULL) {

  #Check they are pifs
  if (!S7::S7_inherits(pif1, pif_class) || !S7::S7_inherits(pif2, pif_class)){
    cli::cli_abort(
      "Variables `pif1` and `pif2` must be of `pif_class`"
    )
  }

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

    #FIXME: Here
    deriv_transformed <- sapply(pif1@weights*pif1@coefs, pif1@pif_deriv_transform)
    return(
      as.numeric(
          (deriv_transformed * pif1@weights) %*% total_cov +
          (deriv_transformed * pif1@coefs) %*% cov_weights
      ) / multiplier
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
#' @inheritParams cov_total_pif
#' @inheritParams cov_ensemble_weights
#' @inheritParams cov_ensemble_atomic
#'
#' @param ... Multiple additional potential impact fraction objects
#' separated by commas.
#'
#'
#' @examples
#' # Get the approximate link_variance of a pif object
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3,
#'               var_p = 0.1, var_beta = 0.2)
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
#' covariance(pif1, pif2, pif3)
#'
#' # The link_covariance between a pif and itself only has the link_variance as entries
#' covariance(pif1, pif1)
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
#' correlation(pif1, pif2, pif3, var_beta = var_beta)
#' @name covcor
NULL

#' @rdname covcor
#' @export
covariance <- S7::new_generic(
  "covariance", "x",
  function(x, ..., var_p = NULL, var_beta = NULL, var_weights = NULL,
           var_pif_weights = NULL, var_pifs = NULL) {
    S7::S7_dispatch()
  }
)
S7::method(covariance,
           S7::new_union(pif_global_ensemble_class,
                         pif_atomic_class)) <- function(x, ...,
                                                        var_p = NULL,
                                                        var_beta = NULL,
                                                        var_weights = NULL,
                                                        var_pif_weights = NULL,
                                                        var_pifs = NULL) {


  #FIXME: Pass var_pifs to covariance
  # Get the list of fractions
  pif_list <- append(list(x), list(...))
  npifs    <- length(pif_list)

  #Get the names of the pif_list
  pif_names <- sapply(pif_list, names)

  #FIXME: Make covariance pass the other structures too
  #FIXME: Check that the names of the pifs given are in var_beta and var_p
  if (!is.null(var_p)){

    if (!is.null(colnames(var_p)) && !all(pif_names %in% colnames(var_p))){
      cli::cli_abort(
        paste0(
          "Covariance for ",
          "{pif_names[which(!(pif_names %in% colnames(var_p)))][1]} ",
          "was not found in `var_p`."
        )
      )
    } else if (is.null(colnames(var_p))){
      if (ncol(var_p) == length(pif_names)){
        colnames(var_p) <- pif_names
      } else {
        cli::cli_abort(
          "No column names were given to `var_p` and they cannot be automatically assigned"
        )
      }
    }

    if (!is.null(rownames(var_p)) && !all(pif_names %in% rownames(var_p))){
        cli::cli_abort(
          paste0(
            "Covariance for ",
            "{pif_names[which(!(pif_names %in% colnames(var_p)))][1]} ",
            "was not found in `var_p`."
          )
        )
    } else if (is.null(rownames(var_p))){
      if (nrow(var_p) == length(pif_names)){
        rownames(var_p) <- pif_names
      } else {
        cli::cli_abort(
          "No row names were given to `var_p` and they cannot be automatically assigned"
        )
      }
    }

    var_p    <- as_covariance_structure(var_p)
  }

  if (!is.null(var_beta)){

    if (!is.null(colnames(var_beta)) && !all(pif_names %in% colnames(var_beta))){
      cli::cli_abort(
        paste0(
          "Covariance for ",
          "{pif_names[which(!(pif_names %in% colnames(var_beta)))][1]} ",
          "was not found in `var_beta`."
        )
      )
    } else if (is.null(colnames(var_beta))){
      if (ncol(var_beta) == length(pif_names)){
        colnames(var_beta) <- pif_names
      } else {
        cli::cli_abort(
          "No column names were given to `var_beta` and they cannot be automatically assigned"
        )
      }
    }

    if (!is.null(rownames(var_beta)) && !all(pif_names %in% rownames(var_beta))){
      cli::cli_abort(
        paste0(
          "Covariance for ",
          "{pif_names[which(!(pif_names %in% colnames(var_beta)))][1]} ",
          "was not found in `var_beta`."
        )
      )
    } else if (is.null(rownames(var_beta))){
      if (nrow(var_beta) == length(pif_names)){
        rownames(var_beta) <- pif_names
      } else {
        cli::cli_abort(
          "No row names were given to `var_beta` and they cannot be automatically assigned"
        )
      }
    }

    var_beta <- as_covariance_structure(var_beta)
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

        cov_mat[i, j] <- cov_total_pif(pif_list[[i]], pif_list[[j]],
                                       var_p = var_p, var_beta = var_beta)
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
  cov_total_pif(x, x)
}

S7::method(variance, cases_class) <- function(x, ...) {
  if (length(list(...)) > 0) {
    cli::cli_warn(
      "Currently this function does not support more than 1 argument. Ignoring the rest."
    )
  }
  x@variance
}

#' @rdname covcor
#' @export
standard_deviation <- S7::new_generic("standard_deviation", "x")
S7::method(standard_deviation, S7::new_union(pif_global_ensemble_class, pif_atomic_class, cases_class)) <- function(x, ...) {
  sqrt(variance(x, ...))
}

#' @rdname covcor
#' @export
correlation <- S7::new_generic(
  "correlation", "x",
  function(x, ..., var_p = NULL, var_beta = NULL) {
    S7::S7_dispatch()
  }
)
S7::method(correlation, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL) {
  cov2cor(
    covariance(x, ...,
        var_p = var_p, var_beta = var_beta
    )
  )
}
