#' Get the colink_variances between two pif classes (atomic pifs)
#'
#' @param pif1 A `pif_class` object
#' @param pif2 A `pif_class` object
#' @param var_p Colink_variance matrix for the prevalences in both `pif1` and `pif2`.
#' @param var_beta Colink_variance matrix for the parameter `beta` in both `pif1`
#' and `pif2`.
#' @param independent_p If `pif1` and `pif2` share the same prevalence data. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param independent_beta If `pif1` and `pif2` share the same `beta` parameter. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param quiet Whether to throw warnings or other messages
#'
#' @section Colink_variance matrices:
#' By default if no `var_p` is specified this assumes the colink_variances
#' between the parameters `p` of `pif1` and `pif2` are uncorrelated.
#' However, if `pif1` and `pif2` share the same prevalence estimates
#' (i.e. share the same `p`s, then the user should set `independent_p = TRUE`
#' to account for that correlation).
#'
#' The same thing happens with `var_beta`. If no `var_beta` is specified
#' then the assumption is that the `beta` parameters from `pif1` and
#' from `pif2` are uncorrelated unless `independent_beta` is set to `TRUE`.
#'
#' If the user provides a colink_variance matrix for `var_p` then `independent_p`
#' is disregarded. Similarly, if the user provides `var_beta` then
#' `independent_beta` is ignored.
#'
#' @keywords internal
pif_class_atomic_variance <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                                      independent_p = "guess", independent_beta = "guess",
                                      quiet = FALSE) {
  # Check they are pif objects
  if (!S7::S7_inherits(pif1, pif_atomic_class) ||
      !S7::S7_inherits(pif2, pif_atomic_class)) {
    cli::cli_abort(
      paste0(
        "One of the arguments passed to colink_variance estimation is not a ",
        "`PIFCI::pif_atomic_class` object."
      )
    )
  }

  # Check that both p's have the same length
  if (length(pif1@p) != length(pif2@p)) {
    cli::cli_abort(
      paste0(
        "pif1 has {.code length(pif1@p) = {length(pif1@p)}} exposure ",
        "categories and pif2 has {.code length(pif2@p) = {length(pif2@p)}}. ",
        "They don't seem to represent the same population. They need to ",
        "be of the same length. Are you missing some categories with 0% ",
        "prevalence perhaps?"
      )
    )
  }

  # Check that both p's have the same length
  if (length(pif1@p_cft) != length(pif2@p_cft)) {
    cli::cli_abort(
      paste0(
        "pif1 has {.code length(pif1@p_cft) = {length(pif1@p_cft)}} counterfactual exposure ",
        "categories and pif2 has {.code length(pif2@p_cft) = {length(pif2@p_cft)}}. ",
        "They don't seem to represent the same population. They need to ",
        "be of the same length. Are you missing some categories with 0% ",
        "prevalence perhaps?"
      )
    )
  }

  # Check that the ps don't appear similar
  if (independent_p == "guess" && is.null(var_p) && all(pif1@p == pif2@p) && all(pif1@var_p == pif2@var_p)) {
    independent_p <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The prevalence parameters p for the potential impact fractions appear ",
          "to be the same. If they are, set `independent_p = FALSE`. Otherwise ",
          "set `independent_p = TRUE"
        )
      )
    }
  }

  if (independent_beta == "guess" && is.null(var_beta) && all(pif1@beta == pif2@beta) && all(pif1@var_beta == pif2@var_beta)) {
    independent_beta <- FALSE
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "The relative risk beta parameters for the potential impact fractions appear ",
          "to be the same. If they are, set `independent_beta = FALSE`. Otherwise ",
          "set `independent_beta = TRUE"
        )
      )
    }
  }

  # If p's are the same set var_p as the colink_variance, otherwise assume independence
  if (is.null(var_p) && (independent_p == "guess" || independent_p)) {
    var_p <- matrix(0, ncol = length(pif1@p), nrow = length(pif1@p))
  } else if (is.null(var_p)) {
    var_p <- pif1@var_p
  }

  # If beta's are the same set var_beta as the colink_variance, otherwise assume independence
  if (is.null(var_beta) && (independent_beta == "guess" || independent_beta)) {
    var_beta <- matrix(0, ncol = length(pif1@beta), nrow = length(pif1@beta))
  } else if (is.null(var_beta)) {
    var_beta <- pif1@var_beta
  }

  from_parameters_pif_covariance(
    p1 = pif1@p,
    p2 = pif2@p,
    p1_cft = pif1@p_cft,
    p2_cft = pif2@p_cft,
    rr1 = pif1@rr,
    rr2 = pif2@rr,
    mu_obs1 = pif1@mu_obs,
    mu_obs2 = pif2@mu_obs,
    mu_cft1 = pif1@mu_cft,
    mu_cft2 = pif2@mu_cft,
    var_p = var_p,
    var_beta = var_beta,
    rr_link_deriv_vals1 = pif1@rr_link_deriv_vals,
    rr_link_deriv_vals2 = pif2@rr_link_deriv_vals,
  )
}


#' Estimate the colink_variance between two pif_atomic objects
#'
#' @param pif1 A `pif_atomic_class` object
#' @param pif2 A `pif_atomic_class` object
#'
#' @inheritParams covcor
#'
#' @keywords internal
cov_atomic_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                           independent_p = "guess", independent_beta = "guess",
                           quiet = FALSE) {

  # If the independent ps are not a matrix nor "guess" then assign the same value
  if (!is.matrix(independent_p)) {
    independent_p <- 0
  }

  if (!is.matrix(independent_beta)) {
    independent_beta <- 0
  }

  # Colink_variance matrix
  cov_mat <- pif_class_atomic_variance(pif1, pif2,
                                       var_p = var_p,
                                       var_beta = var_beta,
                                       independent_p = independent_p,
                                       independent_beta = independent_beta,
                                       quiet = quiet)

  return(cov_mat)
}

#' Colink_variance function for a `pif_total_class`
#'
#' Recursively obtains the colink_variance matrix of a `pif_total_class`
#'
#' @param pif1 Either a `pif_atomic_class` or a `pif_total_class`
#' @param pif2 Either a `pif_atomic_class` or a `pif_total_class`
#'
#' @rdname covcor
#' @keywords internal
cov_total_pif <- function(pif1, pif2, var_p = NULL, var_beta = NULL,
                          independent_p = "guess", independent_beta = "guess",
                          quiet = FALSE) {

  # Base case: both are pif_atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_atomic_class)) {
    return(
      cov_atomic_pif(pif1 = pif1, pif2 = pif2, var_p = var_p,
                     var_beta = var_beta, independent_p = independent_p,
                     independent_beta = independent_beta, quiet = quiet)
    )
  }

  # Initialize total
  total <- 0

  # If pif1 is a list (regardless of p2's type)
  if (S7::S7_inherits(pif1, pif_global_ensemble_class)) {
    for (i in seq_along(pif1@pif_list)) {
      #Compute sum q*covariance(pif1, pif2)
      total <- total +
        pif1@pif_weights[[i]]*cov_total_pif(pif1 = pif1@pif_list[[i]],
                                        pif2 = pif2,
                                        var_p = var_p,
                                        var_beta = var_beta,
                                        independent_p = independent_p,
                                        independent_beta = independent_beta,
                                        quiet = quiet)
    }
    return(total)
  }

  # If p2 is a list (and p1 isn't, from above)
  if (S7::S7_inherits(pif2, pif_global_ensemble_class)) {
    for (i in seq_along(pif2@pif_list)) {
      total <- total +
        pif2@pif_weights[[i]]*cov_total_pif(pif1 = pif1,
                                        pif2 = pif2@pif_list[[i]],
                                        var_p = var_p,
                                        var_beta = var_beta,
                                        independent_p = independent_p,
                                        independent_beta = independent_beta,
                                        quiet = quiet)
    }
    return(total)
  }

  # If we get here, unsupported types
  cli::cli_abort(
    "Unsupported types for covariance calculation"
  )
}


#' Colink_variance matrix, correlation matrix, link_variance and standard deviation
#' for potential impact fractions
#'
#' Computes the colink_variance (`covariance`) or correlation (`correlation`) for multiple
#' potential impact fractions and the link_variance `variance` and standard deviation
#' `standard_deviation`for a potential impact fraction.
#'
#' @param x A potential impact fraction
#'
#' @param ... Multiple additional potential impact fraction objects
#' separated by commas.
#'
#' @param var_p Colink_variance matrix for the prevalences in all `pif1` and
#' the ones included in `...`.
#'
#' @param var_beta Colink_variance matrix for the parameter `beta` in all `pif1`
#' and the ones included in `...`.
#'
#' @param independent_p If all the pifs share the same prevalence data. Either
#' `TRUE`, `FALSE`, `guess` (default) or a matrix. If a matrix is given then
#' `independent_p[i,j] = 1` if the i-th and j-th pifs share the same prevalence data
#' `independent_p[i,j] = 0` if the i-th and j-th pifs don't share the same prevalence data.
#'
#' @param independent_beta If all the pifs share the same `beta` parameter. Either
#' `TRUE`, `FALSE` or `guess` (default) or a matrix. If a matrix is given then
#' `independent_beta[i,j] = 1` if the i-th and j-th pifs share the same relative risk parameters
#' `independent_beta[i,j] = 0` if the i-th and j-th pifs don't share the same relative risk parameters.
#'
#' @param quiet Whether to throw warnings and other messages
#'
#' @examples
#' # Get the approximate link_variance of a pif object
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' variance(my_pif)
#'
#' # This is the same as colink_variance with just 1 pIF
#' covariance(my_pif)
#'
#' # Calculate the colink_variance between 3 fractions with shared relative risk
#' beta <- 0.3
#' var_beta <- 0.1
#' pif1 <- pif(0.5, 0.2, beta, var_p = 0.5 * (1 - 0.5) / 100, var_beta = var_beta)
#' pif2 <- pif(0.3, 0.1, beta, var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta)
#' pif3 <- pif(0.7, 0.3, beta, var_p = 0.7 * (1 - 0.7) / 100, var_beta = var_beta)
#' covariance(pif1, pif2, pif3, independent_beta = FALSE)
#'
#' # The colink_variance between a pif and itself only has the link_variance as entries
#' covariance(pif1, pif1, independent_beta = FALSE, independent_p = FALSE)
#'
#' # Or if there is a colink_variance structure between different betas you can specify with
#' # var_beta in the colink_variance
#' betas <- c(1.3, 1.2, 1.27)
#'
#' # Colink_variance among all betas
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
           independent_p = "guess", independent_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(covariance, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL,
                                                                              independent_p = "guess", independent_beta = "guess",
                                                                              quiet = FALSE) {

  # Get the list of fractions
  pif_list <- append(list(x), list(...))
  npifs    <- length(pif_list)

  # If the independent ps are not a matrix nor "guess" then assign the same value
  if (!is.matrix(independent_p)) {
    independent_p <- matrix(independent_p, ncol = npifs, nrow = npifs)
  }

  if (!is.matrix(independent_beta)) {
    independent_beta <- matrix(independent_beta, ncol = npifs, nrow = npifs)
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

  # Colink_variance matrix
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
                                       independent_p = independent_p[i, j],
                                       independent_beta = independent_beta[i, j],
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
           independent_p = "guess", independent_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(correlation, S7::new_union(pif_global_ensemble_class, pif_atomic_class)) <- function(x, ..., var_p = NULL, var_beta = NULL,
                                                                              independent_p = "guess", independent_beta = "guess",
                                                                              quiet = FALSE) {
  cov2cor(
    covariance(x, ...,
        var_p = var_p, var_beta = var_beta,
        independent_p = independent_p, independent_beta = independent_beta,
        quiet = quiet
    )
  )
}
