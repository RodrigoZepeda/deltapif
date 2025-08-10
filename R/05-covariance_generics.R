#' Get the colink_variances between two pif classes (atomic pifs)
#'
#' @param pif1 A `pif_class` object
#' @param pif2 A `pif_class` object
#' @param sigma_p Colink_variance matrix for the prevalences in both `pif1` and `pif2`.
#' @param sigma_beta Colink_variance matrix for the parameter `beta` in both `pif1`
#' and `pif2`.
#' @param independent_p If `pif1` and `pif2` share the same prevalence data. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param independent_beta If `pif1` and `pif2` share the same `beta` parameter. Either
#' `TRUE`, `FALSE` or `guess` (default).
#' @param quiet Whether to throw warnings or other messages
#'
#' @section Colink_variance matrices:
#' By default if no `sigma_p` is specified this assumes the colink_variances
#' between the parameters `p` of `pif1` and `pif2` are uncorrelated.
#' However, if `pif1` and `pif2` share the same prevalence estimates
#' (i.e. share the same `p`s, then the user should set `independent_p = TRUE`
#' to account for that correlation).
#'
#' The same thing happens with `sigma_beta`. If no `sigma_beta` is specified
#' then the assumption is that the `beta` parameters from `pif1` and
#' from `pif2` are uncorrelated unless `independent_beta` is set to `TRUE`.
#'
#' If the user provides a colink_variance matrix for `sigma_p` then `independent_p`
#' is disregarded. Similarly, if the user provides `sigma_beta` then
#' `independent_beta` is ignored.
#'
#' @keywords internal
pif_class_atomic_variance <- function(pif1, pif2, sigma_p = NULL, sigma_beta = NULL,
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
        "pif1 has {.code length(pif1@p) = {length(pif1@p)}} exposure",
        "categories and pif2 has {.code length(pif2@p) = {length(pif2@p)}}.",
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
        "pif1 has {.code length(pif1@p_cft) = {length(pif1@p_cft)}} counterfactual exposure",
        "categories and pif2 has {.code length(pif2@p_cft) = {length(pif2@p_cft)}}.",
        "They don't seem to represent the same population. They need to ",
        "be of the same length. Are you missing some categories with 0% ",
        "prevalence perhaps?"
      )
    )
  }

  # Check that the ps don't appear similar
  if (independent_p == "guess" && is.null(sigma_p) && all(pif1@p == pif2@p) && all(pif1@sigma_p == pif2@sigma_p)) {
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

  if (independent_beta == "guess" && is.null(sigma_beta) && all(pif1@beta == pif2@beta) && all(pif1@sigma_beta == pif2@sigma_beta)) {
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

  # If p's are the same set sigma_p as the colink_variance, otherwise assume independence
  if (is.null(sigma_p) && (independent_p == "guess" || independent_p)) {
    sigma_p <- matrix(0, ncol = length(pif1@p), nrow = length(pif1@p))
  } else if (is.null(sigma_p)) {
    sigma_p <- pif1@sigma_p
  }

  # If beta's are the same set sigma_beta as the colink_variance, otherwise assume independence
  if (is.null(sigma_beta) && (independent_beta == "guess" || independent_beta)) {
    sigma_beta <- matrix(0, ncol = length(pif1@beta), nrow = length(pif1@beta))
  } else if (is.null(sigma_beta)) {
    sigma_beta <- pif1@sigma_beta
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
    sigma_p = sigma_p,
    sigma_beta = sigma_beta,
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
cov_atomic_pif <- function(pif1, pif2, sigma_p = NULL, sigma_beta = NULL,
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
                                       sigma_p = sigma_p,
                                       sigma_beta = sigma_beta,
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
cov_total_pif <- function(pif1, pif2, sigma_p = NULL, sigma_beta = NULL,
                          independent_p = "guess", independent_beta = "guess",
                          quiet = FALSE) {

  # Base case: both are pif_atomic
  if (S7::S7_inherits(pif1, pif_atomic_class) && S7::S7_inherits(pif2, pif_atomic_class)) {
    return(
      cov_atomic_pif(pif1 = pif1, pif2 = pif2, sigma_p = sigma_p,
                     sigma_beta = sigma_beta, independent_p = independent_p,
                     independent_beta = independent_beta, quiet = quiet)
    )
  }

  # Initialize total
  total <- 0

  # If pif1 is a list (regardless of p2's type)
  if (S7::S7_inherits(pif1, pif_total_class)) {
    for (i in seq_along(pif1@pif_list)) {
      #Compute sum q*cov(pif1, pif2)
      total <- total +
        pif1@weights[[i]]*cov_total_pif(pif1 = pif1@pif_list[[i]],
                                        pif2 = pif2,
                                        sigma_p = sigma_p,
                                        sigma_beta = sigma_beta,
                                        independent_p = independent_p,
                                        independent_beta = independent_beta,
                                        quiet = quiet)
    }
    return(total)
  }

  # If p2 is a list (and p1 isn't, from above)
  if (S7::S7_inherits(pif2, pif_total_class)) {
    for (i in seq_along(pif2@pif_list)) {
      total <- total +
        pif2@weights[[i]]*cov_total_pif(pif1 = pif1,
                                        pif2 = pif2@pif_list[[i]],
                                        sigma_p = sigma_p,
                                        sigma_beta = sigma_beta,
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
#' Computes the colink_variance (`cov`) or correlation (`cor`) for multiple
#' potential impact fractions and the link_variance `var` and standard deviation
#' `sd`for a potential impact fraction.
#'
#' @param x A potential impact fraction
#'
#' @param ... Multiple additional potential impact fraction objects
#' separated by commas.
#'
#' @param sigma_p Colink_variance matrix for the prevalences in all `pif1` and
#' the ones included in `...`.
#'
#' @param sigma_beta Colink_variance matrix for the parameter `beta` in all `pif1`
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
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' var(my_pif)
#'
#' # This is the same as colink_variance with just 1 pIF
#' cov(my_pif)
#'
#' # Calculate the colink_variance between 3 fractions with shared relative risk
#' beta <- 0.3
#' sigma_beta <- 0.1
#' pif1 <- pif(0.5, 0.2, beta, sigma_p = 0.5 * (1 - 0.5) / 100, sigma_beta = sigma_beta)
#' pif2 <- pif(0.3, 0.1, beta, sigma_p = 0.3 * (1 - 0.3) / 100, sigma_beta = sigma_beta)
#' pif3 <- pif(0.7, 0.3, beta, sigma_p = 0.7 * (1 - 0.7) / 100, sigma_beta = sigma_beta)
#' cov(pif1, pif2, pif3, independent_beta = FALSE)
#'
#' # The colink_variance between a pif and itself only has the link_variance as entries
#' cov(pif1, pif1, independent_beta = FALSE, independent_p = FALSE)
#'
#' # Or if there is a colink_variance structure between different betas you can specify with
#' # sigma_beta in the colink_variance
#' betas <- c(1.3, 1.2, 1.27)
#'
#' # Colink_variance among all betas
#' sigma_beta <- matrix(c(
#'   1.0000000, -0.12123053, 0.35429369,
#'   -0.1212305, 1.00000000, -0.04266409,
#'   0.3542937, -0.04266409, 1.00000000
#' ), byrow = TRUE, ncol = 3)
#' pif1 <- pif(0.5, 0.2, betas[1], sigma_p = 0.5 * (1 - 0.5) / 100, sigma_beta = sigma_beta[1, 1])
#' pif2 <- pif(0.3, 0.1, betas[2], sigma_p = 0.3 * (1 - 0.3) / 100, sigma_beta = sigma_beta[2, 2])
#' pif3 <- pif(0.3, 0.1, betas[3], sigma_p = 0.3 * (1 - 0.3) / 100, sigma_beta = sigma_beta[3, 3])
#' cov(pif1, pif2, pif3, sigma_beta = sigma_beta)
#'
#' # Compute the correlation
#' cor(pif1, pif2, pif3, sigma_beta = sigma_beta, quiet = TRUE)
#' @name covcor
NULL

#' @rdname covcor
#' @export
cov <- S7::new_generic(
  "cov", "x",
  function(x, ..., sigma_p = NULL, sigma_beta = NULL,
           independent_p = "guess", independent_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(cov, S7::new_union(pif_total_class, pif_atomic_class)) <- function(x, ..., sigma_p = NULL, sigma_beta = NULL,
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

  if (is.matrix(sigma_beta) && (ncol(sigma_beta) != npifs || nrow(sigma_beta) != npifs)){
    cli::cli_abort(
      paste0(
        "`sigma_beta` has dimensions {nrow(sigma_beta)} x {ncol(sigma_beta)} ",
        "but should be {npifs} x {npifs}"
      )
    )
  }

  if (is.matrix(sigma_p) && (ncol(sigma_p) != npifs || nrow(sigma_p) != npifs)){
    cli::cli_abort(
      paste0(
        "`sigma_p` has dimensions {nrow(sigma_p)} x {ncol(sigma_p)} ",
        "but should be {npifs} x {npifs}"
      )
    )
  }

  # Colink_variance matrix
  if (npifs > 1) {
    cov_mat <- matrix(0, ncol = npifs, nrow = npifs)
    for (i in 1:(npifs - 1)) {
      for (j in (i + 1):npifs) {

        if (!is.null(sigma_p)) {
          nps <-  pif_class_apply_1st(x, length, "p")
          sub_sigma_p <- sigma_p[(nps * (i - 1) + 1):(nps * i), (nps * (j - 1) + 1):(nps * j)]
        } else {
          sub_sigma_p <- NULL
        }

        if (!is.null(sigma_beta)) {
          nbetas <- nps <-  pif_class_apply_1st(x, length, "beta")
          sub_sigma_beta <- sigma_beta[(nbetas * (i - 1) + 1):(nbetas * i), (nbetas * (j - 1) + 1):(nbetas * j)]
        } else {
          sub_sigma_beta <- NULL
        }


        cov_mat[i, j] <- cov_total_pif(pif_list[[i]], pif_list[[j]],
                                       sigma_p = sub_sigma_p, sigma_beta = sub_sigma_beta,
                                       independent_p = independent_p[i, j],
                                       independent_beta = independent_beta[i, j],
                                       quiet = quiet)
      }
    }

    # Add lower triangle
    cov_mat <- cov_mat + t(cov_mat) + diag(sapply(pif_list, var))

  } else {
    cov_mat <- matrix(var(x), ncol = 1, nrow = 1)
  }

  return(cov_mat)
}

#' @rdname covcor
#' @export
var <- S7::new_generic("var", "x")
S7::method(var, S7::new_union(pif_total_class, pif_atomic_class)) <- function(x, ...) {
  if (length(list(...)) > 0) {
    cli::cli_alert_danger(
      "Currently this function does not support more than 1 argument. Ignoring the rest."
    )
  }
  cov_total_pif(x, x)
}

#' @rdname covcor
#' @export
sd <- S7::new_generic("sd", "x")
S7::method(sd, S7::new_union(pif_total_class, pif_atomic_class)) <- function(x, ...) {
  sqrt(var(x, ...))
}

#' @rdname covcor
#' @export
cor <- S7::new_generic(
  "cor", "x",
  function(x, ..., sigma_p = NULL, sigma_beta = NULL,
           independent_p = "guess", independent_beta = "guess",
           quiet = FALSE) {
    S7::S7_dispatch()
  }
)
S7::method(cor, S7::new_union(pif_total_class, pif_atomic_class)) <- function(x, ..., sigma_p = NULL, sigma_beta = NULL,
                                                                              independent_p = "guess", independent_beta = "guess",
                                                                              quiet = FALSE) {
  cov2cor(
    cov(x, ...,
        sigma_p = sigma_p, sigma_beta = sigma_beta,
        independent_p = independent_p, independent_beta = independent_beta,
        quiet = quiet
    )
  )
}
