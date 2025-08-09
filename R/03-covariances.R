#' Covariance matrix, correlation matrix, variance and standard deviation
#' for potential impact fractions
#'
#' Computes the covariance (`cov`) or correlation (`cor`) for multiple
#' potential impact fractions and the variance `var` and standard deviation
#' `sd`for a potential impact fraction.
#'
#' @param x A potential impact fraction
#'
#' @param ... Multiple additional potential impact fraction objects
#' separated by commas.
#'
#' @param sigma_p Covariance matrix for the prevalences in all `pif1` and
#' the ones included in `...`.
#'
#' @param sigma_beta Covariance matrix for the parameter `beta` in all `pif1`
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
#' @examples
#' # Get the approximate variance of a pif object
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' var(my_pif)
#'
#' # This is the same as covariance with just 1 pIF
#' cov(my_pif)
#'
#' # Calculate the covariance between 3 fractions with shared relative risk
#' beta <- 0.3
#' sigma_beta <- 0.1
#' pif1 <- pif(0.5, 0.2, beta, sigma_p = 0.5 * (1 - 0.5) / 100, sigma_beta = sigma_beta)
#' pif2 <- pif(0.3, 0.1, beta, sigma_p = 0.3 * (1 - 0.3) / 100, sigma_beta = sigma_beta)
#' pif3 <- pif(0.7, 0.3, beta, sigma_p = 0.7 * (1 - 0.7) / 100, sigma_beta = sigma_beta)
#' cov(pif1, pif2, pif3, independent_beta = FALSE)
#'
#' # The covariance between a pif and itself only has the variance as entries
#' cov(pif1, pif1, independent_beta = FALSE, independent_p = FALSE)
#'
#' # Or if there is a covariance structure between different betas you can specify with
#' # sigma_beta in the covariance
#' betas <- c(1.3, 1.2, 1.27)
#'
#' # Covariance among all betas
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
#' cov(pif1, pif2, pif_list(pif1, pif3), sigma_beta = sigma_beta, quiet = TRUE)
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
S7::method(cov, S7::new_union(pif_list_class, pif_atomic_class)) <- function(x, ..., sigma_p = NULL, sigma_beta = NULL,
                           independent_p = "guess", independent_beta = "guess",
                           quiet = FALSE) {

  # Get the list of fractions
  pif_list_class <- append(list(x), list(...))
  npifs          <- length(pif_list_class)

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

  # Covariance matrix
  if (npifs > 1) {
    cov_mat <- matrix(0, ncol = npifs, nrow = npifs)
    for (i in 1:(npifs - 1)) {
      for (j in (i + 1):npifs) {

        if (!is.null(sigma_p)) {
          nps <-  pif_list_class_apply_1st(x, length, "p")
          sub_sigma_p <- sigma_p[(nps * (i - 1) + 1):(nps * i), (nps * (j - 1) + 1):(nps * j)]
        } else {
          sub_sigma_p <- NULL
        }

        if (!is.null(sigma_beta)) {
          nbetas <- nps <-  pif_list_class_apply_1st(x, length, "beta")
          sub_sigma_beta <- sigma_beta[(nbetas * (i - 1) + 1):(nbetas * i), (nbetas * (j - 1) + 1):(nbetas * j)]
        } else {
          sub_sigma_beta <- NULL
        }


        cov_mat[i, j] <- cov_pif_list_class(pif_list_class[[i]], pif_list_class[[j]],
                                      sigma_p = sub_sigma_p, sigma_beta = sub_sigma_beta,
                                      independent_p = independent_p[i, j],
                                      independent_beta = independent_beta[i, j],
                                      quiet = quiet)
      }
    }

    # Add lower triangle
    cov_mat <- cov_mat + t(cov_mat) + diag(sapply(pif_list_class, var))

    } else {
      cov_mat <- matrix(var(x), ncol = 1, nrow = 1)
    }

  return(cov_mat)
}

#' Estimate the covariance between two pif_atomic objects
#' @inheritParams covcor
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

  # Covariance matrix
  cov_mat <- pif_class_covariance(pif1, pif2,
    sigma_p = sigma_p,
    sigma_beta = sigma_beta,
    independent_p = independent_p,
    independent_beta = independent_beta,
    quiet = quiet
  )

  return(cov_mat)
}

#' Covariance function for a `pif_list_class`
#'
#' Recursively obtains the covariance matrix of a `pif_list_class`
#'
#' @rdname covcor
#' @keywords internal
cov_pif_list_class <- function(pif1, pif2, sigma_p = NULL, sigma_beta = NULL,
                         independent_p = "guess", independent_beta = "guess",
                         quiet = FALSE) {

  # Base case: both are pif_atomic
  if (S7::S7_inherits(pif1, pif_class) && S7::S7_inherits(pif2, pif_class)) {
    return(
      cov_atomic_pif(pif1, pif2, sigma_p = sigma_p, sigma_beta = sigma_beta,
                     independent_p = independent_p, independent_beta = independent_beta,
                     quiet = quiet)
    )
  }

  # Initialize total
  total <- 0

  # If pif1 is a list (regardless of p2's type)
  if (S7::S7_inherits(pif1, pif_list_class)) {
    for (elem in pif1) {
      total <- total + cov_pif_list_class(elem, pif2, sigma_p = sigma_p, sigma_beta = sigma_beta,
                               independent_p = independent_p, independent_beta = independent_beta,
                               quiet = quiet)
    }
    return(total)
  }

  # If p2 is a list (and p1 isn't, from above)
  if (S7::S7_inherits(pif2, pif_list_class)) {
    for (elem in pif2) {
      total <- total + cov_pif_list_class(pif1, elem, sigma_p = sigma_p, sigma_beta = sigma_beta,
                               independent_p = independent_p, independent_beta = independent_beta,
                               quiet = quiet)
    }
    return(total)
  }

  # If we get here, unsupported types
  cli::cli_abort(
    "Unsupported types for covariance calculation"
  )
}

#' @rdname covcor
#' @export
var <- S7::new_generic("var", "x")
S7::method(var, S7::new_union(pif_list_class, pif_atomic_class)) <- function(x, ...) {
  if (length(list(...)) > 0) {
    cli::cli_alert_danger(
      "Currently this function does not support more than 1 argument. Ignoring the rest."
    )
  }
  cov_pif_list_class(x, x)
}

#' @rdname covcor
#' @export
sd <- S7::new_generic("sd", "x")
S7::method(sd, S7::new_union(pif_list_class, pif_atomic_class)) <- function(x, ...) {
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
S7::method(cor, S7::new_union(pif_list_class, pif_atomic_class)) <- function(x, ..., sigma_p = NULL, sigma_beta = NULL,
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
