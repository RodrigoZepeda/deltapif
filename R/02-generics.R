#' Print or show a potential impact fraction
#'
#' Function to print or show a potential impact fraction object
#'
#' @param x A `pif_class`
#'
#' @param ... Additional arguments to pass to `print` or `show`.
#'
#' @param accuracy The accuracy of the printed value
#'
#' @examples
#' my_pif <- new_pif(p = 0.2, beta = 1.3, sigma_beta = 0.1)
#' print(my_pif)
#'
#' # Change the ammount of digits to show just 1
#' print(my_pif, accuracy = 0.1)
#' @name print
#' @export
S7::method(print, pif_class) <- function(x, ..., accuracy = 0.001) {
  # Putting this value outside otherwise package build throws warning of not using scales
  pif_val <- scales::percent(x@pif, accuracy = accuracy)
  title     <- ifelse(x@type == "PIF", "Potential Impact Fraction", "Population Attributable Fraction")
  cli::cli_h2(title)
  cli::cli_text(
    "{x@type} = {pif_val} ",
    "[{.emph {scales::percent(x@conf_level)} CI}: ",
    "{scales::percent(x@ci[1], accuracy = accuracy)} to  ",
    "{scales::percent(x@ci[2], accuracy = accuracy)}]"
  )
  cli::cli_text(
    "Var({x@type} %) = {scales::comma(100^2*x@variance, accuracy = accuracy)}"
  )
  # cli::cli_h3("Parameters:")
  # cli::cli_ul()
  # cli::cli_li("Observed prevalence ({.code p_obs}) = {x@p}")
  # cli::cli_li("Counterfactual prevalence ({.code p_cft}) = {x@p_cft}")
  # cli::cli_li("Relative risks ({.code rr}) = {x@rr}")
  # cli::cli_end()
}


#' Extract coefficients of a pif object
#'
#' Gets the potential impact fraction value
#'
#' @param object A `pif_class` object.
#' @param ... Additional parameters to pass to `coef` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' coef(my_pif)
#' @name coef
#' @export
S7::method(coef, pif_class) <- function(object, ...) {
  object@pif
}

#' Extract confidence intervals of a pif object
#'
#' Gets the confidence interval for the potential impact fraction
#'
#' @param object A `pif_class` object.
#' @param level Level of confidence desired.
#' @param ... Additional parameters to pass to `confint` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' confint(my_pif)
#' @name confint
#' @export
S7::method(confint, pif_class) <- function(object, ...) {
  #Set the level
  #object@conf_level <- level
  return(object@ci)

}


#' Summary of a pif object
#'
#' Gets the potential impact fraction summary
#'
#' @param object A `pif_class` object.
#' @param ... Additional parameters to pass to `summary` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' summary(my_pif)
#' @name summary
#' @export
S7::method(summary, pif_class) <- function(object, level = 0.95, ...) {
  conf_interval <- confint(object, level = level)
  return(
    c("pif"        = coef(object),
      "sd"         = sd(object),
      "ci_low"     = conf_interval[1],
      "ci_up"      = conf_interval[2],
      "confidence" = level)
  )
}

#' Transform a pif object into a data.frame
#'
#' Gets the potential impact fraction value, the variance and the confidence
#' interval values
#'
#' @param x A `pif_class` object.
#' @param ... Additional parameters (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, sigma_p = 0.1, sigma_beta = 0.2)
#' as.data.frame(my_pif)
#' @name as.data.frame
#' @export
S7::method(as.data.frame, pif_class) <- function(x, ...) {
  as.data.frame(t(summary(x)))
}


