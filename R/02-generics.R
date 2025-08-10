#' Print a `pif_class`
#'
#' Prints a `pif_class` object.
#'
#' @param x A `pif_class`
#'
#' @param accuracy The accuracy parameter for [`scales::percent`].
#'
#' @return Called for its side-effects of printing to the console
#'
#' @keywords internal
print_pif_class <- function(x, accuracy){

  # Printed text looks like:
  #   ── Potential Impact Fraction ──
  #
  # PIF = 4.421% [95% CI: 0.179% to 54.364%]
  # standard_deviation(pif %) = 7.003
  # standard_deviation(link(pif)) = 1.657

  pif_val <- scales::percent(x@pif, accuracy = accuracy)
  cilow   <- scales::percent(x@ci[1], accuracy = accuracy)
  cihigh  <- scales::percent(x@ci[2], accuracy = accuracy)
  title   <- ifelse(x@type == "PIF", "Potential Impact Fraction", "Population Attributable Fraction")

  cli::cli_h2(title)
  cli::cli_text(
    "{x@type} = {pif_val} ",
    "[{.emph {scales::percent(x@conf_level)} CI}: {cilow} to {cihigh}]"
  )
  cli::cli_text(
    "standard_deviation({tolower(x@type)} %) = {scales::comma(100*sqrt(x@variance), accuracy = accuracy)}"
  )
  cli::cli_text(
    "standard_deviation(link({tolower(x@type)})) = {scales::comma(sqrt(x@link_variance), accuracy = accuracy)}"
  )

  return(invisible())
}


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
#' my_pif <- pif(p = 0.2, beta = 1.3, var_beta = 0.1)
#' print(my_pif)
#'
#' # Change the ammount of digits to show just 1
#' print(my_pif, accuracy = 0.1)
#' @name print
#' @export
S7::method(print, pif_class) <- function(x, ..., accuracy = 0.001) {
  print_pif_class(x, accuracy)
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
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
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
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' #Default 95% CI
#' confint(my_pif)
#'
#' #Custom 90% ci:
#' confint(my_pif, level = 0.90)
#' @name confint
#' @export
S7::method(confint, pif_class) <- function(object, ..., level = 0.95) {
  #Set the level
  object@conf_level <- level
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
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' summary(my_pif)
#' @name summary
#' @export
S7::method(summary, pif_class) <- function(object, level = 0.95, ...) {
  conf_interval <- confint(object, level = level)

  #Build the return vector
  return_vec <- c("value"      = coef(object),
                  "standard_deviation" = standard_deviation(object),
                  "ci_low"     = conf_interval[1],
                  "ci_up"      = conf_interval[2],
                  "confidence" = level)

  #Assign the name
  names(return_vec)[1] <- fraction_type(object)

  return(return_vec)
}

#' Transform a pif object into a data.frame
#'
#' Gets the potential impact fraction value, the link_variance and the confidence
#' interval values
#'
#' @param x A `pif_class` object.
#' @param ... Additional parameters (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' as.data.frame(my_pif)
#' @name as.data.frame
#' @export
S7::method(as.data.frame, pif_class) <- function(x, ..., level = 0.95) {
  as.data.frame(t(summary(x, level = level)))
}


