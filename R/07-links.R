#' Link functions
#'
#' A collection of common link functions, for calculating the link_variance of
#' the  potential impact fraction.
#'
#' @note
#' When used, the link_variance is calculated for `linkfun(pif)`.
#'
#' @param pif The value of a potential impact fraction or a
#' population attributable fraction
#'
#' @details
#' The functions programmed are as follows
#' \deqn{
#'  \text{logit}(\text{PIF}) = \ln\Bigg(\dfrac{
#'    \text{PIF}
#'   }{
#'    1 - \text{PIF}
#'   }\Bigg),
#' }
#' \deqn{
#'  \text{log-complement}(\text{PIF}) = \ln\big(1 - \text{PIF}\big),
#' }
#' and
#' \deqn{
#'  \text{Hawkins}(\text{PIF}) = \ln\Big(\text{PIF} + \sqrt{\text{PIF}^2 + 1}\Big).
#' }
#'
#' @seealso [inv_linkfuns] for their inverses and [deriv_linkfuns] for their
#' derivatives
#'
#' @name linkfuns
NULL

#' @rdname linkfuns
#' @export
logit <- function(pif) {
  if (pif > 1 || pif < 0){
    cli::cli_abort("Invalid value for pif = {round(pif,2)}. PIF must be in (0,1)")
  }
  log(pif / (1 - pif))
}

#' @rdname linkfuns
#' @export
log_complement <- function(pif) {
  if (pif > 1){
    cli::cli_abort("Invalid value for pif = {round(pif,2)} > 1.")
  }
  log(1 - pif)
}

#' @rdname linkfuns
#' @export
hawkins <- function(pif) {
  if (pif > 1){
    cli::cli_alert_warning("Invalid value for pif = {round(pif,2)} > 1.")
  }
  log(pif + sqrt(pif^2 + 1))
}

#' Inverses of link functions
#'
#' A collection of the inverses of the link functions
#' for the potential impact fraction.
#'
#' @param x A value such that `inv_link(x)` is a potential impact fraction
#' or a population attributable fraction.
#'
#' @details
#' The functions programmed are as follows
#' \deqn{
#'  \text{inv\_logit}(\text{PIF}) = \dfrac{
#'    1
#'   }{
#'    1 + \exp(-x)
#'   },
#' }
#' \deqn{
#'  \text{inv\_log-complement}(\text{PIF}) = 1 - \exp(x),
#' }
#' and
#' \deqn{
#'  \text{inv\_Hawkins}(\text{PIF}) =
#'   \frac{1}{2} \exp(-x) \cdot \big(\exp(x) - 1\big) \cdot \big(\exp(x) + 1\big)
#' }
#'
#'
#' @seealso [linkfuns] for the definition of the link functions and
#' [deriv_linkfuns] for their derivatives.
#' @name inv_linkfuns
NULL

#' @rdname inv_linkfuns
#' @export
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

#' @rdname inv_linkfuns
#' @export
inv_log_complement <- function(x) {
  1 - exp(x)
}

#' @rdname inv_linkfuns
#' @export
inv_hawkins <- function(x) {
  0.5 * exp(-x) * (exp(2 * x) - 1)
}

#' Derivatives of link functions
#'
#' A collection of the derivatives of the link functions
#' for the potential impact fraction.
#'
#' @inheritParams linkfuns
#'
#' @details
#' The functions programmed are as follows
#' \deqn{
#'  \text{deriv\_logit}(\text{PIF}) = \dfrac{
#'    1
#'   }{
#'    \text{PIF} \cdot (1 - \text{PIF})
#'   },
#' }
#' \deqn{
#'  \text{deriv\_log-complement}(\text{PIF}) = \dfrac{
#'    1
#'   }{
#'    \text{PIF} - 1
#'   },
#' }
#' and
#' \deqn{
#'  \text{deriv\_Hawkins}(\text{PIF}) = \dfrac{
#'    1
#'   }{
#'    \sqrt{\text{PIF}^2 + 1}
#'   },
#' }
#'
#' @seealso [linkfuns] for the definition of the link functions and [inv_linkfuns]
#' for the inverses of the link functions.
#'
#' @name deriv_linkfuns
NULL

#' @rdname deriv_linkfuns
#' @export
deriv_logit <- function(pif) {
  if (pif > 1 || pif < 0){
    cli::cli_abort("Invalid value for pif = {round(pif,2)}. PIF must be in (0,1)")
  }
  1 / (pif * (1 - pif))
}

#' @rdname deriv_linkfuns
#' @export
deriv_log_complement <- function(pif) {
  if (pif > 1){
    cli::cli_abort("Invalid value for pif = {round(pif,2)} > 1.")
  }
  1 / (pif - 1)
}

#' @rdname deriv_linkfuns
#' @export
deriv_hawkins <- function(pif) {
  if (pif > 1){
    cli::cli_alert_warning("Invalid value for pif = {round(pif,2)} > 1.")
  }
  1 / (sqrt(pif^2 + 1))
}

#' Link parsers
#'
#' Functions to parse the link (or inverse link) from a word to a function
#'
#' @param link_name The name of the link or a function.
#'
#' @details
#' The following are valid link names:
#' \describe{
#'   \item{identity}{The function `f(x) = x`.with inverse `finv(x) = x`}
#'   \item{logit}{The function `f(x) = ln(x / (1 - x))` with inverse `finv(x) = 1 / (1 + exp(-x))`}
#'   \item{log-complement}{The function `f(x) = ln(1 - x)` with inverse `finv(x) = 1 - exp(x)`}
#'   \item{hawkins}{The function `f(x) = ln(x + sqrt(x^2 + 1))` with inverse `finv(x) = 0.5 * exp(-x) * (exp(2 * x) - 1)`}
#'   \item{exponential}{The function `f(x) = exp(x)` with inverse `f(x) = ln(x)`}
#' }
#'
#' @return A function corresponding to the `link_name`.
#'
#' @seealso [linkfuns] and [inv_linkfuns]
#'
#' @note If a function is supplied to `link_name` the same function is returned
#' @name link_parsers


#' @rdname link_parsers
#' @export
parse_link <- function(link_name) {
  if (is.function(link_name)) {
    return(link_name)
  }

  link_name <- tolower(link_name)
  link_name <- gsub("[ _-]+", "", link_name)


  if (link_name == "identity") {
    return(identity)
  } else if (link_name == "logit") {
    return(logit)
  } else if (link_name == "logcomplement") {
    return(log_complement)
  } else if (link_name == "hawkins") {
    return(hawkins)
  } else if (link_name == "exponential" || link_name == "exp") {
    return(exp)
  } else {
    cli::cli_abort(
      paste0(
        "Cannot find link {.val {link_name}}. Please specify ",
        "the function using {.code rr_link}"
      )
    )
  }
}

#' @rdname link_parsers
#' @export
parse_inv_link <- function(link_name) {
  if (is.function(link_name)) {
    return(link_name)
  }

  link_name <- tolower(link_name)
  link_name <- gsub("[ _-]+", "", link_name)


  if (link_name == "identity") {
    return(identity)
  } else if (link_name == "logit") {
    return(inv_logit)
  } else if (link_name == "logcomplement") {
    return(inv_log_complement)
  } else if (link_name == "hawkins") {
    return(inv_hawkins)
  } else if (link_name == "exponential" || link_name == "exp") {
    return(log)
  } else {
    cli::cli_abort(
      paste0(
        "Cannot find link {.val {link_name}}. Please specify the",
        "function using {.code rr_link}"
      )
    )
  }
}
