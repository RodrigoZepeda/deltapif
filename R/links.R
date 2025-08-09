#' Link functions
#'
#' A collection of common link functions, for calculating the variance of
#' the  potential impact fraction.
#'
#' @note
#' When used, the variance is calculated for `linkfun(pif)`.
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
  log(pif / (1 - pif))
}

#' @rdname linkfuns
#' @export
log_complement <- function(pif) {
  log(1 - pif)
}

#' @rdname linkfuns
#' @export
hawkins <- function(pif) {
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
  1 / (pif * (1 - pif))
}

#' @rdname deriv_linkfuns
#' @export
deriv_log_complement <- function(pif) {
  1 / (pif - 1)
}

#' @rdname deriv_linkfuns
#' @export
deriv_hawkins <- function(pif) {
  1 / (sqrt(pif^2 + 1))
}
