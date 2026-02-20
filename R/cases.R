#' Attributable cases
#'
#' Calculates the number of attributable cases or the number of cases
#' that would be averted under a counterfactual scenario for a
#' given fraction (either `paf` or `pif`).
#'
#' @param cases The overall number of cases in the population.
#'
#' @param pif A potential impact fraction object created by `pif`, `paf`,
#' `pif_total`, `pif_ensemble`, `paf_total` or `paf_ensemble`.
#'
#' @param paf A population attributable fraction object created by `paf`,
#'  `paf_total` or `paf_ensemble`.
#'
#' @param link Link function such that the case confidence intervals stay
#' within the expected bounds (either `log` or `identity`).
#'
#' @param variance The estimated variance for the cases (default = 0).
#'
#' @inheritParams pif
#'
#' @returns A `cases_class` object with the attributable cases.
#'
#' @section Formulas:
#'
#' The attributable cases are calculated as:
#' \deqn{
#' \text{Attributable cases} = \textrm{PAF} \times \textrm{Cases}
#' }
#' and the averted cases are respectively:
#' \deqn{
#' \text{Averted cases} = \textrm{PIF} \times \textrm{Cases}
#' }
#'
#' The variance is estimated using the product-variance formula:
#' \deqn{
#' \textrm{Var}[\text{Averted cases}] =
#'  \sigma^2_{\textrm{Cases}} \cdot \big( \textrm{PIF}\big)^2 +
#'  \sigma^2_{\textrm{PIF}} \cdot \big( \textrm{Cases} \big)^2 +
#'  \sigma^2_{\textrm{PIF}} \cdot \sigma^2_{\textrm{Cases}}
#' }
#'
#' @details
#' Negative cases are interpreted as cases that would be caused by
#' the intervention.
#' @examples
#' frac <- paf(p = 0.499, beta = log(3.6), var_p = 0.002, var_beta = FALSE)
#' attributable_cases(100, paf = frac)
#'
#' frac <- pif(p = 0.499, beta = log(3.6), p_cft = 0.1, var_p = 0.002, var_beta = FALSE)
#' averted_cases(100, pif = frac)
#'
#'
#' @seealso [pif()], [paf()]
#'
#' @name casecalc
NULL

#' @rdname casecalc
#' @export
averted_cases <- function(cases, pif, variance = 0, conf_level = 0.95,
                          link = "identity", link_inv = NULL,
                          link_deriv = NULL){

  if (!S7::S7_inherits(pif, pif_class)){
    cli::cli_abort(
      "Object `pif` should be a `pif_class`. Use `pif` or `paf` to create a fraction."
    )
  }

  if (!is.function(link) && !(link %in% c("log", "identity"))){
    cli::cli_abort(
      "Invalid link {link}. Use either 'log'  or 'identity'."
    )
  }

  link_name <- link
  check_links(link = link, link_inv = link_inv, link_deriv = link_deriv)

  # Get the inverse function
  if (is.null(link_inv) & is.character(link)) {
    link_inv <- parse_inv_link(link)
  }

  if (is.null(link_deriv) & is.character(link)){
    link_deriv <- parse_deriv_link(link)
  }

  if (!is.function(link_deriv) && is.null(link_deriv)) {
    link_deriv <- Deriv::Deriv(link)
  }

  # Get the functions
  link    <- parse_link(link)

  cases_class(
    overall_cases = cases,
    pif_obj    = pif,
    variance_cases   = variance,
    conf_level = conf_level,
    link       = link,
    link_inv   = link_inv,
    link_deriv = link_deriv
  )

}

#' @rdname casecalc
#' @export
attributable_cases <- function(cases, paf, variance = 0, conf_level = 0.95,
                          link = "identity", link_inv = NULL,
                          link_deriv = NULL){

  #Conceptual safeguard
  if (paf@type != "PAF"){
    cli::cli_abort(
      "Can only estimate `attributable_cases` with a population attributable fraction (paf). Did you mean `averted_cases`?"
    )
  }

  averted_cases(cases = cases, pif = paf, variance = variance,
                conf_level = conf_level, link = link, link_inv = link_inv,
                link_deriv = link_deriv)

}



