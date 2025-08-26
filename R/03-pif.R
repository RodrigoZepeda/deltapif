#' Potential Impact fraction and Population Attributable Fraction
#'
#' Calculates the potential impact fraction `pif` or the population
#' attributable fraction `paf` for a categorical exposure considering
#' an observed prevalence of  `p` and a relative risk
#' (or relative risk parameter) of `beta`.
#'
#' @inheritParams classes
#' @param rr_link_deriv Derivative of the link function for the relative risk.
#' The function tries to build it automatically from `rr_link` using
#' [Deriv::Deriv()].
#'
#' @param link_deriv Derivative of the `link` function. The function tries
#' to build it automatically from `link` using [Deriv::Deriv()].
#'
#' @param quiet Whether to show messages.
#'
#' @note This function assumes `p` and `beta` have been pre-computed from
#' the data and the individual-level data are not accessible to the
#' researchers. If either the data for the individual-level prevalence of
#' exposure `p` or the data for the individual-level risk estimate `beta`
#' can be accessed by the researcher other methods (such as the `pifpaf`
#' package should be preferred).
#'
#' @section Formulas:
#' This function computed the potential impact fraction and its confidence
#' intervals using Walter's formula:
#'
#' \deqn{
#' \textrm{PIF} =
#'  \dfrac{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i -
#'    \sum\limits_{i=1}^N p_i^{\text{cft}} \text{RR}_i
#'   }{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i
#'   }, \quad \text{ and } \quad
#' \textrm{PAF} =
#'  \dfrac{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i - 1
#'   }{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i
#'   }
#' }
#'
#' in the case of `N` exposure categories which is equivalent to Levine's
#' formula when there is only `1` exposure category:
#'
#' \deqn{
#' \textrm{PIF} =
#'  \dfrac{
#'    p (\text{RR} - 1) - p^{\text{cft}} (\text{RR} - 1)
#'   }{
#'    1 + p (\text{RR} - 1)
#'   }
#'   \quad \textrm{ and } \quad
#' \textrm{PAF} =
#'  \dfrac{
#'    p (\text{RR} - 1)
#'   }{
#'    1 + p (\text{RR} - 1)
#'   }
#' }
#'
#' The construction of confidence intervals is done via a link function. Its
#' variance is given by:
#'
#' \deqn{
#' \sigma_f^2 = \text{Var}\Big[ f(\textrm{PIF}) \Big] \approx
#' \Big(f'(\textrm{PIF})\Big)^2 \text{Var}\Big[ \textrm{PIF} \Big]
#' }
#' and the intervals are constructed as:
#' \deqn{
#' \text{CI} = f^{-1}\Big(  f(\textrm{PIF})  \pm Z_{\alpha/2} \cdot \sigma_f \Big)
#' }
#'
#' the function \eqn{f} is called the link function for the PIF, \eqn{f^{-1}}
#' represents its inverse, and \eqn{f'} its derivative.
#'
#' @section Link functions for the PIF:
#'
#' By default the `pif` and `paf` calculations use the `log-complement` link
#' which guarantees the impact fractions' intervals have valid values
#' (between \code{-Inf} and `1`).
#'
#' Depending on the application the following link functions are also
#' implemented:
#'
#' \describe{
#'   \item{log-complement}{To achieve fractions between `-Inf` and `1`. This is
#'      the function `f(x) = ln(1 - x)` with inverse `finv(x) = 1 - exp(x)`}
#'   \item{logit}{To achieve strictly positive fractions between `0` and `1`.
#'      This is the function `f(x) = ln(x / (1 - x))` with inverse
#'      `finv(x) = 1 / (1 + exp(-x))`}
#'   \item{identity}{An approximation for fractions that does not guarantee
#'      confidence intervals in valid regions. This is the function
#'      `f(x) = x`.with inverse `finv(x) = x`}
#'   \item{hawkins}{Hawkins' fraction for controlling the variance. This is the
#'      function `f(x) = ln(x + sqrt(x^2 + 1))` with inverse
#'      `finv(x) = 0.5 * exp(-x) * (exp(2 * x) - 1)`}
#'   \item{custom}{User specified (see below).}
#' }
#'
#' In general, `logit` should be preferred if it is known and certain that the
#' fractions can only be positive (i.e. when all relative risks (including CIs)
#' are > 1 and prevalence > 0 and there is an epidemiological / biological
#' justification).
#'
#' **Custom link** functions can be implemented as long as they are invertible
#' in the range of interest by providing:
#' \describe{
#'   \item{link}{The function}
#'   \item{link_inv}{The inverse function of `link`}
#'   \item{link_deriv}{The derivative of `link`}
#' }
#' If no derivative is provided the package attempts to estimate it
#' symbolically using [Deriv::Deriv()] however there is no guarantee that this
#' will work for non-standard functions (i.e. not logarithm / trigonometric /
#' exponential)
#'
#' As an example considering implementing a square root custom link:
#'
#' ```{r custom_link, eval = FALSE}
#' # The link is square root
#' link      <- sqrt
#'
#' # Inverse of sqrt(x) is x^2
#' link_inv  <- function(x) x^2
#'
#' # The derivative of sqrt(x) is 1/(2 * sqrt(x))
#' link_deriv <- function(pif) 1 / (2 * sqrt(pif))
#'
#' # Then the pif can be calculated as:
#' pif(p = 0.499, beta = 1.6, p_cft = 0.499/2, var_p = 0.1, var_beta = 0.2,
#'      link = link, link_inv = link_inv, link_deriv = link_deriv)
#' ```
#'
#' @section Link functions for beta:
#'
#' By default the `pif` and `paf` use the `exponential` link which means that
#' the values for `beta` are the log-relative risks and the
#' variance `var_beta` corresponds to the log-relative risk's variance.
#' Depending on the relative risk's source the following options are available:
#'
#' \describe{
#'   \item{exponential}{For when the relative risks correspond to the
#'          exponential of another parameter, usually called `beta`. This is
#'          the exponential function `f(beta) = exp(beta)` with inverse
#'          `finv(rr) = log(rr) = beta`}
#'   \item{identity}{For when the relative risk and their variance are
#'          reported directly. This is the function `f(beta) = beta`
#'          with inverse `finv(rr) = rr = beta`}
#'   \item{custom}{User specified (see below).}
#' }
#'
#' **Note** that in most cases including contingency tables, Poisson and
#' Logistic regressions, and Cox Proportional Hazards, the relative risks are
#' estimated by exponentiating a parameter.
#'
#' As in the previous section, custom link functions for beta can be implemented
#' as long as they are invertible in the range of interest by providing the
#' function `rr_link` and its derivative `rr_link_deriv`. If no derivative
#' is provided the package does an attempt to estimate it symbolically
#' using [Deriv::Deriv()] however there is no guarantee that this
#' will work non-standard functions (i.e. not logarithm / trigonometric /
#' exponential)
#'
#' @section Population Attributable Fraction:
#'
#' The population attributable fraction corresponds to the potential impact
#' fraction at the theoretical minimum risk level. It is assumed that the
#' theoretical minimum risk level is a relative risk of 1. If no
#' counterfactual prevalence `p_cft` is specified, the model computes
#' the population attributable fraction.
#'
#' @examples
#' #EXAMPLE 1: ONE EXPOSURE CATEGORY (CLASSIC LEVIN)
#' #---------------------------------------------------------------------------
#' # This example comes from Levin 1953
#' # Relative risk of lung cancer given smoking was 3.6
#' # Proportion of individuals smoking where 49.9%
#' # Calculates PAF (i.e. counterfactual is no smoking)
#' paf(p = 0.499, beta = 3.6)
#'
#' # Assuming that beta and p had a link_variance
#' paf(p = 0.499, beta = 3.6, var_p = 0.001, var_beta = 1)
#'
#' # If the link_variance was to high a logistic transform would be required
#' # Generates incorrect values for the interval:
#' paf(p = 0.499, beta = 3.6, var_p = 0.1, var_beta = 3)
#'
#' # Logit fixes it
#' paf(p = 0.499, beta = 3.6, var_p = 0.1, var_beta = 3,
#'     link = "logit", quiet = TRUE)
#'
#' # If the counterfactual was reducing the smoking population by 1/2
#' pif(p = 0.499, beta = 1.6, p_cft = 0.499/2, var_p = 0.001,
#'     var_beta = 1, link = "logit", quiet = TRUE)
#'
#' #EXAMPLE 2: MULTIPLE EXPOSURE CATEGORIES
#' #---------------------------------------------------------------------------
#' #In "Excess Deaths Associated With Underweight, Overweight, and Obesity"
#' #the PAF of mortality associated to different BMI levels is calculated
#'
#' #These are the relative risks for age-group 25-59 for each BMI level:
#' rr <- c("<18.5" = 1.38, "25 to <30" = 0.83 ,
#'         "30 to <35" = 1.20, ">=35" = 1.83)
#'
#' #While the prevalences are:
#' p <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'         "30 to <35" = 17.3, ">=35" = 13.3) / 100
#'
#' #The variance of the relative risk is obtained from the CIs:
#' var_log_rr <- c("<18.5" = 0.2653156, "25 to <30" =  0.1247604,
#'                 "30 to <35" = 0.1828293, ">=35" = 0.1847374)^2
#'
#' #Note that we are omitting the group "18.5 to < 25" as it is the
#' #reference (RR = 1)
#' paf(p = p, beta = rr, var_beta = var_log_rr, var_p = 0, quiet = TRUE,
#'     link = "logit")
#'
#' #We can compute a potential impact fraction of, for example, reducing the
#' #amount of people over 35 to 0 and having them in the "30 to <35":
#' p_cft <- c("<18.5" = 1.9, "25 to <30" = 34.8,
#'             "30 to <35" = 17.3 + 13.3, ">=35" = 0) / 100
#'
#' #The potential impact fraction is as follows:
#' pif(p = p, p_cft = p_cft, beta = rr, link = "logit",
#'   var_beta = var_log_rr, var_p = 0, quiet = TRUE)
#'
#' @name pifpaf
#' @seealso [pif_total()], [pif_ensemble()], [paf_total()], [paf_ensemble()],
#' @md
NULL

#' Population Attributable Fraction
#' @rdname pifpaf
#' @export
paf <- function(p, beta,
                var_p = NULL,
                var_beta = NULL,
                rr_link = "exponential",
                rr_link_deriv = NULL,
                link = "log-complement",
                link_inv = NULL,
                link_deriv = NULL,
                conf_level = 0.95,
                quiet = FALSE,
                label = NULL) {
  pif(
    p = p, p_cft = rep(0, length(p)), beta = beta, var_p = var_p,
    var_beta = var_beta, rr_link = rr_link,
    rr_link_deriv = rr_link_deriv, link = link, link_inv = link_inv,
    link_deriv = link_deriv, conf_level = conf_level, quiet = quiet,
    type = "PAF", label = label
  )
}

#' Potential impact fraction
#' @rdname pifpaf
#' @export
pif <- function(p,
                p_cft         = rep(0, length(p)),
                beta,
                var_p       = NULL,
                var_beta    = NULL,
                rr_link       = "exponential",
                rr_link_deriv = NULL,
                link          = "log-complement",
                link_inv      = NULL,
                link_deriv    = NULL,
                conf_level    = 0.95,
                type          = "PIF",
                quiet         = FALSE,
                label         = NULL) {

  #Check that type is not PAF if p_cft has values
  if (type == "PAF" & any(p_cft > 0)){
    if (!quiet){
      cli::cli_warn(
        paste0(
          "type `PAF` was selected but the counterfactual prevalence `p_cft` ",
          "has non-zero values. Are you sure you are not calculating a PIF ",
          "instead?"
        )
      )
    }
  }

  #Check that I am not taking a covariance matrix for var_p and var_beta
  if (S7::S7_inherits(var_p, covariance_structure_class)){
    var_p <- as.matrix(var_p)
  }

  if (S7::S7_inherits(var_beta, covariance_structure_class)){
    var_beta <- as.matrix(var_beta)
  }

  # Check that is link is given as character then link inv and link_deriv are not specified
  link_name <- link #Variable saved until the pif is calculated
  check_links(link = link, link_inv = link_inv, link_deriv = link_deriv)
  check_rr_links(rr_link = rr_link, rr_link_deriv = rr_link_deriv)

  # Get the inverse function
  if (is.null(link_inv) & is.character(link)) {
    link_inv <- parse_inv_link(link)
  }

  if (is.null(link_deriv) & is.character(link)){
    link_deriv <- parse_deriv_link(link)
  }

  if (is.null(rr_link_deriv) & is.character(rr_link)){
    rr_link_deriv <- parse_deriv_link(rr_link)
  }

  # Get the functions
  rr_link <- parse_link(rr_link)
  link    <- parse_link(link)

  # Get the derivatives of the links
  if (!is.function(rr_link_deriv) && is.null(rr_link_deriv)) {
    rr_link_deriv <- Deriv::Deriv(rr_link)
  }

  if (!is.function(link_deriv) && is.null(link_deriv)) {
    link_deriv <- Deriv::Deriv(link)
  }

  # Get zero matrices for var_p and var_beta if not given
  if (is.null(var_p)) {
    var_p <- matrix(0, ncol = length(p), nrow = length(p))
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "Assuming parameters `p` have no variance Use `var_p` ",
          "to input their link_variances and/or covariance"
        )
      )
    }
  }

  if (length(var_p) == 1 && var_p == 0){
    var_p <- matrix(0, ncol = length(p), nrow = length(p))
  }

  if (is.null(var_beta)) {
    var_beta <- matrix(0, ncol = length(beta), nrow = length(beta))
    if (!quiet){
      cli::cli_alert_warning(
        paste0(
          "Assuming parameters `beta` have no variance Use `var_beta` ",
          "to input their link_variances and/or covariance"
        )
      )
    }
  }

  if (length(var_beta) == 1 &&  var_beta == 0){
    var_beta <- matrix(0, ncol = length(beta), nrow = length(beta))
  }

  # If vectors provided then transform to approximate matrices
  upper_bound_p <- FALSE
  if (is.vector(var_p) && length(var_p) > 1) {
    #upper_bound_p <- TRUE
    var_p_mat <- matrix(0, ncol = length(var_p), nrow = length(var_p))
    for (i in 1:(nrow(var_p_mat) - 1)){
      for (j in (i + 1):ncol(var_p_mat)){
        var_p_mat[i,j] <- sqrt(var_p[i]*var_p[j])
        var_p_mat[j,i] <- var_p_mat[i,j]
      }
    }
    var_p <- var_p_mat + diag(var_p)
    if (!quiet){
      cli::cli_warn(
        paste0(
          "Assuming parameters `p` are correlated but correlation is unknown. ",
          "If they are uncorrelated redefine `var_p = diag(var_p)` to ",
          "transform them into a matrix with no correlations."
        )
      )
    }
  }

  upper_bound_beta <- FALSE
  if (is.vector(var_beta) && length(var_beta) > 1) {
    #upper_bound_beta <- TRUE
    var_beta_mat <- matrix(0, ncol = length(var_beta), nrow = length(var_beta))
    for (i in 1:(nrow(var_beta_mat) - 1)){
      for (j in (i+1):ncol(var_beta_mat)){
        var_beta_mat[i,j] <- sqrt(var_beta[i]*var_beta[j])
        var_beta_mat[j,i] <- var_beta_mat[i,j]
      }
    }
    var_beta <- var_beta_mat + diag(var_beta)
    if (!quiet){
      cli::cli_warn(
        paste0(
          "Assuming parameters `beta` are correlated but correlation is unknown. ",
          "If they are uncorrelated redefine `var_beta = diag(var_beta)` to ",
          "transform them into a matrix with no correlations."
        )
      )
    }
  }

  if (is.null(label)){
    label <- paste0("deltapif-", sub("\\.", "", as.character(abs(stats::rnorm(1)))))
  }

  pif <- pif_atomic_class(
     p          = p,
     p_cft      = p_cft,
     beta       = beta,
     var_p      = var_p,
     var_beta   = var_beta,
     link       = link,
     link_inv   = link_inv,
     link_deriv = link_deriv,
     rr_link    = rr_link,
     label      = label,
     rr_link_deriv = rr_link_deriv,
     conf_level = conf_level,
     type       = type,
     upper_bound_p = upper_bound_p,
     upper_bound_beta = upper_bound_beta
  )

  if (is.character(link_name) && link_name == "logit" && coef(pif) <= 0){
    cli::cli_warn(
      paste0(
        "Value for {fraction_type(pif)} = {round(coef(pif),2)} <= 0. ",
        "Change link to a different value as `logit` is only ",
        "valid for strictly positive {fraction_type(pif)}s."
      )
    )
  }

  return(pif)


}
