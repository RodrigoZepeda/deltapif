#' Potential Impact Fraction related classes
#'
#' Objects for handling potential impact fractions for a categorical exposure
#' considering an observed prevalence of  `p` and a relative risk
#' (or relative risk parameter) of `beta`.
#'
#' @param conf_level Confidence level for the confidence interval (default 0.95).
#'
#' @param link Link function such that the `pif` confidence intervals
#' stays within the expected bounds.
#'
#' @param link_inv The inverse of `link`. For example if `link`
#' is `logit` this should be `inv_logit`.
#'
#' @param link_deriv The derivative of `link`. For example if `link`
#' is `logit` this should be `deriv_logit` (i.e. `function(pif) 1 / (pif * (1 - pif))`).
#'
#' @param pif Potential Impact Fraction estimate
#'
#' @param type Character either Potential Impact Fraction (`PIF`) or
#' Population Attributable Fraction (`PAF`)
#'
#'
#' @param variance variance estimate for the potential impact fraction (i.e.
#' for `pif`)
#'
#' @param p Prevalence (proportion) of the exposed individuals for
#' each of the `N` exposure levels.
#'
#' @param p_cft Counterfactual prevalence (proportion) of the exposed
#' individuals for each of the `N` exposure levels.
#'
#' @param beta Relative risk parameter for which standard deviation is
#' available (usually its either the relative risk directly or the log
#' of the relative risk as most RRs, ORs and HRs come from exponential
#' models).
#'
#' @param var_p Estimate of the colink_variance matrix of `p` where the entry
#' `var_p[i,j]` represents the colink_variance between `p[i]` and `p[j]`.
#'
#' @param var_beta Estimate of the colink_variance matrix of `beta` where the entry
#' `var_beta[i,j]` represents the colink_variance between `beta[i]` and `beta[j]`.
#'
#' @param rr_link Link function such that the relative risk is given by
#' `rr_link(beta)`.
#'
#' @param rr_link_deriv Derivative of the link function for the relative risk.
#' The constructor tries to build it automatically from `rr_link` using
#' [Deriv::Deriv()].
#'
#' @param upper_bound_p Whether the values for the `p` component
#' of the link_variance should be approximated by an upper bound.
#'
#' @param upper_bound_beta Whether the values for the `beta` component
#' of the link_variance should be approximated by an upper bound.
#'
#' @param weights Weights for calculating the total PIF (respectively PAF)
#' in `pif_total`.
#'
#' @param sigma_weights Colink_variance matrix for the weights when calculating the
#' total PIF (respectively PAF) in `pif_total`.
#'
#' @param pif_list A list of potential impact fractions `pif_class` so that
#' the total can be computed from it.
#'
#' @section Properties of a  `pif_class`:
#' \describe{
#'   \item{`ci`}{`numeric(2)` — Lower and upper confidence limits at level `conf_level`.}
#'   \item{`link_vals`}{`numeric` — Entrywise evaluation of the link function at pif: `link(pif)`.}
#'   \item{`link_deriv_vals`}{`character` — Entrywise evaluation of the derivative of the link function (`link_deriv`) at pif: `link(pif)`.}
#'   \item{`link_variance`}{`numeric` - Estimate for the linked potential impact fraction's variance: `variance(link(pif))`.}
#' }
#'
#' @section Properties of a  `pif_atomic_class`:
#' The `pif_atomic_class` inherits the properties of a `pif_class` as well as:
#' \describe{
#'   \item{`mu_obs`}{`numeric` — Average relative risk in the observed population.}
#'   \item{`mu_cft`}{`numeric` — Average relative risk in the counterfactual population.}
#'   \item{`pif`}{`numeric` — Estimate of the potential impact fraction.}
#'   \item{`rr_link_deriv_vals`}{`character` — Entrywise evaluation of the derivative of the link function (`link_deriv`) at pif: `link(pif)`.}
#' }
#'
#' @section Computation of confidence intervals:
#' Wald-type confidence intervals are calculated for `link(pif)` as follows:
#' \deqn{
#'  \text{CI}_{\text{Link}} = \text{link}\big(\text{PIF}\big) \pm Z_{\alpha/2}\cdot\sqrt{\textrm{link\_variance}}
#' }
#' and then transformed back using the inverse of the link function `inv_link`:
#' \deqn{
#'  \text{CI}_{\text{PIF}} = \text{link}^{-1}\Big(\text{CI}_{\text{Link}}\Big)
#' }
#'
#' @section Class structure:
#' All potential impact fractions inherit from the `pif_class` which
#' provides some of the generics.
#'
#' The `pif_atomic_class` is a type of `pif_class` that contains enough
#' information for it to be computed through the classic formula by
#' Walter:
#' \deqn{
#'  \textrm{PIF} = \dfrac{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i - \sum\limits_{i=1}^N p_i^{\text{cft}} \text{RR}_i
#'   }{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i
#'   }
#' }
#' where the relative risk is a function of a parameter \eqn{\beta_i}
#' \deqn{
#'  \text{RR}_i = g(\beta_i)
#' }
#' and the link_variance is calculated for a function of PIF: \eqn{f(\textrm{PIF})}
#' The `pif_atomic_class` only contains one potential impact fraction
#' and the parameters to estimate it.
#'
#' The `pif_additive_class` is a type of `pif_class` that can be computed
#' as a sum of weighted transformations of potential impact fractions
#' where the weights can be random. Elements of a `pif_additive_class`
#' have a link_variance estimated for the following expression:
#' \deqn{
#'  f\big(\textrm{PIF}_{+}\big) = \sum\limits_{i = 1}^{N} q_i \cdot  f_i\Big(\textrm{PIF}_i\Big)
#' }
#'
#' Examples of calculations that can be added to a `pif_additive_class` are:Ç
#'
#' The total potential impact fraction (combining different subpopulations)
#' \deqn{
#'  \textrm{PIF}_{Total} = \sum\limits_{i = 1}^{N} q_i \cdot \textrm{PIF}_i
#' }
#' with \eqn{q_i} representing the proportions of individuals in each category.
#'
#' The ensemble potential impact fraction (representing different relative risks)
#' for the same outcome
#' \deqn{
#' \textrm{PIF}_{Ensemble} = 1 - \prod\limits_{i = 1}^{N} \Big(1 - \textrm{PIF}_i\Big)
#' }
#' as it can be transformed into
#' \deqn{
#' \ln\Big(1 - \textrm{PIF}_{Ensemble}\Big) =  \sum\limits_{i = 1}^{N} \ln\big(1 - \textrm{PIF}_i\big)
#' }
#'
#' @examples
#' #Create a new pif parent class element
#' pif_class(pif = 0.3, variance = 0.01, conf_level = 0.95, type = "PIF",
#'   link = logit, link_inv = inv_logit, link_deriv = deriv_logit)
#'
#' #Create a new potential impact fraction from the Walter's formula
#' pif_atomic_class(
#'   p = 0.499, p_cft = 0, beta = 3.6, var_p = 0.1, var_beta = 3,
#'   link = logit, link_inv = inv_logit, link_deriv = deriv_logit,
#'   rr_link = identity, rr_link_deriv = function(x) 1,
#'   conf_level = 0.95, type = "PAF",
#'   upper_bound_p = FALSE,
#'   upper_bound_beta = FALSE
#' )
#'
#' #Create a list of pif
#' pif1 <- pif_atomic_class(
#'   p = 0.499, p_cft = 0, beta = 3.6, var_p = 0.01, var_beta = 0.03,
#'   link = logit, link_inv = inv_logit, link_deriv = deriv_logit,
#'   rr_link = identity, rr_link_deriv = function(x) 1,
#'   conf_level = 0.95, type = "PAF",
#'   upper_bound_p = FALSE,
#'   upper_bound_beta = FALSE
#' )
#' pif2 <- pif_atomic_class(
#'   p = 0.79, p_cft = 0, beta = 3.6, var_p = 0.01, var_beta = 0.03,
#'   link = logit, link_inv = inv_logit, link_deriv = deriv_logit,
#'   rr_link = identity, rr_link_deriv = function(x) 1,
#'   conf_level = 0.95, type = "PAF",
#'   upper_bound_p = FALSE,
#'   upper_bound_beta = FALSE
#' )
#' pif3 <- pif_atomic_class(
#'   p = 0.8, p_cft = 0, beta = 3.6, var_p = 0.01, var_beta = 0.03,
#'   link = logit, link_inv = inv_logit, link_deriv = deriv_logit,
#'   rr_link = identity, rr_link_deriv = function(x) 1,
#'   conf_level = 0.95, type = "PAF",
#'   upper_bound_p = FALSE,
#'   upper_bound_beta = FALSE
#' )
#'
#' tp1 <- pif_total_class(pif_list = list(pif1, pif2),
#'   weights = c(0.5, 0.2), sigma_weights = diag(0.001, ncol = 2, nrow = 2),
#'   link = identity, link_inv = identity, link_deriv = identity)
#'
#' pif_total_class(pif_list = list(tp1, pif3),
#'   weights = c(0.7, 0.3), sigma_weights = diag(0.001, ncol = 2, nrow = 2),
#'   link = identity, link_inv = identity, link_deriv = identity)
#' @name classes
NULL

#' @rdname classes
#' @export
#pif_class-----
pif_class <- S7::new_class("pif_class",
   package = "pifes",
   properties = list(
     # > Data inputs-----

     # Potential Impact fraction (number)
     pif           = S7::class_numeric,

     # Variance for the potential impact fraction
     variance      = S7::class_numeric,

     # Confidence level
     conf_level    = S7::new_property(S7::class_numeric, setter = set_conf_level),

     # Type
     type          = S7::new_property(S7::class_character, default = "PIF"),

     # PIF link function
     link          = S7::new_property(S7::class_function),

     # PIF link function's inverse
     link_inv      = S7::new_property(S7::class_function),

     # PIF link function's derivative
     link_deriv    = S7::new_property(S7::class_function),

     # > Dynamic properties----

     # Get the derivative of pif link at pif
     link_deriv_vals    = S7::new_property(S7::class_numeric, getter = get_link_deriv_vals),

     # Transform of the pif via the link function
     link_vals     = S7::new_property(S7::class_numeric, getter = get_link_vals),

     # Variance for link(pif)
     link_variance = S7::new_property(S7::class_numeric, getter = get_link_variance),

     # Get confidence intervals
     ci            = S7::new_property(S7::class_numeric, getter = get_ci)
   ),
   validator = function(self){

     if (!is.null(self@conf_level) && (self@conf_level > 1 | self@conf_level < 0)) {
       cli::cli_abort(
         "Invalid confidence level set to {conf_level}"
       )
     }

     if (length(self@pif) < 1){
       cli::cli_abort(
         "No entries provided for `pif`."
       )
     }

     if (length(self@variance) < 1){
       cli::cli_abort(
         "No entries provided for `variance`."
       )
     }

     if (self@variance < 0){
       cli::cli_abort(
         "Invalid (negative) variance provided"
       )
     }

     if (self@pif > 1) {
       cli::cli_abort(
         paste0(
           "Value for the potential impact fraction PIF > 1. ",
          "This is biologically impossible."
         )
       )
     }

     if (length(self@type) != 1 || (self@type != "PAF" && self@type != "PIF")){
       cli::cli_abort(
         "Invalid `type` should be either `PIF` or `PAF`."
       )
     }
  }
)
#S7::S4_register(pif_class)

#' @rdname classes
#' @export
#pif_atomic_class-------
pif_atomic_class <- S7::new_class("pif_atomic_class",
   package = "pifes",
   parent = pif_class,
   properties = list(
     # > Data inputs-----

     # Proportion of individuals exposed (vector or number)
     p             = S7::class_numeric,

     # Proportion of individuals exposed under counterfactual (vector or number)
     p_cft         = S7::class_numeric,

     # Relative risk parameter
     beta          = S7::class_numeric,

     # Colink_variance matrix for p
     var_p       = S7::class_numeric,

     # Colink_variance matrix for beta
     var_beta    = S7::class_numeric,

     # Relative risk link function RR = rr_link(beta)
     rr_link       = S7::class_function,

     # Relative risk link function's derivative
     rr_link_deriv = S7::class_function,

     # > Flags----
     upper_bound_p    = S7::class_logical,

     upper_bound_beta = S7::class_logical,

     # > Calculated properties----

     #Relative risk values
     rr     = S7::new_property(S7::class_numeric, getter = get_rr),

     # Mean of the relative risk under the observed prevalence
     mu_obs = S7::new_property(S7::class_numeric, getter = get_mu_obs),

     # Mean of the relative risk under the counterfactual scenario
     mu_cft = S7::new_property(S7::class_numeric, getter = get_mu_cft),

     #Potential impact fraction
     pif    = S7::new_property(S7::class_numeric, getter = get_pif),

     # Get the derivative of rr link at pif
     rr_link_deriv_vals = S7::new_property(S7::class_numeric, getter = get_rr_link_deriv_vals),

     # link_variance estimate (with linking function)
     variance           = S7::new_property(S7::class_numeric, get_variance_atomic)


   ),
   validator = function(self) {
     # Validate the length of p vs p_cft
     if (length(self@p) != length(self@p_cft)) {
       cli::cli_abort(
         paste0(
           "Exposure prevelence vector {.code p} and ",
           "counterfactual exposure prevalence vector ",
           "{.code p_cft} must be of the same length."
         )
       )
     }

     # Check the exposures are positive and sum < 1
     if (any(self@p < 0)) {
       cli::cli_abort(
         paste0(
            "Exposure prevelence vector {.code p} ",
            "has values < 0."
         )
       )
     }

     if (any(self@p_cft < 0)) {
       cli::cli_abort(
         paste0(
           "Counterfactual exposure prevelence vector ",
           "{.code p_cft} has values < 0."
         )
       )
     }

     if (sum(self@p) > 1) {
       cli::cli_abort(
         paste0(
            "Exposure prevelence vector {.code p} ",
            "sums to values > 1 ",
            "{.code sum(p) = {sum(self@p)}}"
         )
       )
     }

     if (sum(self@p_cft) > 1) {
       cli::cli_abort(
         paste0(
           "Counterfactual exposure prevelence vector ",
           "{.code p_cft} adds to values > 1 ",
           "{.code sum(p_cft) = {sum(self@p_cft)}}"
         )
       )
     }

     # Validate the length of p vs relative risks (rr_link(beta))
     if (length(self@p) != length(self@rr_link(self@beta))) {
       cli::cli_abort(
         paste0(
            "Exposure prevelence vector {.code p} and ",
            "relative risks {.code rr_link(beta)} have ",
            "different lengths."
         )
       )
     }

     # Check the matrices are positive definite
     if (is.matrix(self@var_beta) && !isSymmetric(self@var_beta, trans = "T")) {
       cli::cli_abort(
         paste0(
            "Matrix {.code var_beta} is not symmetric. ",
            "Entry in row i and column j must be equal ",
            "to entry in row j and column i."
         )
       )
     }

     # Check the matrices are positive definite
     if (is.matrix(self@var_p) && !isSymmetric(self@var_p, trans = "T")) {
       cli::cli_abort(
         paste0(
           "Matrix {.code var_p} is not symmetric. ",
           "Entry in row i and column j must be equal ",
           "to entry in row j and column i."
         )
       )
     }

     if (is.matrix(self@var_p) &&
         (
           (ncol(self@var_p) != length(self@p)) |
           (nrow(self@var_p) != length(self@p))
          )
         ) {
         cli::cli_abort(
           paste0(
             "Exposure prevalence vector {.code p} has ",
             "different length than its colink_variance matrix",
             "{.code var_p}."
           )
         )
     }

     # Validate the length of var_beta
     if (is.matrix(self@var_beta) &&
         (
           (ncol(self@var_beta) != length(self@beta)) |
           (nrow(self@var_beta) != length(self@beta))
         )
        ){
         cli::cli_abort(
           paste0(
             "Exposure prevalence vector {.code p} has ",
             "different length than its colink_variance ",
             "matrix {.code var_p}."
           )
         )
     }

     # Validate the calculations
     if (self@mu_obs <= 0) {
       cli::cli_abort(
         paste0(
           "The average relative risk in the population ",
           "`mu_obs` is <= 0. This is biologically ",
           "impossible."
         )
       )
     }
     if (self@mu_cft < 0) {
       cli::cli_abort(
         paste0(
           "The average relative risk in the ",
           "counterfactual population `mu_cft` is <= 0. ",
           "This is biologically impossible."
         )
       )
     }
   },
   constructor = function(p, p_cft, beta, var_p, var_beta, rr_link,
                          rr_link_deriv, link, link_deriv, link_inv, conf_level,
                          type, upper_bound_p, upper_bound_beta) {


     S7::new_object(S7::S7_object(),

                    #Elements of pif_class
                    conf_level = conf_level,
                    type       = type,
                    link       = link,
                    link_inv   = link_inv,
                    link_deriv = link_deriv,

                    #New added elements of pif_atomic
                    p             = p,
                    p_cft         = p_cft,
                    beta          = beta,
                    var_p       = var_p,
                    var_beta    = var_beta,
                    rr_link       = rr_link,
                    rr_link_deriv = rr_link_deriv,
                    upper_bound_p = upper_bound_p,
                    upper_bound_beta = upper_bound_beta


     )
  }
)
#S7::S4_register(pif_atomic_class)

#' @rdname classes
#' @export
#pif_total_class------------------
pif_total_class <- S7::new_class(
  name      = "pif_total_class",
  package   = "pifes",
  parent    = pif_class,
  properties = list(
    pif_list      = S7::class_list,
    weights       = S7::class_numeric,
    sigma_weights = S7::class_numeric,
    type          = S7::new_property(S7::class_numeric, getter = get_total_type),
    coefs         = S7::new_property(S7::class_numeric, getter = get_total_coefs),
    pif           = S7::new_property(S7::class_numeric, getter = get_total_pif),
    covariance    = S7::new_property(S7::class_numeric, getter = get_covariance_total),
    variance      = S7::new_property(S7::class_numeric, getter = get_variance_total)
  ),
  validator = function(self) {

    for (i in seq_along(self@pif_list)) {
      if (!S7::S7_inherits(self@pif_list[[i]], pif_class)) {
        cli::cli_abort(
          "Element {i} of `pif_list` must be a 'pif_class'."
        )
      }
    }

    if (length(self@weights) != length(self@pif_list)){
      cli::cli_abort(
        paste0(
          "Weights provided have length {length(self@weights)} but",
          "{length(self@pif_list)} fractions were provided."
        )
      )
    }
  },
  constructor = function(pif_list, weights, sigma_weights,
                         conf_level = 0.95, link, link_inv, link_deriv){

    S7::new_object(S7::S7_object(),
                   conf_level = conf_level,
                   link = link,
                   link_inv = link_inv,
                   link_deriv = link_deriv,
                   pif_list = pif_list,
                   weights = weights,
                   sigma_weights = sigma_weights
                   )


  }
)

#' @rdname classes
#' @export
#pif_ensemble_class------------------
pif_ensemble_class <- S7::new_class(
  name      = "pif_ensemble_class",
  package   = "pifes",
  parent    = pif_total_class,
  properties = list(
    pif_list      = S7::class_list,
    type          = S7::new_property(S7::class_numeric, getter = get_ensemble_type),
    coefs         = S7::new_property(S7::class_numeric, getter = get_ensemble_coefs),
    pif           = S7::new_property(S7::class_numeric, getter = get_ensemble_pif),
    covariance    = S7::new_property(S7::class_numeric, getter = get_ensemble_covariance),
    variance      = S7::new_property(S7::class_numeric, getter = get_ensemble_variance)
  ),
  validator = function(self) {

    for (i in seq_along(self@pif_list)) {
      if (!S7::S7_inherits(self@pif_list[[i]], pif_class)) {
        cli::cli_abort(
          "Element {i} of `pif_list` must be a 'pif_class'."
        )
      }
    }
  },
  constructor = function(pif_list, conf_level = 0.95){

    S7::new_object(S7::S7_object(),
                   conf_level = conf_level,
                   link = log_complement,
                   link_inv = inv_log_complement,
                   link_deriv = deriv_log_complement,
                   pif_list = pif_list,
                   weights = rep(1, length(pif_list)),
                   sigma_weights = matrix(0, ncol = length(pif_list), nrow = length(pif_list))
    )


  }
)
