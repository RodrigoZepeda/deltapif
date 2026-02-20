#' Potential Impact Fraction related classes
#'
#' Objects for handling potential impact fractions for a categorical exposure
#' considering an observed prevalence of  `p` and a relative risk
#' (or relative risk parameter) of `beta`.
#'
#' @param pif_obj A `pif_class` to calculate the percentage of cases with.
#'
#' @param overall_cases Overall number of cases so that attributable cases is given by
#' `overall_cases * pif_obj`
#'
#' @param variance_cases Variance of the cases estimate.
#'
#' @param conf_level Confidence level for the confidence interval (default 0.95).
#'
#' @param link Link function such that the `pif` confidence intervals
#' stays within the expected bounds.
#'
#' @param link_inv (Optional). If `link` is a function then
#' yhe inverse of `link`. For example if `link` is `logit` this
#' should be `inv_logit`.
#'
#' @param link_deriv (Optional) If `link` is a function, the derivative of `link`. For example if `link`
#' is `logit` this should be `deriv_logit` (i.e. `function(pif) 1 / (pif * (1 - pif))`).
#'
#' @param pif Potential Impact Fraction point estimate
#'
#' @param type Character either Potential Impact Fraction (`PIF`) or
#' Population Attributable Fraction (`PAF`)
#'
#' @param label Character identifier for the impact fraction. This is for
#  ease of the user.
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
#' @param var_p Estimate of the link_covariance matrix of `p` where the entry
#' `var_p[i,j]` represents the link_covariance between `p[i]` and `p[j]`.
#'
#' @param var_beta Estimate of the link_covariance matrix of `beta` where the entry
#' `var_beta[i,j]` represents the link_covariance between `beta[i]` and `beta[j]`.
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
#' @param weights weights for calculating the total PIF (respectively PAF)
#' in `pif_total`.
#'
#' @param var_weights covariance matrix for the weights when calculating the
#' total PIF (respectively PAF) in `pif_total`.
#'
#' @param var_pif_weights covariance matrix with row `i` and column `j`
#' representing the covariance between the `i`-th potential impact
#' fraction of the list and the `j`-th weight
#'
#' @param pif_list A list of potential impact fractions `pif_class` so that
#' the total can be computed from it.
#'
#' @param pif_transform Transform applied to the `pif` for summation
#' in a `pif_global_ensemble_class` (see section below).
#'
#' @param pif_deriv_transform Derivative of the transform applied to the
#' `pif` for summation in a `pif_global_ensemble_class` (see section below).
#'
#' @param pif_inverse_transform Inverse of the transform applied to the
#' `pif` for summation in a `pif_global_ensemble_class` (see section below).
#'
#'
#' @section Properties of a  `pif_class`:
#' Any object that is a `pif_class` contains a potential impact fraction
#' with intervals estimated as follows:
#' \deqn{
#'  \text{CI}_{\text{Link}} = \text{link}\big(\text{PIF}\big) \pm Z_{\alpha/2}\cdot\sqrt{\textrm{link\_variance}}
#' }
#' and then transformed back using the inverse of the link function `inv_link`:
#' \deqn{
#'  \text{CI}_{\text{PIF}} = \text{link}^{-1}\Big(\text{CI}_{\text{Link}}\Big)
#' }
#'
#' The following are the properties of any `pif_class`
#' \describe{
#'   \item{`ci`}{`numeric(2)` — Lower and upper confidence limits at level `conf_level`.}
#'   \item{`link_vals`}{`numeric` — Entrywise evaluation of the link function at pif: `link(pif)`.}
#'   \item{`link_deriv_vals`}{`character` — Entrywise evaluation of the derivative of the link function (`link_deriv`) at pif: `link(pif)`.}
#'   \item{`link_variance`}{`numeric` - Estimate for the linked potential impact fraction's variance: `variance(link(pif))`.}
#' }
#'
#' @section Properties of a  `pif_atomic_class`:
#'
#' The `pif_atomic_class` is a type of `pif_class` that contains enough
#' information to compute a potential impact fraction
#' through the classic formula by Walter:
#' \deqn{
#'  \textrm{PIF} = \dfrac{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i - \sum\limits_{i=1}^N p_i^{\text{cft}} \text{RR}_i
#'   }{
#'    \sum\limits_{i=1}^N p_i \text{RR}_i
#'   }
#' }
#' where the relative risk is a function of a parameter \eqn{\beta_i}
#' \deqn{
#'  \text{RR}_i = \text{rr\_link}(\beta_i)
#' }
#'
#' The `pif_atomic_class` inherits the properties of a `pif_class` as well as:
#' \describe{
#'   \item{`mu_obs`}{`numeric` — Average relative risk in the observed population.}
#'   \item{`mu_cft`}{`numeric` — Average relative risk in the counterfactual population.}
#'   \item{`pif`}{`numeric` — Estimate of the potential impact fraction.}
#'   \item{`rr_link_deriv_vals`}{`character` — Entrywise evaluation of the derivative of the link function (`link_deriv`) at pif: `link(pif)`.}
#' }
#'
#' Confidence intervals are estimated as with any `pif_class`.
#'
#' @section Properties of a  `pif_global_ensemble_class`:
#'
#' The `pif_global_ensemble_class` creates a new potential impact
#' fraction by summing a weighted combination of potential
#' impact fractions. In general it computes the following expression:
#' \deqn{
#' \textrm{PIF}_{\text{global}} =
#'  g^{-1}\bigg( \sum\limits_{i = 1}^{N} g\big(w_i \cdot \textrm{PIF}_i\big) \bigg)
#' }
#' where \eqn{g} is refered to as the `pif_transform`, its derivative the
#' `pif_deriv_transform`, and its inverse `pif_inverse_transform`.
#'
#' The `pif_global_ensemble_class` inherits the properties of a `pif_class` as well as:
#' \describe{
#'  \item{`weights`}{`numeric` - Vector of weights \eqn{w_i} for weighting the potential impact fraction.}
#'  \item{`var_weights`}{`numeric` - Covariance matrix for the `weights`}
#'  \item{`pif_transform`}{`function` - Function \eqn{g} with which to transform the impact fraction before weighting.}
#'  \item{`pif_deriv_transform`}{`function` - Derivative of the `pif_transform`.}
#'  \item{`pif_inverse_transform`}{`function` - Inverse of the `pif_transform`.}
#'  \item{`type`}{`character` - Whether the quantity represents a `PIF` or a `PAF`}
#'  \item{`coefs`}{`numeric` - Potential impact fractions used for the global ensemble (each of the \eqn{\text{PIF}_i}.}
#'  \item{`sum_transformed_weighted_coefs`}{`numeric` - Sum of the potential impact fractions involved \eqn{\sum g(w_i \text{PIF}_i)}.}
#'  \item{`pif`}{`numeric` — Estimate of the potential impact fraction.}
#'  \item{`covariance`}{`numeric` — Covariance matrix between the potential impact fractions in `coefs` (i.e. each entry is\eqn{\text{Cov}(\text{PIF}_i, \text{PIF}_j))}}
#'  \item{`variance`}{`numeric` — Estimate for the variance of `pif`.}
#' }
#'
#' Confidence intervals are estimated as with any `pif_class`.
#'
#' @section Properties of a  `pif_total_class`:
#'
#' A `pif_total_class` estimated the potential impact fraction of the
#' weighted sum of fractions from different (disjoint) dispopulations:
#' \deqn{
#'  \textrm{PIF}_{Total} = \sum\limits_{i = 1}^{N} w_i \cdot \textrm{PIF}_i
#' }
#' with \eqn{w_i} representing the proportions of individuals in each category.
#' This is a type of `pif_global_ensemble_class` with `pif_transform = identity`.
#'
#' @section Properties of a  `pif_ensemble_class`:
#'
#' The ensemble potential impact fraction (representing different relative risks)
#' for the same outcome is given by the weighted product:
#' \deqn{
#' \textrm{PIF}_{Ensemble} = 1 - \prod\limits_{i = 1}^{N} \Big(1 - w_i \textrm{PIF}_i\Big)
#' }
#'
#' However it can be transformed into a `pif_global_ensemble_class` by
#' taking the log-complement:
#' \deqn{
#' \ln\Big(1 - \textrm{PIF}_{Ensemble}\Big) =  \sum\limits_{i = 1}^{N} \ln\big(1 - w_i \textrm{PIF}_i\big)
#' }
#' hence it is a  `pif_global_ensemble_class` with `pif_transform = log_complement`.
#'
#' @name classes
NULL




#' @rdname classes
#pif_class-----
pif_class <- S7::new_class("pif_class",
   package = "deltapif",
   properties = list(
     # > Data inputs-----

     # Potential Impact fraction (number)
     pif           = S7::class_numeric,

     # Variance for the potential impact fraction
     variance      = S7::class_numeric,

     # Confidence level
     conf_level    = S7::new_property(S7::class_numeric, setter = set_conf_level),

     #Label for easy identification
     label         = S7::class_character,

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
       cli::cli_warn(
         paste0(
           "Value for the potential impact fraction PIF > 1. ",
          "This is biologically impossible. Are you calculating correctly? (in some cases this might make statistical sense)"
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
#cases_class-----
cases_class <- S7::new_class("cases_class",
   package = "deltapif",
   properties = list(
     # > Data inputs-----
     overall_cases = S7::class_numeric,

     # Fraction
     pif_obj       = S7::class_any,

     # Variance for the potential impact fraction
     variance_cases = S7::class_numeric,


     # Confidence level
     conf_level    = S7::new_property(S7::class_numeric, setter = set_conf_level),

     # PIF link function
     link          = S7::new_property(S7::class_function),

     # PIF link function's inverse
     link_inv      = S7::new_property(S7::class_function),

     # PIF link function's derivative
     link_deriv    = S7::new_property(S7::class_function),

     # > Dynamic properties----
     cases         = S7::new_property(S7::class_numeric, getter = get_cases),

     # Variance
     variance      = S7::new_property(S7::class_numeric, getter = get_variance_cases),

     # Get the derivative of pif link at pif
     link_deriv_vals = S7::new_property(S7::class_numeric, getter = get_link_deriv_vals_cases),

     # Transform of the pif via the link function
     link_vals     = S7::new_property(S7::class_numeric, getter = get_link_vals_cases),

     # Variance for link(pif)
     link_variance = S7::new_property(S7::class_numeric, getter = get_link_variance),

     # Get confidence intervals
     ci            = S7::new_property(S7::class_numeric, getter = get_ci)

   ),
   constructor = function(overall_cases, pif_obj, variance_cases, link, link_deriv, link_inv, conf_level) {

     S7::new_object(S7::S7_object(),
                    overall_cases = overall_cases,
                    pif_obj    = pif_obj,
                    variance_cases   = variance_cases,
                    conf_level = conf_level,
                    link       = link,
                    link_inv   = link_inv,
                    link_deriv = link_deriv)
   }
)

#' @rdname classes
#pif_atomic_class-------
pif_atomic_class <- S7::new_class("pif_atomic_class",
   package = "deltapif",
   parent = pif_class,
   properties = list(
     # > Data inputs-----

     # Proportion of individuals exposed (vector or number)
     p             = S7::class_numeric,

     # Proportion of individuals exposed under counterfactual (vector or number)
     p_cft         = S7::class_numeric,

     # Relative risk parameter
     beta          = S7::class_numeric,

     # link_covariance matrix for p
     var_p       = S7::class_numeric,

     # link_covariance matrix for beta
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
       cli::cli_alert_warning(
         paste0(
            "Matrix {.code var_beta} is not symmetric. ",
            "Entry in row i and column j must be equal ",
            "to entry in row j and column i."
         )
       )
     }

     # Check the matrices are positive definite
     if (is.matrix(self@var_p) && !isSymmetric(self@var_p, trans = "T")) {
       cli::cli_alert_warning(
         paste0(
           "Matrix {.code var_p} is not symmetric. ",
           "Entry in row i and column j must be equal ",
           "to entry in row j and column i."
         )
       )
     }

     if (is.matrix(self@var_p) &&
         (
           (ncol(self@var_p) != length(self@p)) ||
           (nrow(self@var_p) != length(self@p))
          )
         ) {
         cli::cli_abort(
           paste0(
             "Exposure prevalence vector {.code p} has ",
             "different length than its covariance matrix",
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
             "Relative risk parameter vector {.code beta} has ",
             "different length than its covariance ",
             "matrix {.code var_beta}."
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
                          type, label, upper_bound_p, upper_bound_beta) {


     S7::new_object(S7::S7_object(),

                    #Elements of pif_class
                    conf_level = conf_level,
                    type       = type,
                    label      = label,
                    link       = link,
                    link_inv   = link_inv,
                    link_deriv = link_deriv,

                    #New added elements of pif_atomic
                    p             = p,
                    p_cft         = p_cft,
                    beta          = beta,
                    var_p         = var_p,
                    var_beta      = var_beta,
                    rr_link       = rr_link,
                    rr_link_deriv = rr_link_deriv,
                    upper_bound_p = upper_bound_p,
                    upper_bound_beta = upper_bound_beta


     )
  }
)


#' @rdname classes
#pif_global_ensemble_class-------
pif_global_ensemble_class <- S7::new_class(
  name = "pif_global_ensemble_class",
  package = "deltapif",
  parent = pif_class,
  properties = list(
    pif_list                = S7::class_any,
    weights                 = S7::class_numeric,
    var_weights             = S7::class_any,
    var_pif_weights         = S7::class_any,
    pif_transform           = S7::class_function,
    pif_deriv_transform     = S7::class_function,
    pif_inverse_transform   = S7::class_function,
    type                    = S7::new_property(S7::class_character, getter = get_ensemble_type),
    coefs                   = S7::new_property(S7::class_numeric, getter = get_ensemble_coefs),
    sum_transformed_weighted_coefs = S7::new_property(S7::class_numeric, getter = get_sum_transformed_weighted_coefs),
    pif                            = S7::new_property(S7::class_numeric, getter = get_global_ensemble_pif),
    covariance            = S7::new_property(S7::class_numeric, getter = get_covariance),
    variance              = S7::new_property(S7::class_numeric, getter = get_variance)
  ),
  validator = validate_global_ensemble,
  constructor = function(pif_list, weights, var_weights, var_pif_weights,
                         pif_transform,
                         pif_deriv_transform, pif_inverse_transform,
                         link, link_inv, link_deriv, conf_level = 0.95,  label){

    S7::new_object(S7::S7_object(),
                   conf_level = conf_level,
                   pif_transform = pif_transform,
                   label      = label,
                   pif_inverse_transform = pif_inverse_transform,
                   pif_deriv_transform = pif_deriv_transform,
                   link = link,
                   link_inv = link_inv,
                   link_deriv = link_deriv,
                   pif_list = pif_list,
                   weights = weights,
                   var_weights = var_weights,
                   var_pif_weights = var_pif_weights,
                   variance = 0
    )}
)

#' @rdname classes
#pif_total_class------------------
pif_total_class <- S7::new_class(
  name      = "pif_total_class",
  package   = "deltapif",
  parent    = pif_global_ensemble_class,
  validator = validate_global_ensemble,
  constructor = function(pif_list, weights, var_weights, var_pif_weights,
                         link, link_inv, link_deriv,
                         conf_level = 0.95, label){

    S7::new_object(S7::S7_object(),
                   conf_level = conf_level,
                   label      = label,
                   link = link,
                   link_inv = link_inv,
                   link_deriv = link_deriv,
                   pif_list = pif_list,
                   weights = weights,
                   var_weights = var_weights,
                   var_pif_weights = var_pif_weights,
                   pif_transform = identity,
                   pif_deriv_transform = function(x) rep(1, length(x)),
                   pif_inverse_transform = identity,
                   variance = 0
                   )


  }
)

#' @rdname classes
#pif_ensemble_class------------------
pif_ensemble_class <- S7::new_class(
  name      = "pif_ensemble_class",
  package   = "deltapif",
  parent    = pif_global_ensemble_class,
  validator = validate_global_ensemble,
  constructor = function(pif_list,
                         weights,
                         var_weights,
                         var_pif_weights,
                         link,
                         link_inv,
                         link_deriv,
                         conf_level = 0.95,
                         label){

    S7::new_object(S7::S7_object(),
                   conf_level = conf_level,
                   label = label,
                   link = link,
                   link_inv = link_inv,
                   link_deriv = link_deriv,
                   weights = weights,
                   var_weights = var_weights,
                   var_pif_weights = var_pif_weights,
                   pif_list = pif_list,
                   pif_transform = log_complement,
                   pif_inverse_transform = inv_log_complement,
                   pif_deriv_transform = deriv_log_complement,
                   variance = 0
    )


  }
)
