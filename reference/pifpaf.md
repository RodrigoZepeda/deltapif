# Potential Impact fraction and Population Attributable Fraction

Calculates the potential impact fraction `pif` or the population
attributable fraction `paf` for a categorical exposure considering an
observed prevalence of `p` and a relative risk (or relative risk
parameter) of `beta`.

## Usage

``` r
paf(
  p,
  beta,
  var_p = NULL,
  var_beta = NULL,
  rr_link = "exponential",
  rr_link_deriv = NULL,
  link = "log-complement",
  link_inv = NULL,
  link_deriv = NULL,
  conf_level = 0.95,
  quiet = FALSE,
  label = NULL
)

pif(
  p,
  p_cft = rep(0, length(p)),
  beta,
  var_p = NULL,
  var_beta = NULL,
  rr_link = "exponential",
  rr_link_deriv = NULL,
  link = "log-complement",
  link_inv = NULL,
  link_deriv = NULL,
  conf_level = 0.95,
  type = "PIF",
  quiet = FALSE,
  label = NULL
)
```

## Arguments

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- beta:

  Relative risk parameter for which standard deviation is available
  (usually its either the relative risk directly or the log of the
  relative risk as most RRs, ORs and HRs come from exponential models).

- var_p:

  Estimate of the link_covariance matrix of `p` where the entry
  `var_p[i,j]` represents the link_covariance between `p[i]` and `p[j]`.

- var_beta:

  Estimate of the link_covariance matrix of `beta` where the entry
  `var_beta[i,j]` represents the link_covariance between `beta[i]` and
  `beta[j]`.

- rr_link:

  Link function such that the relative risk is given by `rr_link(beta)`.

- rr_link_deriv:

  Derivative of the link function for the relative risk. The function
  tries to build it automatically from `rr_link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- quiet:

  Whether to show messages.

- label:

  Character identifier for the impact fraction. This is for

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

- type:

  Character either Potential Impact Fraction (`PIF`) or Population
  Attributable Fraction (`PAF`)

## Note

This function assumes `p` and `beta` have been pre-computed from the
data and the individual-level data are not accessible to the
researchers. If either the data for the individual-level prevalence of
exposure `p` or the data for the individual-level risk estimate `beta`
can be accessed by the researcher other methods (such as the `pifpaf`
package should be preferred).

## Formulas

This function computed the potential impact fraction and its confidence
intervals using Walter's formula:

\$\$ \textrm{PIF} = \dfrac{ \sum\limits\_{i=1}^N p_i \text{RR}\_i -
\sum\limits\_{i=1}^N p_i^{\text{cft}} \text{RR}\_i }{
\sum\limits\_{i=1}^N p_i \text{RR}\_i }, \quad \text{ and } \quad
\textrm{PAF} = \dfrac{ \sum\limits\_{i=1}^N p_i \text{RR}\_i - 1 }{
\sum\limits\_{i=1}^N p_i \text{RR}\_i } \$\$

in the case of `N` exposure categories which is equivalent to Levine's
formula when there is only `1` exposure category:

\$\$ \textrm{PIF} = \dfrac{ p (\text{RR} - 1) - p^{\text{cft}}
(\text{RR} - 1) }{ 1 + p (\text{RR} - 1) } \quad \textrm{ and } \quad
\textrm{PAF} = \dfrac{ p (\text{RR} - 1) }{ 1 + p (\text{RR} - 1) } \$\$

The construction of confidence intervals is done via a link function.
Its variance is given by:

\$\$ \sigma_f^2 = \text{Var}\Big\[ f(\textrm{PIF}) \Big\] \approx
\Big(f'(\textrm{PIF})\Big)^2 \text{Var}\Big\[ \textrm{PIF} \Big\] \$\$
and the intervals are constructed as: \$\$ \text{CI} = f^{-1}\Big(
f(\textrm{PIF}) \pm Z\_{\alpha/2} \cdot \sigma_f \Big) \$\$

the function \\f\\ is called the link function for the PIF, \\f^{-1}\\
represents its inverse, and \\f'\\ its derivative.

## Link functions for the PIF

By default the `pif` and `paf` calculations use the `log-complement`
link which guarantees the impact fractions' intervals have valid values
(between `-Inf` and `1`).

Depending on the application the following link functions are also
implemented:

- log-complement:

  To achieve fractions between `-Inf` and `1`. This is the function
  `f(x) = ln(1 - x)` with inverse `finv(x) = 1 - exp(x)`

- logit:

  To achieve strictly positive fractions between `0` and `1`. This is
  the function `f(x) = ln(x / (1 - x))` with inverse
  `finv(x) = 1 / (1 + exp(-x))`

- identity:

  An approximation for fractions that does not guarantee confidence
  intervals in valid regions. This is the function `f(x) = x`.with
  inverse `finv(x) = x`

- hawkins:

  Hawkins' fraction for controlling the variance. This is the function
  `f(x) = ln(x + sqrt(x^2 + 1))` with inverse
  `finv(x) = 0.5 * exp(-x) * (exp(2 * x) - 1)`

- custom:

  User specified (see below).

In general, `logit` should be preferred if it is known and certain that
the fractions can only be positive (i.e. when all relative risks
(including CIs) are \> 1 and prevalence \> 0 and there is an
epidemiological / biological justification).

**Custom link** functions can be implemented as long as they are
invertible in the range of interest by providing:

- link:

  The function

- link_inv:

  The inverse function of `link`

- link_deriv:

  The derivative of `link`

If no derivative is provided the package attempts to estimate it
symbolically using
[`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html) however
there is no guarantee that this will work for non-standard functions
(i.e. not logarithm / trigonometric / exponential)

As an example considering implementing a square root custom link:

    # The link is square root
    link      <- sqrt

    # Inverse of sqrt(x) is x^2
    link_inv  <- function(x) x^2

    # The derivative of sqrt(x) is 1/(2 * sqrt(x))
    link_deriv <- function(pif) 1 / (2 * sqrt(pif))

    # Then the pif can be calculated as:
    pif(p = 0.499, beta = 1.6, p_cft = 0.499/2, var_p = 0.1, var_beta = 0.2,
         link = link, link_inv = link_inv, link_deriv = link_deriv)

## Link functions for beta

By default the `pif` and `paf` use the `exponential` link which means
that the values for `beta` are the log-relative risks and the variance
`var_beta` corresponds to the log-relative risk's variance. Depending on
the relative risk's source the following options are available:

- exponential:

  For when the relative risks correspond to the exponential of another
  parameter, usually called `beta`. This is the exponential function
  `f(beta) = exp(beta)` with inverse `finv(rr) = log(rr) = beta`

- identity:

  For when the relative risk and their variance are reported directly.
  This is the function `f(beta) = beta` with inverse
  `finv(rr) = rr = beta`

- custom:

  User specified (see below).

**Note** that in most cases including contingency tables, Poisson and
Logistic regressions, and Cox Proportional Hazards, the relative risks
are estimated by exponentiating a parameter.

As in the previous section, custom link functions for beta can be
implemented as long as they are invertible in the range of interest by
providing the function `rr_link` and its derivative `rr_link_deriv`. If
no derivative is provided the package does an attempt to estimate it
symbolically using
[`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html) however
there is no guarantee that this will work non-standard functions (i.e.
not logarithm / trigonometric / exponential)

## Population Attributable Fraction

The population attributable fraction corresponds to the potential impact
fraction at the theoretical minimum risk level. It is assumed that the
theoretical minimum risk level is a relative risk of 1. If no
counterfactual prevalence `p_cft` is specified, the model computes the
population attributable fraction.

## See also

[`pif_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
[`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
[`paf_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
[`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),

## Examples

``` r
#EXAMPLE 1: ONE EXPOSURE CATEGORY (CLASSIC LEVIN)
#---------------------------------------------------------------------------
# This example comes from Levin 1953
# Relative risk of lung cancer given smoking was 3.6
# Proportion of individuals smoking where 49.9%
# Calculates PAF (i.e. counterfactual is no smoking)
paf(p = 0.499, beta = log(3.6))
#> ! Assuming parameters `p` have no variance Use `var_p` to input their link_variances and/or covariance
#> ! Assuming parameters `beta` have no variance Use `var_beta` to input their link_variances and/or covariance
#> 
#> ── Population Attributable Fraction: [deltapif-112383884326229] ──
#> 
#> PAF = 56.473% [95% CI: 56.473% to 56.473%]
#> standard_deviation(paf %) = 0.000

# Assuming that beta and p had a link_variance
paf(p = 0.499, beta = log(3.6), var_p = 0.001, var_beta = 0.1)
#> 
#> ── Population Attributable Fraction: [deltapif-0397001493059625] ──
#> 
#> PAF = 56.473% [95% CI: 28.972% to 73.326%]
#> standard_deviation(paf %) = 10.875

# If the link_variance was to high a logistic transform would be required
# Generates incorrect values for the interval:
paf(p = 0.499, beta = log(3.6), var_p = 0.1, var_beta = 0.3)
#> 
#> ── Population Attributable Fraction: [deltapif-0823261150567127] ──
#> 
#> PAF = 56.473% [95% CI: -29.969% to 85.422%]
#> standard_deviation(paf %) = 24.294

# Logit fixes it
paf(p = 0.499, beta = log(3.6), var_p = 0.1, var_beta = 0.3,
    link = "logit", quiet = TRUE)
#> 
#> ── Population Attributable Fraction: [deltapif-05788846247121] ──
#> 
#> PAF = 56.473% [95% CI: 15.753% to 90.002%]
#> standard_deviation(paf %) = 24.294

# If the counterfactual was reducing the smoking population by 1/2
pif(p = 0.499, beta = log(3.6), p_cft = 0.499/2, var_p = 0.001,
    var_beta = 0.1, link = "logit", quiet = TRUE)
#> 
#> ── Potential Impact Fraction: [deltapif-17637893775715] ──
#> 
#> PIF = 28.236% [95% CI: 18.101% to 41.192%]
#> standard_deviation(pif %) = 5.963

#EXAMPLE 2: MULTIPLE EXPOSURE CATEGORIES
#---------------------------------------------------------------------------
#In "Excess Deaths Associated With Underweight, Overweight, and Obesity"
#the PAF of mortality associated to different BMI levels is calculated

#These are the relative risks for age-group 25-59 for each BMI level:
rr <- c("<18.5" = 1.38, "25 to <30" = 0.83 ,
        "30 to <35" = 1.20, ">=35" = 1.83)

#While the prevalences are:
p <- c("<18.5" = 1.9, "25 to <30" = 34.8,
        "30 to <35" = 17.3, ">=35" = 13.3) / 100

#The variance of the relative risk is obtained from the CIs:
var_log_rr <- c("<18.5" = 0.2653156, "25 to <30" =  0.1247604,
                "30 to <35" = 0.1828293, ">=35" = 0.1847374)^2

#Note that we are omitting the group "18.5 to < 25" as it is the
#reference (RR = 1)
paf(p = p, beta = rr, var_beta = var_log_rr, var_p = 0, quiet = TRUE,
    link = "logit")
#> 
#> ── Population Attributable Fraction: [deltapif-0132992146404077] ──
#> 
#> PAF = 61.599% [95% CI: 50.274% to 71.792%]
#> standard_deviation(paf %) = 5.571

#We can compute a potential impact fraction of, for example, reducing the
#amount of people over 35 to 0 and having them in the "30 to <35":
p_cft <- c("<18.5" = 1.9, "25 to <30" = 34.8,
            "30 to <35" = 17.3 + 13.3, ">=35" = 0) / 100

#The potential impact fraction is as follows:
pif(p = p, p_cft = p_cft, beta = rr, link = "logit",
  var_beta = var_log_rr, var_p = 0, quiet = TRUE)
#> 
#> ── Potential Impact Fraction: [deltapif-0376499328478288] ──
#> 
#> PIF = 14.882% [95% CI: 13.702% to 16.144%]
#> standard_deviation(pif %) = 0.623
```
