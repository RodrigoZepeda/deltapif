# Components of the variance of a potential impact fraction

Splits the delta-method variance of a potential impact fraction (or
population attributable fraction) into the part due to the uncertainty
in the exposure prevalences `p`, the part due to the uncertainty in the
relative-risk parameters `beta`, and, for totals and ensembles, the part
due to the uncertainty in the `weights`.

## Usage

``` r
variance_components(x, ...)
```

## Arguments

- x:

  A potential impact fraction created with
  [`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md),
  [`paf()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md),
  [`pif_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  or their `paf` counterparts.

- ...:

  Not used.

## Value

A named numeric vector with entries `p`, `beta` and `weights`.

## Details

As the prevalences and the relative risks are estimated independently,
the variance of the fraction is the sum of a prevalence component and a
relative-risk component: \$\$ \text{Var}\big(\widehat{\textrm{PIF}}\big)
\approx \nabla_p \textrm{PIF}^{\top} \Sigma_p \nabla_p \textrm{PIF} +
\nabla\_\beta \textrm{PIF}^{\top} \Sigma\_\beta \nabla\_\beta
\textrm{PIF}. \$\$ Comparing both components shows whether a more
precise exposure survey or a more precise relative risk would most
reduce the uncertainty of the fraction. For totals and ensembles the
remaining variance (from `var_weights` and `var_pif_weights`) is
reported as the `weights` component. Components are on the scale of the
fraction (not of the link) and add up to
[`variance()`](https://rodrigozepeda.github.io/deltapif/reference/covcor.md).

## See also

[`variance()`](https://rodrigozepeda.github.io/deltapif/reference/covcor.md),
[`standard_deviation()`](https://rodrigozepeda.github.io/deltapif/reference/covcor.md)

## Examples

``` r
my_pif <- pif(p = 0.085, p_cft = function(p) 0.85 * p, beta = log(1.59),
              var_p = 0.01^2, var_beta = 0.165^2, quiet = TRUE)
variance_components(my_pif)
#>            p         beta      weights 
#> 6.439931e-07 9.199774e-06 0.000000e+00 

# Share of the variance due to each source
variance_components(my_pif) / variance(my_pif)
#>         p      beta   weights 
#> 0.0654214 0.9345786 0.0000000 
```
