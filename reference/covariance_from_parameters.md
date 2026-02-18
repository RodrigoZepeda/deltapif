# Covariance components

Calculates the covariance `p` or `beta` component between two potential
impact fractions as well as the complete variance and covariance from
parameters.

## Usage

``` r
from_parameters_covariance_p_component(
  p1,
  p2,
  p1_cft,
  p2_cft,
  rr1,
  rr2,
  var_p,
  upper_bound
)

from_parameters_covariance_beta_component(
  p1,
  p2,
  p1_cft,
  p2_cft,
  rr1,
  rr2,
  rr_link_deriv_vals1,
  rr_link_deriv_vals2,
  var_beta,
  upper_bound
)

from_parameters_pif_covariance(
  p1,
  p2,
  p1_cft,
  p2_cft,
  rr1,
  rr2,
  rr_link_deriv_vals1,
  rr_link_deriv_vals2,
  var_p,
  var_beta,
  upper_bound_p = FALSE,
  upper_bound_beta = FALSE
)

from_parameters_pif_variance(
  p,
  p_cft,
  rr,
  rr_link_deriv_vals,
  var_p,
  var_beta,
  upper_bound_p = FALSE,
  upper_bound_beta = FALSE
)
```

## Arguments

- p1:

  Observed proportion exposed in the first potential impact fraction

- p2:

  Observed proportion exposed in the first potential impact fraction

- p1_cft:

  Counterfactual proportion exposed in the first potential impact
  fraction

- p2_cft:

  Counterfactual proportion exposed in the first potential impact
  fraction

- rr1:

  Relative risk in the first potential impact fraction

- rr2:

  Relative risk in the second potential impact fraction

- var_p:

  Covariance matrix with the entry `var_p[i,j]` corresponding to the
  covariance between `p1[i]` and `p2[j]`.

- upper_bound:

  Whether the variance should be calculated or an upper bound for the
  variance (assuming perfect correlation).

- rr_link_deriv_vals1:

  Values for the derivative of the relative risk function evaluated at
  `theta` (the relative risk parameter).

- rr_link_deriv_vals2:

  Values for the derivative of the relative risk function evaluated at
  `theta` (the relative risk parameter).

- var_beta:

  Estimate of the link_covariance matrix of `beta` where the entry
  `var_beta[i,j]` represents the link_covariance between `beta[i]` and
  `beta[j]`.

- upper_bound_p:

  Whether the values for the `p` component of the link_variance should
  be approximated by an upper bound.

- upper_bound_beta:

  Whether the values for the `beta` component of the link_variance
  should be approximated by an upper bound.

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

- rr:

  The relative risk for each of the exposure levels.

- rr_link_deriv_vals:

  The derivative of the relative risk function `g` with respect to the
  parameter `beta` evaluated at `beta`.

## Value

The covariance component for either `p` or `beta`

## Note

This estimation does not require a `pif_class` object but instead is
built from parameters.

## Formulas

The following represents the covariance components. For `p`:

\$\$ \text{CC}\_{p} = \dfrac{ \partial \textrm{PIF}\_i}{\partial
p_i}^{\top} \textrm{covariance}\big( \hat{p}\_i, \hat{p}\_j\big)
\dfrac{\partial \textrm{PIF}\_j}{\partial p_j} \$\$

and for `beta`:

\$\$ \text{CC}\_{\beta} = \dfrac{ \partial \textrm{PIF}\_i}{\partial
\beta_i}^{\top} \textrm{covariance}\big( \hat{\beta}\_i,
\hat{\beta}\_j\big) \dfrac{\partial \textrm{PIF}\_j}{\partial \beta_j}
\$\$

The total covariance is approximated by the sum: \$\$
\text{Cov}(\text{PIF}\_i, \text{PIF}\_j) \approx \text{CC}\_{p} +
\text{CC}\_{\beta} \$\$

## See also

[derivatives](https://rodrigozepeda.github.io/deltapif/reference/derivatives.md)
for the definition of the derivatives involved in the formulas.
