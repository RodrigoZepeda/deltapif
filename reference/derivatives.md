# Partial derivatives of PIF

Calculates the partial derivatives of a potential impact fraction with
respect to the parameters `p` or `beta`.

## Usage

``` r
deriv_pif_p(p, p_cft, rr, mu_obs = NULL, mu_cft = NULL)

deriv_pif_beta(p, p_cft, rr, rr_link_deriv_vals, mu_obs = NULL, mu_cft = NULL)
```

## Arguments

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

- rr:

  The relative risk for each of the exposure levels.

- mu_obs:

  The average value of the relative risk in the observed population.

- mu_cft:

  The average value of the counterfactual relative risk in the
  population.

- rr_link_deriv_vals:

  The derivative of the relative risk function `g` with respect to the
  parameter `beta` evaluated at `beta`.

## Value

The partial derivative (usually a vector)

## Note

As `p` and `beta` are usually vectors these are vector-valued
derivatives.

## Formulas

The partial derivative of `PIF` with respect to `p` is: \$\$
\dfrac{\partial \textrm{PIF}}{\partial p} =
\dfrac{\mu^{\text{cft}}}{\big(\mu^{\text{obs}}\big)^2} \cdot \big(
\text{RR}(\beta) - 1\big) \$\$ The partial derivative of `PIF` with
respect to `beta` is: \$\$ \dfrac{\partial \textrm{PIF}}{\partial \beta}
= \Bigg( \dfrac{ \mu^{\text{cft}} \cdot p - \mu^{\text{obs}} \cdot
p\_{\*} }{ \Big( \mu^{\text{obs}}\Big)^2 } \Bigg)\odot
\text{RR}^{\prime}(\beta) \$\$ with \\\odot\\ representing the Hadamard
(elementwise) product.
