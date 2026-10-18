# Partial derivatives of PIF

Calculates the partial derivatives of a potential impact fraction with
respect to the parameters `p` or `beta`.

## Usage

``` r
deriv_pif_p(p, p_cft, rr, mu_obs = NULL, mu_cft = NULL, p_cft_jacobian = NULL)

deriv_pif_beta(p, p_cft, rr, rr_link_deriv_vals, mu_obs = NULL, mu_cft = NULL)
```

## Arguments

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels. Either a numeric vector, for a fixed
  target that does not depend on the observed prevalence, or a function
  of `p` returning the counterfactual prevalence, for scenarios defined
  relative to the current exposure (e.g. `function(p) 0.85 * p` for a
  15% reduction). In the latter case the uncertainty in `p` is also
  propagated through the counterfactual. See the section *Fixed and
  relative counterfactuals*.

- rr:

  The relative risk for each of the exposure levels.

- mu_obs:

  The average value of the relative risk in the observed population.

- mu_cft:

  The average value of the counterfactual relative risk in the
  population.

- p_cft_jacobian:

  (Optional) Jacobian of the counterfactual prevalence `p_cft` with
  respect to `p` when the counterfactual is a function of the observed
  prevalence (for example a relative reduction `p_cft = 0.85 * p`).
  Entry `[j, i]` is the derivative of `p_cft[j]` with respect to `p[i]`.
  `NULL` (default) assumes `p_cft` is a fixed target.

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
\text{RR}(\beta) - 1\big) \$\$ If the counterfactual is a function of
the prevalence, \\p\_{\*} = g(p)\\, the total derivative adds the term
\\- J_g^{\top} \big( \text{RR}(\beta) - 1\big) / \mu^{\text{obs}}\\
where \\J_g\\ is the Jacobian of \\g\\ evaluated at \\p\\. The partial
derivative of `PIF` with respect to `beta` is: \$\$ \dfrac{\partial
\textrm{PIF}}{\partial \beta} = \Bigg( \dfrac{ \mu^{\text{cft}} \cdot
p - \mu^{\text{obs}} \cdot p\_{\*} }{ \Big( \mu^{\text{obs}}\Big)^2 }
\Bigg)\odot \text{RR}^{\prime}(\beta) \$\$ with \\\odot\\ representing
the Hadamard (elementwise) product.
