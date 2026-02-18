# Potential Impact Fraction related classes

Objects for handling potential impact fractions for a categorical
exposure considering an observed prevalence of `p` and a relative risk
(or relative risk parameter) of `beta`.

## Usage

``` r
pif_class(
  pif = integer(0),
  variance = integer(0),
  conf_level = integer(0),
  label = character(0),
  type = "PIF",
  link = function() NULL,
  link_inv = function() NULL,
  link_deriv = function() NULL
)

cases_class(
  overall_cases,
  pif_obj,
  variance_cases,
  link,
  link_deriv,
  link_inv,
  conf_level
)

pif_atomic_class(
  p,
  p_cft,
  beta,
  var_p,
  var_beta,
  rr_link,
  rr_link_deriv,
  link,
  link_deriv,
  link_inv,
  conf_level,
  type,
  label,
  upper_bound_p,
  upper_bound_beta
)

pif_global_ensemble_class(
  pif_list,
  weights,
  var_weights,
  var_pif_weights,
  pif_transform,
  pif_deriv_transform,
  pif_inverse_transform,
  link,
  link_inv,
  link_deriv,
  conf_level = 0.95,
  label
)

pif_total_class(
  pif_list,
  weights,
  var_weights,
  var_pif_weights,
  link,
  link_inv,
  link_deriv,
  conf_level = 0.95,
  label
)

pif_ensemble_class(
  pif_list,
  weights,
  var_weights,
  var_pif_weights,
  link,
  link_inv,
  link_deriv,
  conf_level = 0.95,
  label
)
```

## Arguments

- pif:

  Potential Impact Fraction estimate

- variance:

  variance estimate for the potential impact fraction (i.e. for `pif`)

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- label:

  Character identifier for the impact fraction. This is for

- type:

  Character either Potential Impact Fraction (`PIF`) or Population
  Attributable Fraction (`PAF`)

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_inv:

  The inverse of `link`. For example if `link` is `logit` this should be
  `inv_logit`.

- link_deriv:

  The derivative of `link`. For example if `link` is `logit` this should
  be `deriv_logit` (i.e. `function(pif) 1 / (pif * (1 - pif))`).

- overall_cases:

  Overall number of cases so that attributable cases is given by
  `overall_cases * pif_obj`

- pif_obj:

  A `pif_class` to calculate the percentage of cases with.

- variance_cases:

  Variance of the cases estimate.

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

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

  Derivative of the link function for the relative risk. The constructor
  tries to build it automatically from `rr_link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- upper_bound_p:

  Whether the values for the `p` component of the link_variance should
  be approximated by an upper bound.

- upper_bound_beta:

  Whether the values for the `beta` component of the link_variance
  should be approximated by an upper bound.

- pif_list:

  A list of potential impact fractions `pif_class` so that the total can
  be computed from it.

- weights:

  weights for calculating the total PIF (respectively PAF) in
  `pif_total`.

- var_weights:

  covariance matrix for the weights when calculating the total PIF
  (respectively PAF) in `pif_total`.

- var_pif_weights:

  covariance matrix with row `i` and column `j` representing the
  covariance between the `i`-th potential impact fraction of the list
  and the `j`-th weight

- pif_transform:

  Transform applied to the `pif` for summation in a
  `pif_global_ensemble_class` (see section below).

- pif_deriv_transform:

  Derivative of the transform applied to the `pif` for summation in a
  `pif_global_ensemble_class` (see section below).

- pif_inverse_transform:

  Inverse of the transform applied to the `pif` for summation in a
  `pif_global_ensemble_class` (see section below).

## Properties of a `pif_class`

Any object that is a `pif_class` contains a potential impact fraction
with intervals estimated as follows: \$\$ \text{CI}\_{\text{Link}} =
\text{link}\big(\text{PIF}\big) \pm
Z\_{\alpha/2}\cdot\sqrt{\textrm{link\\variance}} \$\$ and then
transformed back using the inverse of the link function `inv_link`: \$\$
\text{CI}\_{\text{PIF}} =
\text{link}^{-1}\Big(\text{CI}\_{\text{Link}}\Big) \$\$

The following are the properties of any `pif_class`

- `ci`:

  `numeric(2)` — Lower and upper confidence limits at level
  `conf_level`.

- `link_vals`:

  `numeric` — Entrywise evaluation of the link function at pif:
  `link(pif)`.

- `link_deriv_vals`:

  `character` — Entrywise evaluation of the derivative of the link
  function (`link_deriv`) at pif: `link(pif)`.

- `link_variance`:

  `numeric` - Estimate for the linked potential impact fraction's
  variance: `variance(link(pif))`.

## Properties of a `pif_atomic_class`

The `pif_atomic_class` is a type of `pif_class` that contains enough
information to compute a potential impact fraction through the classic
formula by Walter: \$\$ \textrm{PIF} = \dfrac{ \sum\limits\_{i=1}^N p_i
\text{RR}\_i - \sum\limits\_{i=1}^N p_i^{\text{cft}} \text{RR}\_i }{
\sum\limits\_{i=1}^N p_i \text{RR}\_i } \$\$ where the relative risk is
a function of a parameter \\\beta_i\\ \$\$ \text{RR}\_i =
\text{rr\\link}(\beta_i) \$\$

The `pif_atomic_class` inherits the properties of a `pif_class` as well
as:

- `mu_obs`:

  `numeric` — Average relative risk in the observed population.

- `mu_cft`:

  `numeric` — Average relative risk in the counterfactual population.

- `pif`:

  `numeric` — Estimate of the potential impact fraction.

- `rr_link_deriv_vals`:

  `character` — Entrywise evaluation of the derivative of the link
  function (`link_deriv`) at pif: `link(pif)`.

Confidence intervals are estimated as with any `pif_class`.

## Properties of a `pif_global_ensemble_class`

The `pif_global_ensemble_class` creates a new potential impact fraction
by summing a weighted combination of potential impact fractions. In
general it computes the following expression: \$\$
\textrm{PIF}\_{\text{global}} = g^{-1}\bigg( \sum\limits\_{i = 1}^{N}
g\big(w_i \cdot \textrm{PIF}\_i\big) \bigg) \$\$ where \\g\\ is refered
to as the `pif_transform`, its derivative the `pif_deriv_transform`, and
its inverse `pif_inverse_transform`.

The `pif_global_ensemble_class` inherits the properties of a `pif_class`
as well as:

- `weights`:

  `numeric` - Vector of weights \\w_i\\ for weighting the potential
  impact fraction.

- `var_weights`:

  `numeric` - Covariance matrix for the `weights`

- `pif_transform`:

  `function` - Function \\g\\ with which to transform the impact
  fraction before weighting.

- `pif_deriv_transform`:

  `function` - Derivative of the `pif_transform`.

- `pif_inverse_transform`:

  `function` - Inverse of the `pif_transform`.

- `type`:

  `character` - Whether the quantity represents a `PIF` or a `PAF`

- `coefs`:

  `numeric` - Potential impact fractions used for the global ensemble
  (each of the \\\text{PIF}\_i\\.

- `sum_transformed_weighted_coefs`:

  `numeric` - Sum of the potential impact fractions involved \\\sum
  g(w_i \text{PIF}\_i)\\.

- `pif`:

  `numeric` — Estimate of the potential impact fraction.

- `covariance`:

  `numeric` — Covariance matrix between the potential impact fractions
  in `coefs` (i.e. each entry is\\\text{Cov}(\text{PIF}\_i,
  \text{PIF}\_j))\\

- `variance`:

  `numeric` — Estimate for the variance of `pif`.

Confidence intervals are estimated as with any `pif_class`.

## Properties of a `pif_total_class`

A `pif_total_class` estimated the potential impact fraction of the
weighted sum of fractions from different (disjoint) dispopulations: \$\$
\textrm{PIF}\_{Total} = \sum\limits\_{i = 1}^{N} w_i \cdot
\textrm{PIF}\_i \$\$ with \\w_i\\ representing the proportions of
individuals in each category. This is a type of
`pif_global_ensemble_class` with `pif_transform = identity`.

## Properties of a `pif_ensemble_class`

The ensemble potential impact fraction (representing different relative
risks) for the same outcome is given by the weighted product: \$\$
\textrm{PIF}\_{Ensemble} = 1 - \prod\limits\_{i = 1}^{N} \Big(1 - w_i
\textrm{PIF}\_i\Big) \$\$

However it can be transformed into a `pif_global_ensemble_class` by
taking the log-complement: \$\$ \ln\Big(1 -
\textrm{PIF}\_{Ensemble}\Big) = \sum\limits\_{i = 1}^{N} \ln\big(1 - w_i
\textrm{PIF}\_i\big) \$\$ hence it is a `pif_global_ensemble_class` with
`pif_transform = log_complement`.
