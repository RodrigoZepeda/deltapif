# Covariance function for a `pif_global_ensemble_class`

Recursively obtains the covariance matrix of a
`pif_global_ensemble_class`

## Usage

``` r
cov_total_pif(
  pif1,
  pif2,
  var_p = NULL,
  var_beta = NULL,
  var_weights = NULL,
  var_pif_weights = NULL,
  var_pifs = NULL,
  warning = TRUE
)
```

## Arguments

- pif1:

  Either a `pif_atomic_class` or a `pif_global_ensemble_class`

- pif2:

  Either a `pif_atomic_class` or a `pif_global_ensemble_class`

- var_p:

  covariance matrix for the prevalences in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

- var_beta:

  covariance matrix for the parameter `beta` in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

- var_weights:

  Covariance between the weights of all the fractions contained in
  `pif1` and all the fractions contained in `pif2`.

## Computation

This computes: \$\$ \operatorname{Cov}(\widehat{\textrm{PIF}}\_{A},
\widehat{\textrm{PIF}}\_{B}) \approx
\frac{1}{g'\big(\widehat{\textrm{PIF}}\_A\big)}\sum\limits\_{i=1}^{M_1}g'(\hat{q}\_i
\cdot \widehat{\textrm{PIF}}\_{A,i}) \Bigg\[ \mathbb{E}\[\hat{q}\_i\]
\operatorname{Cov}\Big( \textrm{PIF}\_{A,i},
\widehat{\textrm{PIF}}\_{B}\Big) + \mathbb{E}\[\textrm{PIF}\_{A,i}\]
\operatorname{Cov}\Big( \hat{q}\_i, \widehat{\textrm{PIF}}\_{B}\Big)
\Bigg\] \$\$ where \\\operatorname{Cov}\Big( \textrm{PIF}\_{A,i},
\widehat{\textrm{PIF}}\_{B}\Big)\\ is computed as in
[`cov_ensemble_atomic()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_atomic.md)
in the case \\\widehat{\textrm{PIF}}\_{A,i}\\ is atomic or with the same
formula in the case it is an ensemble. The expression
\\\operatorname{Cov}\Big( \hat{q}\_i, \widehat{\textrm{PIF}}\_{B}\Big)
\\ is computed as in
[`cov_ensemble_weights()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_weights.md)
for each weight \\q_j\\.

## See also

[`from_parameters_covariance_p_component()`](https://rodrigozepeda.github.io/deltapif/reference/covariance_from_parameters.md),
[`cov_atomic_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_atomic_pif.md),
[`cov_ensemble_atomic()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_atomic.md),
[`cov_ensemble_weights()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_weights.md)
