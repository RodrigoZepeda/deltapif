# Covariance between two atomic pifs

Calculates the approximate covariance between two atomic pif class
expressions using the delta method.

## Usage

``` r
cov_atomic_pif(pif1, pif2, var_p = NULL, var_beta = NULL)
```

## Arguments

- pif1:

  A `pif_atomic_class` object.

- pif2:

  A second `pif_atomic_class` object.

- var_p:

  covariance matrix for the prevalences in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

- var_beta:

  covariance matrix for the parameter `beta` in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

## Computation

The `cov_atomic_pif` computes the approximate covariance between two
potential impact fractions \\\text{PIF}\_1\\ and \\\text{PIF}\_2\\ as
follows: \$\$ \text{Cov}\Big(\text{PIF}\_1, \text{PIF}\_2 \Big) \approx
\Bigg( \dfrac{\partial \text{PIF}\_1}{\partial p_1}^{\top} \Sigma_p
\dfrac{\partial \text{PIF}\_2}{\partial p_2} + \dfrac{\partial
\text{PIF}\_1}{\partial \beta_1}^{\top} \Sigma\_{\beta} \dfrac{\partial
\text{PIF}\_2}{\partial \beta_2} \Bigg) \$\$ where \\\Sigma_p\\ is
`var_p` and \\\Sigma\_{\beta}\\ is `var_beta`. The parameters \\p_1\\
and \\\beta_1\\ refer to the prevalence (`p`) and relative risk exposure
parameter `beta` of the first potential impact fraction `pif1` and
\\p_2\\ and \\\beta_2\\ to the respective parameters of `pif2`.

## See also

[`from_parameters_covariance_p_component()`](https://rodrigozepeda.github.io/deltapif/reference/covariance_from_parameters.md)
to calculate the covariance from first principles,
[cov_ensemble_atomic](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_atomic.md)
and
[`cov_ensemble_weights()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_weights.md)
for the covariance between an ensemble fraction and an atomic fraction
and the covariance between the weights of an ensemble fraction and
another fraction.

## Examples

``` r
if (FALSE) { # \dontrun{
p1 <- pif(p = 0.5, p_cft = 0.1, beta = 0.2, 0.003, 0.04, label = "Population 1")
p2 <- pif(p = 0.4, p_cft = 0.1, beta = 0.2, 0.003, 0.04, label = "Population 2")

#Positive covariance as program automatically detects same beta and p
cov_atomic_pif(p1, p2)

#Zero covariance as they have no covariates in common
p3 <- pif(p = 0.6, p_cft = 0.1, beta = 0.12, 0.003, 0.04, label = "Population 3")
cov_atomic_pif(p1, p3)

#Covariance has to be given to specify the covariance between p and/or beta parameters
cov_atomic_pif(p1, p3, var_p = 0.12, var_beta = 0.11)

#Works with parameters of different dimensions with entry i-j being the
#covariance between entry p[i] of pif1 and entry p[j] of pif2
#(same with the betas)
p_dim_3 <- paf(p = c(0.5, 0.2, 0.1), beta = c(0.2, 0.1, 0.4),
  label = "Population 1", quiet = TRUE)
p_dim_2 <- pif(p = c(0.4, 0.1), p_cft = c(0.1, 0.05), beta = c(0.3, 0.8),
  label = "Population 2", quiet = TRUE)
cov_atomic_pif(p_dim_2, p_dim_3,
  var_p = matrix(c(0.1, 0.2, 0.4, 0.5, 0.6, 0.08), ncol = 3))
} # }
```
