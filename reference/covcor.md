# Covariance matrix, correlation matrix, variance and standard deviation for potential impact fractions

Computes the covariance (`covariance`) or correlation (`correlation`)
for multiple potential impact fractions and the variance `variance` and
standard deviation `standard_deviation`for a potential impact fractions.

## Usage

``` r
covariance(
  x,
  ...,
  var_p = NULL,
  var_beta = NULL,
  var_weights = NULL,
  var_pif_weights = NULL,
  var_pifs = NULL,
  warning = FALSE
)

variance(x, ...)

standard_deviation(x, ...)

correlation(x, ..., var_p = NULL, var_beta = NULL)
```

## Arguments

- x:

  A potential impact fraction

- ...:

  Multiple additional potential impact fraction objects separated by
  commas.

- var_p:

  covariance matrix for the prevalences in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

- var_beta:

  covariance matrix for the parameter `beta` in both `pif1` and `pif2`.
  Default `NULL` (no correlation).

- var_weights:

  Covariance between the weights of all the fractions contained in
  `pif1` and all the fractions contained in `pif2`.

- var_pif_weights:

  Covariance vector between the weights in `pif_ensemble` and the
  `pif_atomic`. This refers to the term \\\operatorname{Cov}\Big(
  \hat{q}\_i,\widehat{\textrm{PIF}}\_{B,j}\Big)\\ in the equation below.
  If set to `NULL` its automatically calculated.

- var_pifs:

  Covariance vector between the potential impact fractions in
  `pif_ensemble` and `pif_atomic`. This refers to the term
  \\\operatorname{Cov}\Big( \textrm{PIF}\_{A,i},
  \widehat{\textrm{PIF}}\_{B,j}\Big)\\ in the equation below. If set to
  `NULL` its automatically calculated.

- warning:

  Boolean indicating whether to throw a warning if the labels on the
  fractions involved are not unique.

## Examples

``` r
# Get the approximate link_variance of a pif object
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3,
              var_p = 0.1, var_beta = 0.2)
variance(my_pif)
#> [1] 0.07233754

# This is the same as link_covariance with just 1 pIF
covariance(my_pif)
#>            [,1]
#> [1,] 0.07233754

# Calculate the link_covariance between 3 fractions with shared relative risk
beta <- 0.3
var_beta <- 0.1
pif1 <- pif(0.5, 0.2, beta, var_p = 0.5 * (1 - 0.5) / 100, var_beta = var_beta)
pif2 <- pif(0.3, 0.1, beta, var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta)
pif3 <- pif(0.7, 0.3, beta, var_p = 0.7 * (1 - 0.7) / 100, var_beta = var_beta)
covariance(pif1, pif2, pif3)
#>             [,1]        [,2]        [,3]
#> [1,] 0.008789254 0.006486541 0.010220324
#> [2,] 0.006486541 0.005074095 0.007703812
#> [3,] 0.010220324 0.007703812 0.012268945

# The link_covariance between a pif and itself only has the link_variance as entries
covariance(pif1, pif1)
#>             [,1]        [,2]
#> [1,] 0.008789254 0.008789254
#> [2,] 0.008789254 0.008789254

# Or if there is a link_covariance structure between different betas you can specify with
# var_beta in the link_covariance
betas <- c(1.3, 1.2, 1.27)

# link_covariance among all betas
var_beta <- matrix(c(
  1.0000000, -0.12123053, 0.35429369,
  -0.1212305, 1.00000000, -0.04266409,
  0.3542937, -0.04266409, 1.00000000
), byrow = TRUE, ncol = 3)
pif1 <- pif(0.5, 0.2, betas[1], var_p = 0.5 * (1 - 0.5) / 100, var_beta = var_beta[1, 1])
pif2 <- pif(0.3, 0.1, betas[2], var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta[2, 2])
pif3 <- pif(0.3, 0.1, betas[3], var_p = 0.3 * (1 - 0.3) / 100, var_beta = var_beta[3, 3])
covariance(pif1, pif2, pif3, var_beta = var_beta)
#>              [,1]          [,2]          [,3]
#> [1,]  0.042197697 -5.651801e-03  1.629740e-02
#> [2,] -0.005651801  5.536138e-02 -9.642887e-05
#> [3,]  0.016297402 -9.642887e-05  5.410105e-02

# Compute the correlation
correlation(pif1, pif2, pif3, var_beta = var_beta)
#>            [,1]         [,2]         [,3]
#> [1,]  1.0000000 -0.116933525  0.341091700
#> [2,] -0.1169335  1.000000000 -0.001761979
#> [3,]  0.3410917 -0.001761979  1.000000000
```
