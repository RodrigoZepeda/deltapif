# Covariance between `pif_global_ensemble_class` and a `pif_atomic`

Calculates the covariance of a potential impact fraction of a
`pif_global_ensemble_class` with a `pif_atomic`

## Usage

``` r
cov_ensemble_atomic(
  pif_ensemble,
  pif_atomic,
  var_p = NULL,
  var_beta = NULL,
  var_pifs = NULL,
  var_pif_weights = NULL,
  recursive = !is.null(var_pifs),
  warning = TRUE
)
```

## Arguments

- pif_ensemble:

  A `pif_global_ensemble_class`

- pif_atomic:

  A `pif_atomic_class`

- var_pifs:

  Covariance vector between the potential impact fractions in
  `pif_ensemble` and `pif_atomic`. This refers to the term
  \\\operatorname{Cov}\Big( \textrm{PIF}\_{A,i},
  \widehat{\textrm{PIF}}\_{B,j}\Big)\\ in the equation below. If set to
  `NULL` its automatically calculated.

- var_pif_weights:

  Covariance vector between the weights in `pif_ensemble` and the
  `pif_atomic`. This refers to the term \\\operatorname{Cov}\Big(
  \hat{q}\_i,\widehat{\textrm{PIF}}\_{B,j}\Big)\\ in the equation below.
  If set to `NULL` its automatically calculated.

## Formula

Given a `pif_global_ensemble`: \$\$ \widehat{\textrm{PIF}}\_{1} =
g^{-1}\Bigg(\sum\limits\_{i=1}^{M_1} g(\hat{q}\_i \cdot
\widehat{\textrm{PIF}}\_{1,i})\Bigg) \$\$ and a second
`pif_global_ensemble`: \$\$ \widehat{\textrm{PIF}}\_{2} =
h^{-1}\Bigg(\sum\limits\_{j=1}^{M_2} h(\hat{w}\_j \cdot
\widehat{\textrm{PIF}}\_{2,j})\Bigg) \$\$ This function computes the
covariance between the first impact fraction and the coefficients of the
second fraction by computing: \$\$
\operatorname{Cov}(\widehat{\textrm{PIF}}\_A,
\widehat{\textrm{PIF}}\_{B,j}) \approx
\frac{1}{g'\big(\widehat{\textrm{PIF}}\_A\big)}\sum\limits\_{i=1}^{M_1}g'(\hat{q}\_i
\cdot \widehat{\textrm{PIF}}\_{A,i}) \Bigg\[ \hat{q}\_i
\operatorname{Cov}\Big( \textrm{PIF}\_{A,i},
\widehat{\textrm{PIF}}\_{B,j}\Big) + \widehat{\textrm{PIF}}\_{A,i}
\operatorname{Cov}\Big( \hat{q}\_i, \widehat{\textrm{PIF}}\_{B,j}\Big)
\Bigg\] \$\$ where \\\widehat{\textrm{PIF}}\_{1,:} =
(\widehat{\textrm{PIF}}\_{1,1}, \widehat{\textrm{PIF}}\_{1,2}, \dots,
\widehat{\textrm{PIF}}\_{1,M_1})^{\top}\\

## See also

[`cov_ensemble_weights()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_weights.md),
[`cov_atomic_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_atomic_pif.md),
[`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md).

## Examples

``` r
if (FALSE) { # \dontrun{
#Works for two atomic fractions
p1 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
  var_beta = diag(c(0.51, 0.3)), label ="1")
p2 <- paf(c(0.2, 0.3), beta = c(0.01, 0.14), var_p = diag(c(0.06, 0.052)),
  var_beta = diag(c(0.51, 0.3)), label = "2")
cov_ensemble_atomic(p1, p2)

#Works for fractions with ensembles
pt1 <- paf_total(p1, p2, weights = c(0.1, 0.9), var_weights = diag(c(0.6, 0.2)),
  label ="total")
p3 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
  var_beta = 0, label = "3")
cov_ensemble_atomic(pt1, p3)

#Changes if it changes to ensemble
p4 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
  var_beta = 0, label = "4")
p5 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
  var_beta = 0, label = "5")
pe1 <- paf_ensemble(p3, p4, p5, weights = c(0.1, 0.5, 0.4),
  var_weights = diag(c(0.6, 0.2, 0.1)), label = "ensemble")

#This returns the correlation of pt1 with its weights
cov_ensemble_atomic(pe1, p1)

#You can also specify var_weights as either a matrix
var_pif_weights <- matrix(c(0.01, 0.03, 0.021), ncol = 1)
cov_ensemble_atomic(pe1, p1, var_pif_weights = var_pif_weights)

#or a covariance structure:
var_pif_weights2     <- default_weight_covariance_structure2(pe1, p1)
var_pif_weights2[[p3@label]][[p1@label]]  <- var_pif_weights[1]
var_pif_weights2[[p1@label]][[p3@label]]  <- var_pif_weights[1]
var_pif_weights2[[p4@label]][[p1@label]]  <- var_pif_weights[2]
var_pif_weights2[[p1@label]][[p4@label]]  <- var_pif_weights[2]
var_pif_weights2[[p5@label]][[p1@label]]  <- var_pif_weights[3]
var_pif_weights2[[p1@label]][[p5@label]]  <- var_pif_weights[3]

cov_ensemble_atomic(pe1, p1, var_pif_weights = var_pif_weights2)

} # }


```
