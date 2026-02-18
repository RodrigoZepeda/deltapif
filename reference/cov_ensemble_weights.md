# Covariance between a pif and a weight

Calculates the covariance of the weight of a potential impact fraction
of `pif_global_ensemble_class` with a second `pif_global_ensemble_class`
or `pif_atomic_class`.

## Usage

``` r
cov_ensemble_weights(
  pif1,
  pif2,
  var_weights = NULL,
  var_pif_weights = NULL,
  recursive = !is.null(var_pif_weights)
)
```

## Arguments

- pif1:

  A `pif_global_ensemble_class` from which the weight is taken.

- pif2:

  A `pif_global_ensemble_class` or `pif_atomic_class` to compute the
  covariance

- var_weights:

  Covariance matrix between the weights of `pif1` and the weights of
  `pif2`. Entry `var_weights[i,j]` is the covariance between the `i`th
  weight of `pif1` and the `j`th weight of `pif2`. This refers to the
  term \\\operatorname{Cov}\Big( \hat{q}\_i,\hat{w}\_j\Big)\\ in the
  equation before.

- var_pif_weights:

  Covariance structure between the potential impact fractions in `pif1`
  and the weights in `pif2`. Entry `var_pif_weights[i,j]` is the
  covariance between the `i`th fraction of `pif1` and the `j`th weight
  of `pif2`. This refers to the term
  \\\operatorname{Cov}\Big(\widehat{\textrm{PIF}}\_{A,i},
  \hat{w}\_j\Big)\\ in the equation before.

## Note

The model currently works under the assumption that all
`pif_atomic_class` are independent from weights. Hence will return `0`
if `pif2` is an atomic pif.

## Formula

Given a `pif_global_ensemble`: \$\$ \widehat{\textrm{PIF}}\_{1} =
g^{-1}\Bigg(\sum\limits\_{i=1}^{M_1} g(\hat{q}\_i \cdot
\widehat{\textrm{PIF}}\_{1,i})\Bigg) \$\$ and a second
`pif_global_ensemble`: \$\$ \widehat{\textrm{PIF}}\_{2} =
h^{-1}\Bigg(\sum\limits\_{j=1}^{M_2} h(\hat{w}\_j \cdot
\widehat{\textrm{PIF}}\_{2,j})\Bigg) \$\$ This function computes the
covariance between the first impact fraction and the weights of the
second fraction. This returns a vector where the `j`th entry is given
by: \$\$ \operatorname{Cov}(\widehat{\textrm{PIF}}\_A, \hat{w}\_j)
\approx
\frac{1}{g'\big(\widehat{\textrm{PIF}}\_A\big)}\sum\limits\_{i=1}^{M_1}g'(\hat{q}\_i
\cdot \widehat{\textrm{PIF}}\_{A,i}) \Bigg\[ \hat{q}\_i
\operatorname{Cov}\Big( \widehat{\textrm{PIF}}\_{A,i}, \hat{w}\_j\Big) +
\widehat{\textrm{PIF}}\_{A,i} \operatorname{Cov}\Big( \hat{q}\_i,
\hat{w}\_j\Big) \Bigg\] \$\$ where \\\widehat{\textrm{PIF}}\_{1,:} =
(\widehat{\textrm{PIF}}\_{1,1}, \widehat{\textrm{PIF}}\_{1,2}, \dots,
\widehat{\textrm{PIF}}\_{1,M_1})^{\top}\\

## See also

[`cov_atomic_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_atomic_pif.md),
[`cov_ensemble_atomic()`](https://rodrigozepeda.github.io/deltapif/reference/cov_ensemble_atomic.md)
[`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md)

## Examples

``` r
if (FALSE) { # \dontrun{
p1 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
  var_beta = 0, label ="1")
p2 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
  var_beta = 0, label = "2")
var_weights <- matrix(c(0.1, 0.01, -0.03, 0.01, 0.4, 0.02, -0.03, 0.02, 0.01), ncol = 3)
pt1 <- paf_total(p1, p2, weights = c(0.1, 0.9), var_weights = var_weights[1:2,1:2],
  label ="total")

#This gives 0 as the assumption is that atomic pifs don't correlate with anything
cov_ensemble_weights(p1, p2, var_weights = NULL, var_pif_weights = NULL)

#This returns the correlation of pt1 with its weights
cov_ensemble_weights(pt1, pt1, var_weights = NULL, var_pif_weights = NULL)

#Changes if it changes to ensemble
p3 <- paf(c(0.1, 0.2), beta = c(0.01, 0.14), var_p = diag(c(0.01, 0.011)),
  var_beta = 0, label = "3")
p4 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
  var_beta = 0, label = "4")
p5 <- paf(c(0.2, 0.3), beta = c(0.01, 0.15), var_p = diag(c(0.06, 0.052)),
  var_beta = 0, label = "5")
pe1 <- paf_ensemble(p3, p4, p5, weights = c(0.1, 0.5, 0.4),
  var_weights = var_weights, label = "ensemble")

#This returns the correlation of pt1 with its weights
cov_ensemble_weights(pe1, pe1, var_weights = NULL, var_pif_weights = NULL)

#You can also specify var_weights as either a matrix
var_weights     <- matrix(c(0.1, 0.2, 0.1, 0.3, 0.4, -0.2), ncol = 2)
var_pif_weights <- matrix(c(0.01, 0.03, 0.021, -0.02, 0.004, -0.5), ncol = 2)
cov_ensemble_weights(pe1, pt1, var_weights = var_weights,
  var_pif_weights = var_pif_weights)

#or a covariance structure:
var_weights2     <- default_weight_covariance_structure2(pe1, pt1)
var_weights2[[pe1@label]][[pt1@label]] <- var_weights
var_weights2[[pt1@label]][[pe1@label]] <- var_weights

var_pif_weights2 <- default_weight_pif_covariance_structure2(pe1, pt1)
var_pif_weights2[[pe1@label]][[pt1@label]] <- var_pif_weights
var_pif_weights2[[pt1@label]][[pe1@label]] <- var_pif_weights
cov_ensemble_weights(pe1, pt1, var_weights = var_weights2,
  var_pif_weights = var_pif_weights2)

} # }

```
