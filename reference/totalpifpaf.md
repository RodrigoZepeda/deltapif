# Combine Potential Impact Fractions and Population Attributable Fractions

Combine potential impact fractions or the population attributable
fractions to either generate the total fraction from the fractions of
subpopulations (`pif_total/paf_total`) or the ensemble fraction of a
population from different (independent) exposures.

## Usage

``` r
paf_total(
  paf1,
  ...,
  weights,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  conf_level = 0.95,
  link = "log-complement",
  link_inv = NULL,
  link_deriv = NULL,
  quiet = FALSE,
  label = NULL
)

pif_total(
  pif1,
  ...,
  weights,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  conf_level = 0.95,
  link = "log-complement",
  link_inv = NULL,
  link_deriv = NULL,
  quiet = FALSE,
  label = NULL,
  is_paf = FALSE,
  weights_sum_to_1 = TRUE
)

paf_ensemble(
  paf1,
  ...,
  weights = NULL,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  link = "identity",
  link_inv = NULL,
  link_deriv = NULL,
  conf_level = 0.95,
  label = NULL,
  quiet = FALSE
)

pif_ensemble(
  pif1,
  ...,
  weights = NULL,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  link = "identity",
  link_inv = NULL,
  link_deriv = NULL,
  conf_level = 0.95,
  quiet = FALSE,
  label = NULL,
  is_paf = FALSE
)
```

## Arguments

- paf1:

  A population attributable fraction (class `pif_class`)

- ...:

  The remaining potential impact fractions or (respectively) population
  attributable fractions.

- weights:

  A vector containing the proportion of the population for each of the
  categories (for each of the pifs given).

- var_weights:

  link_covariance structure for the `weights`. Can be `0` (default) if
  the weights are not random, a vector if only the link_variances of the
  weights are available or a link_covariance matrix.

- var_pif_weights:

  covariance matrix with row `i` and column `j` representing the
  covariance between the `i`-th potential impact fraction of the list
  and the `j`-th weight

- var_p, var_beta:

  (Optional) Covariance matrices between the estimates of `p`
  (respectively `beta`) of different fractions. Row and column names are
  the identifiers given to the parameters as named lists in
  [`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)
  or
  [`paf()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md);
  a vector-valued parameter repeats its identifier once per element. The
  diagonal blocks must equal the variances given to each fraction.
  Fractions that share an identifier are always correlated, so these
  matrices are only needed for correlations between different estimates
  (e.g. relative risks from the same meta-analysis).

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- quiet:

  Whether to show messages.

- label:

  Character identifier for the impact fraction. This is for

- pif1:

  A potential impact fraction (class `pif_class`)

- is_paf:

  Whether the computed quantity is a population attributable fraction or
  not

- weights_sum_to_1:

  Boolean flag indicating if the weights sum to 1 (normalized weights)
  or if they are not (unnormalized).

## Value

A `pif_class` object containing the estimate for the total fraction or
the ensemble fraction that aggregates several fractions (either
potential impact or population attributable).

## Total potential impact fraction

Assuming the overall population can be subdivided into \\N\\ distinct
subpopulations each of them with a different potential impact fraction
(or population attributable fraction) we can estimate the total
population attributable fraction or potential impact fraction of the
whole population as:

\$\$ \text{PIF}\_{\text{Total}} = \sum\limits\_{i = 1}^{N} \pi_i \cdot
\text{PIF}\_i \$\$

where each \\\text{PIF}\_i\\ corresponds to the potential impact
fraction of the i-th subpopulation and \\\pi_i\\ correspond to the
proportion of the total population occupied by \\\text{PIF}\_i\\. The
weights are such that \\\sum\_{i=1}^{N} \pi_i = 1\\.

## Ensemble potential impact fraction

If a population is exposed to \\K\\ different independent risk factors
then the ensemble impact fraction of the combination of those factors
can be written as:

\$\$ \text{PIF}\_{\text{Ensemble}} = 1 - \prod\limits\_{\ell = 1}^{K}
\Big(1 - \pi\_{\ell} \cdot \text{PIF}\_{\ell}\Big) \$\$

where each \\\text{PIF}\_{\ell}\\ corresponds to the potential impact
fraction of the \\\ell\\-th risk factor for the same population.

## Examples

``` r
#Potential impact fraction for women
pif_women <- pif(0.32, 0.1, 1.2, quiet = TRUE, var_p = 0.1)

#Potential impact fraction for men
pif_men <- pif(0.27, 0.1, 1.3, quiet = TRUE, var_p = 0.1)

#Population potential impact fraction with 49% men and 51% women
pif_total(pif_men, pif_women, weights = c(0.49, 0.51), link = "logit")
#> 
#> ── Potential Impact Fraction: [deltapif-114269487842917] ──
#> 
#> PIF = 27.862% [95% CI: 3.823% to 78.962%]
#> standard_deviation(pif %) = 23.319
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 26.372% (sd %: 36.119) --- [deltapif-106199914454589]
#> • 29.294% (sd %: 29.772) --- [deltapif-256440833799741]
#> ────────────────────────────────────────────────────────────────────────────────

#Population attributable  fraction for women
paf_women <- paf(0.32, 1.3, quiet = TRUE, var_p = 0.1)

#Population attributable  fraction for men
paf_men <- paf(0.27, 1.3, quiet = TRUE, var_p = 0.1)
paf_total(paf_men, paf_women, weights = c(0.49, 0.51), link = "logit")
#> 
#> ── Population Attributable Fraction: [deltapif-0415406398921884] ──
#> 
#> PAF = 44.018% [95% CI: 15.026% to 77.759%]
#> standard_deviation(paf %) = 18.760
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 41.884% (sd %: 28.509) --- [deltapif-0558514421664496]
#> • 46.068% (sd %: 24.552) --- [deltapif-00751059001682063]
#> ────────────────────────────────────────────────────────────────────────────────

# Calculate the ensemble from lead and radiation exposure
paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001)
paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001)
pif_ensemble(paf_lead, paf_rad)
#> 
#> ── Population Attributable Fraction: [deltapif-138291086057034] ──
#> 
#> PAF = 68.841% [95% CI: 62.782% to 74.901%]
#> standard_deviation(paf %) = 3.092
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 61.612% (sd %: 3.740) --- [deltapif-120376783938691]
#> • 18.832% (sd %: 1.529) --- [deltapif-143127077688999]
#> ────────────────────────────────────────────────────────────────────────────────

# Fractions that use the same relative risk are correlated. Declare it
# by giving `beta` the same identifier (name) in a list:
paf_women <- paf(0.32, beta = list("rr_meta" = log(1.3)), var_beta = 0.01,
                 var_p = 0.001, label = "Women")
paf_men   <- paf(0.27, beta = list("rr_meta" = log(1.3)), var_beta = 0.01,
                 var_p = 0.001, label = "Men")
paf_total(paf_men, paf_women, weights = c(0.49, 0.51))
#> 
#> ── Population Attributable Fraction: [deltapif-217757602827549] ──
#> 
#> PAF = 8.139% [95% CI: 1.465% to 14.361%]
#> standard_deviation(paf %) = 3.287
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 7.493% (sd %: 3.111) --- [Men]
#> • 8.759% (sd %: 3.552) --- [Women]
#> ────────────────────────────────────────────────────────────────────────────────

# Correlations between different estimates are given by identifier
paf_lead <- paf(0.2, beta = list("lead" = log(2.2)), var_beta = 0.04,
                var_p = 0.001, label = "Lead")
paf_rad  <- paf(0.1, beta = list("radiation" = log(1.2)), var_beta = 0.02,
                var_p = 0.0001, label = "Radiation")
cov_rr <- matrix(c(0.04, 0.01, 0.01, 0.02), ncol = 2,
                 dimnames = list(c("lead", "radiation"), c("lead", "radiation")))
paf_ensemble(paf_lead, paf_rad, var_beta = cov_rr)
#> 
#> ── Population Attributable Fraction: [deltapif-0876777326556075] ──
#> 
#> PAF = 20.936% [95% CI: 7.889% to 33.983%]
#> standard_deviation(paf %) = 6.657
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 19.355% (sd %: 6.233) --- [Lead]
#> • 1.961% (sd %: 1.642) --- [Radiation]
#> ────────────────────────────────────────────────────────────────────────────────

# Totals and ensembles can be combined
pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001)
pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001)
pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women)
pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001)
pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001)
pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men)

pif_total(pif_men, pif_women, weights = c(0.49, 0.51))
#> 
#> ── Population Attributable Fraction: [deltapif-0654550251367236] ──
#> 
#> PAF = 75.731% [95% CI: 72.232% to 78.788%]
#> standard_deviation(paf %) = 1.668
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 76.180% (sd %: 2.271) --- [deltapif-0360278296382459]
#> • 75.299% (sd %: 2.435) --- [deltapif-147619690566501]
#> ────────────────────────────────────────────────────────────────────────────────
```
