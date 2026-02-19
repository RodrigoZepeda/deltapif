# Weighted Adjusted PAF

Convenience wrapper around
[`weighted_adjusted_fractions()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_fractions.md)
for population attributable fractions (PAF).

## Usage

``` r
weighted_adjusted_paf(
  paf1,
  ...,
  weights = NULL,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  conf_level = 0.95,
  label_ensemble = NULL,
  label_sum = NULL,
  quiet = FALSE
)

weighted_adjusted_pif(
  pif1,
  ...,
  weights = NULL,
  var_weights = 0,
  var_pif_weights = NULL,
  var_p = NULL,
  var_beta = NULL,
  pif_total_link = "log-complement",
  pif_total_link_inv = NULL,
  pif_total_link_deriv = NULL,
  pif_ensemble_link = "identity",
  pif_ensemble_link_inv = NULL,
  pif_ensemble_link_deriv = NULL,
  conf_level = 0.95,
  label_ensemble = NULL,
  label_sum = NULL,
  quiet = FALSE
)
```

## Arguments

- paf1:

  A population attributable fraction

- ...:

  The remaining potential impact fractions (class `pif_class`). All
  fractions must be of the same type (all `PIF` or all `PAF`).

- weights:

  Weights for the ensemble (`pif_ensemble`). Passed directly to
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  /
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md).
  Defaults to `NULL` (equal weights of 1 for each fraction).

- var_weights:

  Covariance structure for `weights`. Passed directly to
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  /
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md).
  Defaults to `0`.

- var_pif_weights:

  Covariance matrix between individual fractions and `weights`. Passed
  to
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  /
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md).
  Defaults to `NULL`.

- var_p:

  Covariance matrix for the prevalence parameters `p` across all
  fractions. Passed to
  [`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md)
  when computing cross-covariances. Defaults to `NULL`.

- var_beta:

  Covariance matrix for the relative-risk parameters `beta` across all
  fractions. Passed to
  [`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md)
  when computing cross-covariances. Defaults to `NULL`.

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- label_ensemble:

  Character label for the internally constructed ensemble fraction.
  Defaults to `NULL` (auto-generated).

- label_sum:

  Character label for the internally constructed sum (total) fraction.
  Defaults to `NULL` (auto-generated).

- quiet:

  Whether to show messages.

- pif1:

  A potential impact fraction

- pif_total_link:

  Link to pass to `pif_total` in the denominator of the adjustment

- pif_total_link_inv:

  Inverse of the link to pass to `pif_total` in the denominator of the
  adjustment

- pif_total_link_deriv:

  Derivative of the link to pass to `pif_total` in the denominator of
  the adjustment

- pif_ensemble_link:

  Link to pass to `pif_ensemble` in the numerator of the adjustment

- pif_ensemble_link_inv:

  Inverse of the link to pass to `pif_ensemble` in the numerator of the
  adjustment

- pif_ensemble_link_deriv:

  Derivative of the link to pass to `pif_ensemble` in the numerator of
  the adjustment

## Value

A named list of `pif_class` objects, one per input fraction, each being
the weighted adjusted PAF (or PIF, respectively).

## See also

[`weighted_adjusted_fractions()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_fractions.md),
[`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
[`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)

[`weighted_adjusted_fractions()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_fractions.md),
[`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)

## Examples

``` r
paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001, label = "Lead")
paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001, label = "Radiation")

weighted_adjusted_paf(paf_lead, paf_rad)
#> $Lead
#> 
#> ── Population Attributable Fraction: [Lead_adj] ──
#> 
#> PAF = 52.726% [95% CI: 46.604% to 58.847%]
#> standard_deviation(paf %) = 3.123
#> 
#> $Radiation
#> 
#> ── Population Attributable Fraction: [Radiation_adj] ──
#> 
#> PAF = 16.116% [95% CI: 13.560% to 18.671%]
#> standard_deviation(paf %) = 1.304
#> 

pif_lead <- pif(0.2, p_cft = 0.1, beta = log(2.2), quiet = TRUE,
                var_p = 0.001, label = "Lead")
pif_rad  <- pif(0.1, p_cft = 0.05, beta = log(1.2), quiet = TRUE,
                var_p = 0.0001, label = "Radiation")

weighted_adjusted_pif(pif_lead, pif_rad)
#> $Lead
#> 
#> ── Potential Impact Fraction: [Lead_adj] ──
#> 
#> PIF = 9.591% [95% CI: 4.226% to 14.956%]
#> standard_deviation(pif %) = 2.737
#> 
#> $Radiation
#> 
#> ── Potential Impact Fraction: [Radiation_adj] ──
#> 
#> PIF = 0.972% [95% CI: 0.595% to 1.349%]
#> standard_deviation(pif %) = 0.192
#> 
```
