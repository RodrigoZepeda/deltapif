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
```

## Arguments

- paf1:

  A population attributable fraction (class `pif_class` with
  `type = "PAF"`).

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

## Value

A named list of `pif_class` objects (type `"PAF"`), one per input
fraction, each being the weighted adjusted PAF.

## See also

[`weighted_adjusted_fractions()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_fractions.md),
[`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)

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
```
