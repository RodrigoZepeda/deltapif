# Weighted Adjusted Fractions

Calculates the weighted adjusted potential impact fractions (or
population attributable fractions). Each individual fraction
\\\widehat{\text{PIF}}\_i\\ is rescaled proportionally so that together
they are consistent with the ensemble fraction.

## Usage

``` r
weighted_adjusted_fractions(
  pif1,
  ...,
  weights = NULL,
  pif_total_link = "log-complement",
  pif_total_link_inv = NULL,
  pif_total_link_deriv = NULL,
  pif_ensemble_link = "identity",
  pif_ensemble_link_inv = NULL,
  pif_ensemble_link_deriv = NULL,
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

- pif1:

  A potential impact fraction (class `pif_class`). This is the first of
  the individual fractions to be adjusted.

- ...:

  The remaining potential impact fractions (class `pif_class`). All
  fractions must be of the same type (all `PIF` or all `PAF`).

- weights:

  Weights for the ensemble (`pif_ensemble`). Passed directly to
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  /
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md).
  Defaults to `NULL` (equal weights of 1 for each fraction).

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

A named list of `pif_class` objects, one per input fraction, where each
element is \\\widehat{\text{PIF}}\_i^{\text{adj}}\\. Names are taken
from the `label` slots of the input fractions.

## Formula

The weighted adjusted fraction for the \\i\\-th exposure is: \$\$
\widehat{\text{PIF}}\_i^{\text{adj}} =
\frac{\widehat{\text{PIF}}\_i}{\sum\_{j=1}^{n} \widehat{\text{PIF}}\_j}
\cdot \widehat{\text{PIF}}\_{\text{Ensemble}} \$\$

Using the log-transform, the variance is: \$\$ \begin{aligned}
\operatorname{Var}\\\Big\[\ln \widehat{\text{PIF}}\_i^{\text{adj}}\Big\]
&= \operatorname{Var}\\\left\[\ln \widehat{\text{PIF}}\_i\right\] +
\operatorname{Var}\\\left\[\ln \textstyle\sum_j
\widehat{\text{PIF}}\_j\right\] + \operatorname{Var}\\\left\[\ln
\widehat{\text{PIF}}\_{\text{Ensemble}}\right\] \\ &\quad + 2\Bigg\[
\operatorname{Cov}\\\Big(\ln \widehat{\text{PIF}}\_i, \ln
\widehat{\text{PIF}}\_{\text{Ensemble}}\Big) -
\operatorname{Cov}\\\Big(\ln \widehat{\text{PIF}}\_i, \ln
\textstyle\sum_j \widehat{\text{PIF}}\_j\Big) -
\operatorname{Cov}\\\Big(\ln \widehat{\text{PIF}}\_{\text{Ensemble}},
\ln \textstyle\sum_j \widehat{\text{PIF}}\_j\Big) \Bigg\] \end{aligned}
\$\$

where each log-covariance is approximated via the delta method: \$\$
\operatorname{Cov}\\\big(\ln X, \ln Y\big) \approx
\frac{\operatorname{Cov}(X, Y)}{X \cdot Y} \$\$

The covariances \\\operatorname{Cov}(\widehat{\text{PIF}}\_i, \cdot)\\
are computed by
[`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md)
applied to the individual fraction and the internally constructed
sum/ensemble objects, which automatically propagates the full
uncertainty structure already embedded in those objects.

## Internal construction

Internally the function builds:

- A **sum** fraction via
  [`pif_total_class()`](https://rodrigozepeda.github.io/deltapif/reference/classes.md)
  with weights all equal to `1`, so its `pif` slot equals \\\sum_i
  \widehat{\text{PIF}}\_i\\.

- An **ensemble** fraction via
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  /
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md).

All cross-covariances are then computed automatically through the
recursive
[`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md)
machinery that already handles atomic, total, and ensemble fractions.

## See also

[`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
[`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
[`cov_total_pif()`](https://rodrigozepeda.github.io/deltapif/reference/cov_total_pif.md),
[`weighted_adjusted_paf()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_paf.md),
[`weighted_adjusted_pif()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted_pif.md)

## Examples

``` r
paf_lead <- paf(0.2, 2.2, quiet = TRUE, var_p = 0.001, label = "Lead")
paf_rad  <- paf(0.1, 1.2, quiet = TRUE, var_p = 0.0001, label = "Radiation")

# Adjusted fractions (covariances computed automatically)
adj <- weighted_adjusted_fractions(paf_lead, paf_rad)
adj$Lead
#> 
#> ── Population Attributable Fraction: [Lead_adj] ──
#> 
#> PAF = 52.726% [95% CI: 46.604% to 58.847%]
#> standard_deviation(paf %) = 3.123
adj$Radiation
#> 
#> ── Population Attributable Fraction: [Radiation_adj] ──
#> 
#> PAF = 16.116% [95% CI: 13.560% to 18.671%]
#> standard_deviation(paf %) = 1.304

# With correlated prevalences
adj2 <- weighted_adjusted_fractions(paf_lead, paf_rad, var_p = 0.0001)
adj2$Lead
#> 
#> ── Population Attributable Fraction: [Lead_adj] ──
#> 
#> PAF = 52.726% [95% CI: 44.538% to 60.913%]
#> standard_deviation(paf %) = 4.177
```
