# Changelog

## deltapif 0.5.0

### Breaking changes

- **Shared parameters are declared explicitly.** `p` and `beta` in
  [`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)
  and
  [`paf()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)
  can be given as a named list with one element, e.g.
  `beta = list("smoking_dementia" = 0.46)` or
  `p = list("447" = c(0.12, 0.87))`. Fractions using the same name
  (identifier) share the estimate and are correlated; using one
  identifier with different values or variances is an error. Parameters
  given as plain numbers are now treated as independent. Previously,
  fractions whose `beta` and `var_beta` (or `p` and `var_p`) were
  numerically identical were silently assumed to share the estimate,
  which missed shared estimates that were rounded differently and
  invented sharing between independent estimates that happened to
  coincide. Code that relied on this should give the shared parameters a
  common identifier.

- [`weighted_adjusted_paf()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md)
  and
  [`weighted_adjusted_pif()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md)
  build their confidence intervals on the log scale by default (new
  argument `link`), consistent with the log-scale variance of the
  adjusted fractions; use `link = "identity"` for the previous symmetric
  intervals.

- The `dementiacov` dataset is replaced by `dementiacor`. `dementiacov`
  was built from the pairwise tetrachoric correlations between the
  dementia risk factors (exposures) in ARIC (eTable 2 of Lee et al.)
  multiplied by the standard errors of the log relative risks, so it was
  not a covariance between relative-risk estimates. `dementiacor`
  contains the published correlations themselves, documented as
  correlations between exposures.

### New features

- [`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)
  accepts the counterfactual prevalence `p_cft` as a function of `p`
  (e.g. `function(p) 0.85 * p`) for scenarios defined relative to the
  current exposure. The uncertainty in `p` is then propagated through
  the counterfactual (numerically or with the new argument
  `p_cft_deriv`). Numeric `p_cft` values are still treated as fixed
  targets.

- [`pif_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  their `paf` versions and the weighted adjusted fractions accept
  `var_p` and `var_beta`: covariance matrices between different
  estimates whose row and column names are the identifiers of the
  parameters (e.g. relative risks from the same meta-analysis).

- New
  [`variance_components()`](https://rodrigozepeda.github.io/deltapif/reference/variance_components.md)
  splits the variance of a fraction into the parts due to the
  prevalences (`p`), the relative risks (`beta`) and the weights.

- Compatibility with the development version of S7 (\> 0.2.2):
  constructors of the fraction classes build their parent instances, and
  the package uses S7’s load, unload and build hooks when available
  (@t-kalinowski, \#2).

### Bug fixes

- [`pif_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  [`paf_total()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  [`pif_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md),
  [`paf_ensemble()`](https://rodrigozepeda.github.io/deltapif/reference/totalpifpaf.md)
  (and hence
  [`weighted_adjusted_pif()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md)
  and
  [`weighted_adjusted_paf()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md))
  used the inverse of the link instead of its derivative, producing
  confidence intervals that were too narrow. Point estimates and
  variances were not affected.

- [`weighted_adjusted_pif()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md)
  and
  [`weighted_adjusted_paf()`](https://rodrigozepeda.github.io/deltapif/reference/weighted_adjusted.md)
  now pass `pif_total_link` and `pif_ensemble_link` (and their inverses
  and derivatives) instead of ignoring them.

### Data

- `dementiarisk` gains a `communality` column (%, Table 2 of Lee et
  al.). Its documentation now states that the prevalence columns are
  percentages and that `sdlog` is the standard error of the log relative
  risk.

### Documentation

- The vignettes used the ARIC exposure prevalences of smoking, excessive
  alcohol and physical inactivity (5.9%, 3.8% and 31.1%) as if they were
  communalities; they now use the published communalities from
  `dementiarisk$communality`, with weights `1 - communality` in the same
  order as the fractions.

- Fixes the variance of the number of cases derived from its confidence
  interval and the alcohol example, which passed relative risks instead
  of log relative risks.

- The README and the Introduction vignette now state that Lee et
  al. obtained their intervals by substituting the confidence limits of
  the relative risk, and the vignettes show relative counterfactuals,
  [`variance_components()`](https://rodrigozepeda.github.io/deltapif/reference/variance_components.md)
  and shared identifiers.

## deltapif 0.4.6

- Fixes bug in printing standard deviation of
  [`attributable_cases()`](https://rodrigozepeda.github.io/deltapif/reference/casecalc.md)
  and
  [`averted_cases()`](https://rodrigozepeda.github.io/deltapif/reference/casecalc.md)
  that multiplied by 100.

## deltapif 0.4.4

- Fixes typo that impeded pdf manual creation.

## deltapif 0.4.3

- Initial release of package
- Improved Readme to include the adjusted fraction.
- Separated introduction from examples vignette
