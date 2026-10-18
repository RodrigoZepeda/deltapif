# deltapif 0.4.6

* **`dementiacov` is replaced by `dementiacor`.** The former was built from the
pairwise tetrachoric correlations between the dementia risk factors (exposures)
in ARIC (eTable 2 of Lee et al.) multiplied by the standard errors of the log
relative risks, so it was not a covariance between relative-risk estimates.
`dementiacor` contains the published correlations themselves, documented as
correlations between exposures.

* `dementiarisk` gains a `communality` column (%, Table 2 of Lee et al.).
Its documentation now states that the prevalence columns are percentages. The
vignettes previously used the ARIC exposure prevalences of smoking, excessive
alcohol and physical inactivity (5.9%, 3.8% and 31.1%) as if they were
communalities; they now use the published communalities.

* **Shared parameters are now declared explicitly (behaviour change).** `p` and
`beta` in `pif()`/`paf()` can be given as a named list with one element, e.g.
`beta = list("smoking_dementia" = 0.46)` or `p = list("447" = c(0.12, 0.87))`.
Fractions using the same name (identifier) share the estimate and are
correlated; using one identifier with different values or variances is an
error. Parameters given as plain numbers are now treated as independent:
previously, fractions whose `beta` and `var_beta` (or `p` and `var_p`) were
numerically identical were silently assumed to share the estimate, which
missed shared estimates that were rounded differently and invented sharing
between independent estimates that happened to coincide. Code that relied on
this should give the shared parameters a common identifier.

* `pif_total()`, `pif_ensemble()`, their `paf` versions and the weighted adjusted
fractions accept `var_p` and `var_beta`: covariance matrices between different
estimates whose row and column names are the identifiers (e.g. relative risks
from the same meta-analysis).

* `pif()` accepts the counterfactual prevalence `p_cft` as a function of `p`
(e.g. `function(p) 0.85 * p`) for scenarios defined relative to the current
exposure. The uncertainty in `p` is then propagated through the counterfactual
(numerically or with the new argument `p_cft_deriv`). Numeric `p_cft` values
keep being treated as fixed targets.

* New `variance_components()` splits the variance of a fraction into the parts
due to the prevalences (`p`), the relative risks (`beta`) and the weights.

* `weighted_adjusted_paf()` and `weighted_adjusted_pif()` build their confidence
intervals on the log scale by default (new argument `link`), consistent with the
log-scale variance of the adjusted fractions; use `link = "identity"` for the
previous symmetric intervals. `weighted_adjusted_pif()` and `weighted_adjusted_paf()`
now pass `pif_total_link` and `pif_ensemble_link` (and their inverses and
derivatives) instead of ignoring them.

* Fixes the vignettes and README: the variance of the number of cases derived
from its confidence interval, the alcohol example (which passed relative risks
instead of log relative risks), and the ensemble and adjusted examples (which used
communalities instead of `1 - communality` as weights, in the wrong order). Fixes
the description of `sdlog` in `dementiarisk`.

* Fixes bug in `pif_total()`, `paf_total()`, `pif_ensemble()`, `paf_ensemble()`
(and hence `weighted_adjusted_pif()` and `weighted_adjusted_paf()`) where the
inverse of the link was used instead of its derivative, producing confidence
intervals that were too narrow. Point estimates and variances were not affected.

* Fixes bug in printing standard deviation of `attributable_cases()` and
`averted_cases()` that multiplied by 100. 

# deltapif 0.4.4

* Fixes typo that impeded pdf manual creation. 

# deltapif 0.4.3

* Initial release of package
* Improved Readme to include the adjusted fraction. 
* Separated introduction from examples vignette
