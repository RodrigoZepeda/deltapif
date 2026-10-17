# deltapif 0.4.6

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
