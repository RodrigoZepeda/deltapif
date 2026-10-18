# Correlation between dementia risk factors from Lee et al

Pairwise tetrachoric correlations between the 12 dementia risk factors
(exposures) in the Atherosclerosis Risk in Communities (ARIC) study,
from eTable 2 of the supplement of Lee et al. These correlations
describe how the exposures co-occur within individuals and were used to
compute the communality of each risk factor (see the `communality`
column of
[dementiarisk](https://rodrigozepeda.github.io/deltapif/reference/dementiarisk.md)).

## Usage

``` r
dementiacor
```

## Format

A symmetric 12 x 12 correlation matrix with row and column names given
by the risk factors (in the same order as `dementiarisk$risk_factor`).

## Note

These are correlations between exposures, not between the estimates of
the relative risks. They should not be used as a covariance for `beta`.

## References

Lee, Mark, et al. "Variation in population attributable fraction of
dementia associated with potentially modifiable risk factors by race and
ethnicity in the US." JAMA network open 5.7 (2022): e2219672-e2219672.

## See also

[dementiarisk](https://rodrigozepeda.github.io/deltapif/reference/dementiarisk.md)
for the relative risks, prevalence estimates and communalities
