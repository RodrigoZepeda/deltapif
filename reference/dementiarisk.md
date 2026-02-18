# Exposure and Relative Risk data from Lee et Al

Relative risk and exposure data for dementia risk-factors

## Usage

``` r
dementiarisk
```

## Format

A data frame with 12 rows and 11 columns:

- risk_factor:

  The risk factor for dementia

- RR:

  The relative risk of dementia associated to the risk factor

- lower_CI, upper_CI:

  Lower and upper bounds for the 95% confidence interval

- logrr:

  The logarithm of the relative risk `log(RR)`

- sdlog:

  The variance of the logarithm of the relative risk `variance(log(RR))`

- total:

  The proportion of individuals exposed in the overall population

- hispanic:

  The proportion of hispanic individuals exposed

- asian:

  The proportion of non-hispanic asian individuals exposed

- black:

  The proportion of non-hispanic black individuals exposed

- white:

  The proportion of non-hispanic white individuals exposed

## References

Lee, Mark, et al. "Variation in population attributable fraction of
dementia associated with potentially modifiable risk factors by race and
ethnicity in the US." JAMA network open 5.7 (2022): e2219672-e2219672.

## See also

[dementiacov](https://rodrigozepeda.github.io/deltapif/reference/dementiacov.md)
for the covariance between the risk factors
