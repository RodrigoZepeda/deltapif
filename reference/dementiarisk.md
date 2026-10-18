# Exposure and Relative Risk data from Lee et Al

Relative risk and exposure data for dementia risk-factors

## Usage

``` r
dementiarisk
```

## Format

A data frame with 12 rows and 12 columns:

- risk_factor:

  The risk factor for dementia

- RR:

  The relative risk of dementia associated to the risk factor

- lower_CI, upper_CI:

  Lower and upper bounds for the 95% confidence interval

- logrr:

  The logarithm of the relative risk `log(RR)`

- sdlog:

  The standard deviation (standard error) of the logarithm of the
  relative risk `sd(log(RR))`

- total:

  The percentage (%) of individuals exposed in the overall population

- hispanic:

  The percentage (%) of hispanic individuals exposed

- asian:

  The percentage (%) of non-hispanic asian individuals exposed

- black:

  The percentage (%) of non-hispanic black individuals exposed

- white:

  The percentage (%) of non-hispanic white individuals exposed

- communality:

  Communality (%) of the risk factor with the other risk factors (Table
  2 of Lee et al.). The communality weight used to combine fractions is
  `1 - communality / 100`.

## References

Lee, Mark, et al. "Variation in population attributable fraction of
dementia associated with potentially modifiable risk factors by race and
ethnicity in the US." JAMA network open 5.7 (2022): e2219672-e2219672.

## See also

[dementiacor](https://rodrigozepeda.github.io/deltapif/reference/dementiacor.md)
for the correlation between the risk factors
