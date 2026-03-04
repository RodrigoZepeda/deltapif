# Extract confidence intervals

Gets the confidence interval for any fraction or cases.

## Arguments

- object:

  A `pif_class` or a `cases_class`.

- level:

  Level of confidence desired.

- ...:

  Additional parameters to pass to `confint` (ignored)

## Value

A vector containing the lower and upper bounds of the confidence
interval for the object.

## Examples

``` r
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#Default 95% CI
confint(my_pif)
#> [1] -0.4940436  0.6586233

#Custom 90% ci:
confint(my_pif, level = 0.90)
#> [1] -0.3268596  0.6156099
```
