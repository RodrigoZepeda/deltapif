# Subset a `covariance_structure`

Obtains a smaller `covariance_structure` with the entries given by the
`select` option as a vector.

## Arguments

- x:

  A `covariance_structure`

- select:

  A vector of covariate names to keep in the `covariance_structure`.

- cols:

  A vector of covariate column names to keep in the
  `covariance_structure`.

- rows:

  A vector of covariate row names to keep in the `covariance_structure`.

- negate:

  If `TRUE` subsets the variables that have not been specified

- ...:

  Additional parameters (ignored)

## Examples

``` r
pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
                      label = "Women lead")
pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
                      label = "Women radiation")
covstr <- default_parameter_covariance_structure2(pif_lead_women, pif_rad_women, parameter = "beta")
subset(covstr, "Women lead")
#>            Women lead
#> Women lead 0.015     
subset(covstr, "Women lead", negate = TRUE)
#>                 Women radiation
#> Women radiation 0.022          
subset(covstr, c("Women radiation", "Women lead"))
#>                 Women lead Women radiation
#> Women lead      0.015      .              
#> Women radiation .          0.022          
subset(covstr, 2)
#>                 Women radiation
#> Women radiation 0.022          
subset(covstr, 1:2)
#>                 Women lead Women radiation
#> Women lead      0.015      .              
#> Women radiation .          0.022          
```
