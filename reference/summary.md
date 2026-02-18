# Summary of a pif object

Gets the potential impact fraction summary

## Arguments

- object:

  A `pif_class` object.

- ...:

  Additional parameters to pass to `summary` (ignored)

## Examples

``` r
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
summary(my_pif)
#>                PIF standard_deviation             ci_low              ci_up 
#>          0.2858350          0.2689564         -0.4940436          0.6586233 
#>         confidence 
#>          0.9500000 
```
