# Extract weights of a pif_global_ensemble

Gets the weights of a `pif_global_ensemble` object

## Arguments

- object:

  A `pif_global_ensemble` object.

- ...:

  Additional parameters to pass to `weights` (ignored)

## Examples

``` r
my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
my_pif2 <- pif(p = 0.3, p_cft = 0.1, beta = 1.5, var_p = 0.1, var_beta = 0.2)
my_pif  <- pif_total(my_pif1, my_pif2, weights = c(0.8, 0.2))
weights(my_pif)
#> [1] 0.8 0.2
```
