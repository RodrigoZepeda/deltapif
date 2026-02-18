# Extract coefficients of a pif object

Gets the potential impact fraction value

## Arguments

- object:

  A `pif_class` object.

- ...:

  Additional parameters to pass to `coef` (ignored)

## Examples

``` r
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
coef(my_pif)
#> [1] 0.285835
```
