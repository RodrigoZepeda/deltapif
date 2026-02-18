# Get the label of a PIF or PAF

Gets the label of a potential impact fraction or a population
attributable fraction

## Arguments

- x:

  A `pif_class` object.

- ...:

  Additional parameters (ignored)

## Examples

``` r
#A simple example
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
                var_beta = 0.2, label = "Test")
names(my_pif)
#> [1] "Test"

#A pif composed of others
my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
                var_beta = 0.2, label = "Test 1")
my_pif2 <- pif(p = 0.4, p_cft = 0.1, beta = 1.3, var_p = 0.1,
                var_beta = 0.2, label = "Test 2")
pif_tot <- pif_total(my_pif1, my_pif2, weights = c(0.2, 0.8),
                label = "Parent")
names(pif_tot)
#> [1] "Test 1" "Test 2"
```
