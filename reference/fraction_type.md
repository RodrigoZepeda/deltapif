# Get the type of the fraction

Obtain whether a fraction is a potential impact fraction (PIF) or a
population attributable fraction (PAF)

## Usage

``` r
fraction_type(x)
```

## Arguments

- x:

  A `pif_class` object

## Value

A character either `PIF` or `PAF` depending on the object

## Examples

``` r
#A potential impact fraction
pif1 <- pif(p = 0.2, p_cft = 0.1, beta = 1.2, quiet = TRUE)
fraction_type(pif1)
#> [1] "PIF"

#A population attributable fraction
paf1 <- paf(p = 0.2, beta = 1.2, quiet = TRUE)
fraction_type(paf1)
#> [1] "PAF"
```
