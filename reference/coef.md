# Extract coefficients of a pif object

Gets the potential impact fraction value

## Usage

``` r
## S7 method for class <deltapif::pif_class>
coef(object, ...)

## S7 method for class <deltapif::cases_class>
coef(object, ...)
```

## Arguments

- object:

  A `pif_class` object.

- ...:

  Additional parameters to pass to `coef` (ignored)

## Value

A `double` containing point-estimate of the object.

## Examples

``` r
my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
coef(my_pif)
#> [1] 0.285835
```
