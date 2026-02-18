# Convert a `covariance_structure` to `matrix`

Transforms a `covariance_structure` into a `matrix`.

## Arguments

- x:

  A `covariance_structure`

- ...:

  Additional parameters (ignored)

## Examples

``` r
as.matrix(covariance_structure_class(list(b = list(a = 1:3))))
#>      [,1]
#> [1,]    1
#> [2,]    2
#> [3,]    3

```
