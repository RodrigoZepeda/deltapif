# Convert a `covariance_structure` to `matrix`

Transforms a `covariance_structure` into a `matrix`.

## Arguments

- x:

  A `covariance_structure`

- default:

  How to fill the empty values (default = NA)

- ...:

  Additional parameters (ignored)

## Value

A matrix with the flattened covariance among all the involved fractions

## Details

Because each entry of a covariance structure class can be a matrix
containing the covariances of different fractions, the `as.matrix`
command flattens that matrix into a single object.

## Examples

``` r
#Simple covariance structure to matrix
my_cov <- covariance_structure_class(
  list(
    "pif1" = list("pif1" = 0.21, "pif2" = 0.12, "pif3" = 0.31),
    "pif2" = list("pif1" = 0.12, "pif2" = 0.33, "pif3" = -0.01),
    "pif3" = list("pif1" = 0.31, "pif2" = -0.01, "pif3" = 0.80)
  )
)
as.matrix(my_cov)
#>      pif1  pif2  pif3
#> pif1 0.21  0.12  0.31
#> pif2 0.12  0.33 -0.01
#> pif3 0.31 -0.01  0.80

#More complicated example: the covariance matrix is flattened
mat <- matrix(c(0.1, 0.21, 0.47, -0.3), ncol = 2)
vec <- c(0.22, -0.9, 0.01)

cov2 <- covariance_structure_class(
  list(
    "pif1" = list("pif1" = 0.21, "pif2" = mat, "pif3" = 0.31),
    "pif2" = list("pif1" = mat, "pif2" = 0.33, "pif3" = vec),
    "pif3" = list("pif1" = 0.31, "pif2" = vec, "pif3" = 0.80)
  )
)
as.matrix(cov2)
#>      pif1  pif1 pif1 pif2  pif2 pif2 pif3 pif3 pif3
#> pif1 0.21    NA   NA 0.10  0.47   NA 0.31   NA   NA
#> pif1   NA    NA   NA 0.21 -0.30   NA   NA   NA   NA
#> pif1   NA    NA   NA   NA    NA   NA   NA   NA   NA
#> pif2 0.10  0.47   NA 0.33    NA   NA 0.22 -0.9 0.01
#> pif2 0.21 -0.30   NA   NA    NA   NA   NA   NA   NA
#> pif2   NA    NA   NA   NA    NA   NA   NA   NA   NA
#> pif3 0.31    NA   NA 0.22 -0.90 0.01 0.80   NA   NA
#> pif3   NA    NA   NA   NA    NA   NA   NA   NA   NA
#> pif3   NA    NA   NA   NA    NA   NA   NA   NA   NA

```
