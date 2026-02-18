# Covariance structure class

The covariance structure class represents a collection of matrices and
0's such that `cov[[i]][[j]]` contains the covariance of elements of the
i-th potential impact fraction and the j-th potential impact fraction.

## Usage

``` r
covariance_structure_class(cov_list = list())
```

## Arguments

- cov_list:

  Named list such that `list[[i]][[j]]` represents the covariance
  between elements of the i-th potential impact fraction and the j-th
  potential impact fraction.

## See also

[`as_covariance_structure()`](https://rodrigozepeda.github.io/deltapif/reference/as_covstr.md)
to transform matrices to covariance structures and
[`covariance_structures()`](https://rodrigozepeda.github.io/deltapif/reference/covariance_structures.md)
for default covariance structures

## Examples

``` r
#A simple covariance structure
cov <- covariance_structure_class(
  list(
    "pif1" = list("pif1" = 0.21, "pif2" = 0.12, "pif3" = 0.31),
    "pif2" = list("pif1" = 0.12, "pif2" = 0.33, "pif3" = -0.01),
    "pif3" = list("pif1" = 0.31, "pif2" = -0.01, "pif3" = 0.80)
  )
)
cov
#>      pif1 pif2  pif3 
#> pif1 0.21 0.12  0.31 
#> pif2 0.12 0.33  -0.01
#> pif3 0.31 -0.01 0.8  

#Values can be extracted as in matrices
cov[[1]][[3]]
#> [1] 0.31
cov[["pif1"]][["pif1"]]
#> [1] 0.21

#And assignment also works
cov[[1]][[3]] <- 100
cov[["pif1"]][["pif1"]] <- 500
cov
#>      pif1 pif2  pif3 
#> pif1 500  0.12  100  
#> pif2 0.12 0.33  -0.01
#> pif3 0.31 -0.01 0.8  

#Covariance structures are designed to contain the covariance between
#numbers or vectors (or between numbers and vectors). Hence they can
#contain matrices or vectors too:
mat <- matrix(c(0.1, 0.21, 0.47, -0.3), ncol = 2)
vec <- c(0.22, -0.9, 0.01)

cov2 <- covariance_structure_class(
  list(
    "pif1" = list("pif1" = 0.21, "pif2" = mat, "pif3" = 0.31),
    "pif2" = list("pif1" = mat, "pif2" = 0.33, "pif3" = vec),
    "pif3" = list("pif1" = 0.31, "pif2" = vec, "pif3" = 0.80)
  )
)
cov2
#>      pif1 pif2 pif3
#> pif1 0.21 2x2  0.31
#> pif2 2x2  0.33 1x3 
#> pif3 0.31 1x3  0.8 
cov2[["pif1"]][["pif2"]]
#>      [,1]  [,2]
#> [1,] 0.10  0.47
#> [2,] 0.21 -0.30

```
