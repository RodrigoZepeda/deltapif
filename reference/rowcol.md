# Row and Column names

Retrieve the row and column names of a `covariance_structure`

## Usage

``` r
row_names(x)

col_names(x)
```

## Arguments

- x:

  A `covariance_structure`

## Value

The names of the rows or columns of a `covariance_structure`

## Examples

``` r
#' #A pif composed of others
my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
                var_beta = 0.2, label = "pif 1")
my_pif2 <- pif(p = 0.4, p_cft = 0.1, beta = 1.3, var_p = 0.1,
                var_beta = 0.2, label = "pif 2")
covst <- covariance_structure2(my_pif1, my_pif2)
row_names(covst)
#> [1] "pif 1" "pif 2"
col_names(covst)
#> [1] "pif 1" "pif 2"

```
