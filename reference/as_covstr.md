# Transform into a covariance structure

Transforms either a matrix or a vector into a covariance structure.

## Usage

``` r
as_covariance_structure(x, col_names = NULL, row_names = NULL, ...)
```

## Arguments

- x:

  A matrix or a vector

- col_names:

  Names to assign to the columns

- row_names:

  Names to assign to the rows

- ...:

  Additional arguments (currently ignored)

## Examples

``` r
mat <- matrix(c(1,3,2,4), ncol = 2,
          dimnames = list(list("pif1", "pif2"), list("pif1", "pif2")))
as_covariance_structure(mat)
#>      pif1 pif2
#> pif1 1    2   
#> pif2 3    4   

#Different colnames than dimnames
as_covariance_structure(mat, col_names = c("first", "second"))
#>      first second
#> pif1 1     2     
#> pif2 3     4     

#Also with a number
as_covariance_structure(2, col_names = "col", row_names = "row")
#>     col
#> row 2  

#Or with a vector
as_covariance_structure(seq(0.1, 0.2, length.out = 4),
    row_names = c("r1","r2","r3","r4"), col_names = "col")
#>    col              
#> r1 0.1              
#> r2 0.133333333333333
#> r3 0.166666666666667
#> r4 0.2              

#As well as a data.frame
data_mat <- as.data.frame(mat)
as_covariance_structure(data_mat, row_names = c("r1","r2"))
#>    pif1 pif2
#> r1 1    2   
#> r2 3    4   
```
