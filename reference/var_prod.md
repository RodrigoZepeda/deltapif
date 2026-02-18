# Variance product formula

Implements the variance product formula for the product of two random
variables. This is the formula for `var(X*Y)` given by: \$\$
\text{Var}\[XY\] = Y^2\cdot\text{Var}\[X\] + X^2\cdot \text{Var}\[Y\] +
\text{Var}\[X\] \cdot \text{Var}\[Y\] \$\$

## Usage

``` r
var_prod(x, y, var_x, var_y)
```

## Arguments

- x:

  Numeric value representing the mean of the first random variable X.

- y:

  Numeric value representing the mean of the second random variable Y.

- var_x:

  Numeric value representing the variance of X. Must be non-negative.

- var_y:

  Numeric value representing the variance of Y. Must be non-negative.

## Value

Numeric value representing the variance of the product X\*Y.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
var_prod(x = 2, y = 3, var_x = 0.5, var_y = 0.8)

# When one variable has zero variance
var_prod(x = 5, y = 4, var_x = 0, var_y = 1.2)

# When both variables have zero variance
var_prod(x = 10, y = 20, var_x = 0, var_y = 0)

# With larger variances
var_prod(x = 100, y = 200, var_x = 25, var_y = 36)
} # }
```
