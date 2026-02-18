# Link parsers

Functions to parse the link (or inverse link) from a word to a function

## Usage

``` r
parse_link(link_name)

parse_inv_link(link_name)

parse_deriv_link(link_name)
```

## Arguments

- link_name:

  The name of the link or a function.

## Value

A function corresponding to the `link_name`.

## Details

The following are valid link names:

- identity:

  The function `f(x) = x`.with inverse `finv(x) = x`

- logit:

  The function `f(x) = ln(x / (1 - x))` with inverse
  `finv(x) = 1 / (1 + exp(-x))`

- log-complement:

  The function `f(x) = ln(1 - x)` with inverse `finv(x) = 1 - exp(x)`

- hawkins:

  The function `f(x) = ln(x + sqrt(x^2 + 1))` with inverse
  `finv(x) = 0.5 * exp(-x) * (exp(2 * x) - 1)`

- exponential:

  The function `f(x) = exp(x)` with inverse `f(x) = ln(x)`

## Note

If a function is supplied to `link_name` the same function is returned

## See also

[linkfuns](https://rodrigozepeda.github.io/deltapif/reference/linkfuns.md)
and
[inv_linkfuns](https://rodrigozepeda.github.io/deltapif/reference/inv_linkfuns.md)
