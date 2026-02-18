# Confidence interval for a potential impact fraction

Confidence interval for a potential impact fraction

## Usage

``` r
pif_atomic_ci(link_vals, link_variance, conf_level, link_inv)
```

## Arguments

- link_vals:

  Values of the link function evaluated at `pif` (i.e. `link(pif)`)

- link_variance:

  Link_variance estimate for the linked potential impact fraction (i.e.
  for `link(pif)`)

- conf_level:

  Confidence level for the interval

- link_inv:

  Inverse of the link function used to compute `link_vals` and
  `link_variance`.

## Value

A vector with the lower and upper bounds of the confidence interval

## See also

[confint](https://rodrigozepeda.github.io/deltapif/reference/confint.md)
