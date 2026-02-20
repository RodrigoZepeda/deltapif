# Change the link

Change the link function for the potential impact fraction or population
attributable fraction to a different link.

## Usage

``` r
change_link(x, link = "identity", link_inv = NULL, link_deriv = NULL)
```

## Arguments

- x:

  A `pif_class` object

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

## Value

A `pif_class` object with a different `link`.

## Examples

``` r
#A potential impact fraction
pif1 <- pif(p = 0.2, p_cft = 0.1, beta = 1.2, var_p = 0.01,
  var_beta = 0.2)
pif1
#> 
#> ── Potential Impact Fraction: [deltapif-0118194874409683] ──
#> 
#> PIF = 15.848% [95% CI: -19.419% to 40.699%]
#> standard_deviation(pif %) = 15.028

#Now change the pif to logit to control the negatives
pif1_logit <- change_link(pif1, link = "logit")
pif1_logit
#> 
#> ── Potential Impact Fraction: [deltapif-0118194874409683] ──
#> 
#> PIF = 15.848% [95% CI: 2.027% to 63.158%]
#> standard_deviation(pif %) = 15.028
```
