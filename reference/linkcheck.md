# Function that checks the link

Takes the link or the relative risk link of a `pif` or `paf` and checks
that it has been correctly specified.

## Usage

``` r
check_links(link, link_deriv, link_inv)

check_rr_links(rr_link, rr_link_deriv)
```

## Arguments

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- rr_link:

  Link function such that the relative risk is given by `rr_link(beta)`.

- rr_link_deriv:

  Derivative of the link function for the relative risk. The function
  tries to build it automatically from `rr_link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

## Value

Invisible. Called for its side effects
