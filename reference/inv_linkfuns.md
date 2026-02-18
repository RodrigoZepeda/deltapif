# Inverses of link functions

A collection of the inverses of the link functions for the potential
impact fraction.

## Usage

``` r
inv_logit(x)

inv_log_complement(x)

inv_hawkins(x)
```

## Arguments

- x:

  A value such that `inv_link(x)` is a potential impact fraction or a
  population attributable fraction.

## Details

The functions programmed are as follows

\$\$ \text{inv\\logit}(\text{PIF}) = \dfrac{ 1 }{ 1 + \exp(-x) }, \$\$

\$\$ \text{inv\\log-complement}(\text{PIF}) = 1 - \exp(x), \$\$

and

\$\$ \text{inv\\Hawkins}(\text{PIF}) = \frac{1}{2} \exp(-x) \cdot
\big(\exp(x) - 1\big) \cdot \big(\exp(x) + 1\big) \$\$

## See also

[linkfuns](https://rodrigozepeda.github.io/deltapif/reference/linkfuns.md)
for the definition of the link functions and
[deriv_linkfuns](https://rodrigozepeda.github.io/deltapif/reference/deriv_linkfuns.md)
for their derivatives.
