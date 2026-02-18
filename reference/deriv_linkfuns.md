# Derivatives of link functions

A collection of the derivatives of the link functions for the potential
impact fraction.

## Usage

``` r
deriv_logit(pif)

deriv_log_complement(pif)

deriv_hawkins(pif)

deriv_identity(pif)
```

## Arguments

- pif:

  The value of a potential impact fraction or a population attributable
  fraction

## Details

The functions programmed are as follows

\$\$ \text{deriv\\logit}(\text{PIF}) = \dfrac{ 1 }{ \text{PIF} \cdot
(1 - \text{PIF}) }, \$\$

\$\$ \text{deriv\\log-complement}(\text{PIF}) = \dfrac{ 1 }{
\text{PIF} - 1 }, \$\$

and

\$\$ \text{deriv\\Hawkins}(\text{PIF}) = \dfrac{ 1 }{
\sqrt{\text{PIF}^2 + 1} }, \$\$

## See also

[linkfuns](https://rodrigozepeda.github.io/deltapif/reference/linkfuns.md)
for the definition of the link functions and
[inv_linkfuns](https://rodrigozepeda.github.io/deltapif/reference/inv_linkfuns.md)
for the inverses of the link functions.
