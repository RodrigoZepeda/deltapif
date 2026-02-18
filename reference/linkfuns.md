# Link functions

A collection of common link functions, for calculating the link_variance
of the potential impact fraction.

## Usage

``` r
logit(pif)

log_complement(pif)

hawkins(pif)

identity_link(pif)
```

## Arguments

- pif:

  The value of a potential impact fraction or a population attributable
  fraction

## Details

The functions programmed are as follows

\$\$ \text{logit}(\text{PIF}) = \ln\Bigg(\dfrac{ \text{PIF} }{ 1 -
\text{PIF} }\Bigg), \$\$

\$\$ \text{log-complement}(\text{PIF}) = \ln\big(1 - \text{PIF}\big),
\$\$

and

\$\$ \text{Hawkins}(\text{PIF}) = \ln\Big(\text{PIF} +
\sqrt{\text{PIF}^2 + 1}\Big). \$\$

## Note

When used, the link_variance is calculated for `linkfun(pif)`.

## See also

[inv_linkfuns](https://rodrigozepeda.github.io/deltapif/reference/inv_linkfuns.md)
for their inverses and
[deriv_linkfuns](https://rodrigozepeda.github.io/deltapif/reference/deriv_linkfuns.md)
for their derivatives
