# Attributable cases

Calculates the number of attributable cases or the number of cases that
would be averted under a counterfactual scenario for a given fraction
(either `paf` or `pif`).

## Usage

``` r
averted_cases(
  cases,
  pif,
  variance = 0,
  conf_level = 0.95,
  link = "identity",
  link_inv = NULL,
  link_deriv = NULL
)

attributable_cases(
  cases,
  paf,
  variance = 0,
  conf_level = 0.95,
  link = "identity",
  link_inv = NULL,
  link_deriv = NULL
)
```

## Arguments

- cases:

  The overall number of cases in the population.

- pif:

  A potential impact fraction object created by `pif`, `paf`,
  `pif_total`, `pif_ensemble`, `paf_total` or `paf_ensemble`.

- variance:

  The estimated variance for the cases (default = 0).

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- link:

  Link function such that the case confidence intervals stay within the
  expected bounds (either `log` or `identity`).

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- paf:

  A population attributable fraction object created by `paf`,
  `paf_total` or `paf_ensemble`.

## Value

A `cases_class` object with the attributable cases.

## Details

Negative cases are interpreted as cases that would be caused by the
intervention.

## Formulas

The attributable cases are calculated as: \$\$ \text{Attributable cases}
= \textrm{PAF} \times \textrm{Cases} \$\$ and the averted cases are
respectively: \$\$ \text{Averted cases} = \textrm{PIF} \times
\textrm{Cases} \$\$

The variance is estimated using the product-variance formula: \$\$
\textrm{Var}\[\text{Averted cases}\] = \sigma^2\_{\textrm{Cases}} \cdot
\big( \textrm{PIF}\big)^2 + \sigma^2\_{\textrm{PIF}} \cdot \big(
\textrm{Cases} \big)^2 + \sigma^2\_{\textrm{PIF}} \cdot
\sigma^2\_{\textrm{Cases}} \$\$

## See also

[`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md),
[`paf()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)

## Examples

``` r
frac <- paf(p = 0.499, beta = log(3.6), var_p = 0.002, var_beta = FALSE)
attributable_cases(100, paf = frac)
#> 
#> ── Attributable cases: [deltapif-275541757533686] ──
#> 
#> Attributable cases = 56.473 [95% CI: 52.155 to 60.790]
#> standard_deviation(attributable cases) = 220.300

frac <- pif(p = 0.499, beta = log(3.6), p_cft = 0.1, var_p = 0.002, var_beta = FALSE)
averted_cases(100, pif = frac)
#> 
#> ── Averted cases: [deltapif-00465313804441933] ──
#> 
#> Averted cases = 45.155 [95% CI: 39.715 to 50.596]
#> standard_deviation(averted cases) = 277.578

```
