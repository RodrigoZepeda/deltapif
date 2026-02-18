# Compute the potential impact fraction

Calculates the potential impact fraction following the formula: \$\$
\text{PIF} = \dfrac{ \sum\limits\_{i=1}^n p_i \cdot \text{RR}\_i -
\sum\limits\_{i=1}^n p_i^{\text{cft}} \cdot \text{RR}\_i }{
\sum\limits\_{i=1}^n p_i \cdot \text{RR}\_i } \$\$

## Usage

``` r
pif_fun(p, p_cft, rr)
```

## Arguments

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

- rr:

  The relative risk for each of the exposure levels.

## Value

A potential impact fraction (numeric)

## See also

[`mu_cft_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_cft_fun.md),
[`mu_obs_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_obs_fun.md),
and
[`pif_fun2()`](https://rodrigozepeda.github.io/deltapif/reference/pif_fun2.md)
to calculate from the average relative risks. For the main function in
the package see
[`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#This example comes from Levin 1953
#Relative risk of lung cancer given smoking was 3.6
#Proportion of individuals smoking where 49.9%
#Counterfactual is 0% smoking
pif_fun(0.499, 0.0, 3.6)

#You can also use multiple exposure categories. As an example
#These are the relative risks for age-group 25-59 for each BMI level:
rr <- c("<18.5" = 1.38, "25 to <30" = 0.83,
  "30 to <35" = 1.20, ">=35" = 1.83)

#While the prevalences are:
p <- c("<18.5" = 1.9, "25 to <30" = 34.8,
  "30 to <35" = 17.3, ">=35" = 13.3) / 100

#A counterfactual of reducing the prevalence of >=35 by half and
#having them be on the "30 to <35" category instead
p_cft <- c("<18.5" = 1.9, "25 to <30" = 34.8,
  "30 to <35" = 17.3 + 13.3/2, ">=35" = 13.3/2) / 100

pif_fun(p = p, p_cft = p_cft, rr = rr)
} # }
```
