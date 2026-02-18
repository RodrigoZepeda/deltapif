# Compute the potential impact fraction

Compute the potential impact fraction

## Usage

``` r
pif_fun2(mu_obs, mu_cft)
```

## Arguments

- mu_obs:

  The average value of the relative risk in the observed population.

- mu_cft:

  The average value of the counterfactual relative risk in the
  population.

## Value

A potential impact fraction (numeric)

## See also

[`mu_cft_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_cft_fun.md),
[`mu_obs_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_obs_fun.md),
and
[`pif_fun()`](https://rodrigozepeda.github.io/deltapif/reference/pif_fun.md)
to calculate from the prevalence and relative risks. For the main
function in the package see
[`pif()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#This example comes from Levin 1953
#Relative risk of lung cancer given smoking was 3.6
#Proportion of individuals smoking where 49.9%
#Counterfactual is 0% smoking
mu_obs <- mu_obs_fun(p = 0.499, rr = 3.6)
mu_cft <- mu_cft_fun(p_cft = 0.0, rr = 3.6)
pif_fun2(mu_obs, mu_cft)

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

mu_obs <- mu_obs_fun(p = p, rr = rr)
mu_cft <- mu_cft_fun(p_cft = p_cft, rr = rr)
pif_fun2(mu_obs, mu_cft)
} # }
```
