# Compute the mean under the counterfactual prevalence

Calculates the mean average relative risk of a population under the
counterfactual prevalence: \$\$ E\[\text{RR}\] = \sum\limits\_{i=1}^{n}
p_i^{\text{cft}} \cdot \text{RR}\_i \$\$

## Usage

``` r
mu_cft_fun(p_cft, rr)
```

## Arguments

- p_cft:

  Counterfactual prevalence (proportion) of the exposed individuals for
  each of the `N` exposure levels.

- rr:

  The relative risk for each of the exposure levels.

## Value

The mean relative risk under the counterfactual prevalence

## See also

[`mu_obs_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_obs_fun.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#Consider a popualtion with 3 exposure categories with
#relative risks 1.0, 1.1, 1.3 and prevalences 0.6, 0.3, 0.1
#Notice that the reference relative risk (1) is not added
pval  <- c(0.4, 0.0) #Prevalences on counterfactual scenarios
rrval <- c(1.1, 1.3) #Risks
mu_cft_fun(p_cft = pval, rr = rrval) #Average relative risk
} # }
```
