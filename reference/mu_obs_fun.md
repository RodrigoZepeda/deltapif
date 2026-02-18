# Compute the mean under the observed prevalence

Calculates the mean average relative risk of a population under the
observed prevalence: \$\$ E\[\text{RR}\] = \sum\limits\_{i=1}^{n} p_i
\cdot \text{RR}\_i \$\$

## Usage

``` r
mu_obs_fun(p, rr)
```

## Arguments

- p:

  Prevalence (proportion) of the exposed individuals for each of the `N`
  exposure levels.

- rr:

  The relative risk for each of the exposure levels.

## Value

The mean relative risk under the observed prevalence

## See also

[`mu_cft_fun()`](https://rodrigozepeda.github.io/deltapif/reference/mu_cft_fun.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#Consider a popualtion with 3 exposure categories with
#relative risks 1.0, 1.1, 1.3 and prevalences 0.6, 0.3, 0.1
#Notice that the reference relative risk (1) is not added
pval  <- c(0.3, 0.1) #Prevalences
rrval <- c(1.1, 1.3) #Risks
mu_obs_fun(p = pval, rr = rrval) #Average relative risk
} # }
```
