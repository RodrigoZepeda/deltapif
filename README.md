
<!-- README.md is generated from README.Rmd. Please edit that file -->

# deltapif <img src="man/figures/logo.png" align="right" height="127" alt="The logo of the deltapif method showing an observed and a counterfactual population with coloured squares. A line partitions them in the middle. The image reads 'deltapif'." />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/deltapif)](https://CRAN.R-project.org/package=deltapif)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/RodrigoZepeda/deltapif/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RodrigoZepeda/deltapif/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/RodrigoZepeda/pifes/graph/badge.svg)](https://app.codecov.io/gh/RodrigoZepeda/deltapif)
<!-- badges: end -->

Calculate both **Potential Impact Fractions (PIF)** and **Population
Attributable Fractions (PAF)** for aggregated data and their confidence
intervals using the delta method.

## Installation

You can install the development version of pifes like so:

``` r
remotes::install_github("RodrigoZepeda/deltapif")
```

## Calculating a population attributable fraction

To estimate a population attributable fraction two ingredients are
required:

- `beta`: the relative risk or a value from which to calculate the
  relative risk (with variance `var_beta`).
- `p`: the exposure prevalence (with variance `var_p`).

and in the case of potential impact fractions we also require: -
`p_cft`: the counterfactual prevalence.

> **Note** An important hypothesis of the current method is that the
> relative risk estimate and the prevalence of exposure are both
> independent in the sense that they were estimated in different
> populations and studies.

For example [Lee et
al](https://doi.org/10.1001/jamanetworkopen.2022.19672) estimate the
attributable fraction of smoking on dementia. For that purpose they
report the following:

- A relative risk of `1.59 (1.15, 2.20)`
- An exposure prevalence of `8.5`

We can calculate a point estimate of the PAF as follows:

``` r
library(deltapif)
paf(p = 0.085, beta = 1.59, quiet = TRUE)
#> 
#> ── Population Attributable Fraction ──
#> 
#> PAF = 4.776% [95% CI: 4.776% to 4.776%]
#> standard_deviation(paf %) = 0.000
#> standard_deviation(link(paf)) = 0.000
```

Note that this follows the formula by
[Levin](https://doi.org/10.1016/j.gloepi.2021.100062):

$$
\textrm{PAF} = \frac{p \cdot (\text{RR} - 1)}{1 + p \cdot (\text{RR} - 1)}
$$

Additional examples show how to calculate the PAF for multiple
categories.

### Adding uncertainty

Note that there is no uncertainty in the fraction as we have not inputed
the standard deviation of `beta`. To do so we notice that the variance
for the relative risk is actually for the log of beta. We follow the
formula from [Cochrane’s
handbook](https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm):

$$
\text{variance}_{\ln(\beta)} = \Bigg(\frac{\ln(\text{upper limit}) - \ln(\text{lower limit})}{2\times 1.95}\Bigg)^2
$$

``` r
((log(2.20) - log(1.15)) / (2*1.95))^2
#> [1] 0.02766639
```

And because the uncertainty was specified for the log of the relative
risk we use the natural logarithm `ln` of `1.59` as our beta where
`log(1.59) = 0.4637` and specify that the way to get our relative risk
is through the exponential function. Finally as no uncertainty was given
for `p` we set `var_p = 0`.

``` r
paf(p = 0.085, beta = 0.4637, var_beta = 0.02766639, var_p = 0, rr_link = exp)
#> 
#> ── Population Attributable Fraction ──
#> 
#> PAF = 4.775% [95% CI: 0.695% to 8.688%]
#> standard_deviation(paf %) = 2.038
#> standard_deviation(link(paf)) = 0.021
```

## Calculating a potential impact fraction

The paper also considers the potential impact fraction of a 15%
reduction of smoking this can be achieved with the `pif` function by
specifying the counterfactual distribution (in this case the 15%
reduction results in $0.085 \times (1 - 0.15) = 0.07225$)

``` r
pif(p = 0.085, p_cft = 0.07225, beta = 0.4637, var_beta = 0.02766639, var_p = 0, rr_link = exp)
#> 
#> ── Potential Impact Fraction ──
#> 
#> PIF = 0.716% [95% CI: 0.115% to 1.314%]
#> standard_deviation(pif %) = 0.306
#> standard_deviation(link(pif)) = 0.003
```

Note that both the PIF and PAF result in similar point estimates and CIs
as in *Lee et al* who report:

- **PAF**: 4.9 (1.3-9.3)
- **PIF**: 0.7 (0.2-1.4)

## More

Visit the [deltapif](https://rodrigozepeda.github.io/deltapif/) website
for more examples and additional operations you can do with `PIF`.
