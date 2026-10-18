# Introduction

``` r

library(deltapif)
```

## Introduction

The `deltapif` R package calculates **Potential Impact Fractions (PIF)
and Population Attributable Fractions (PAF) for aggregated data**. It
uses the delta method to derive confidence intervals, providing a robust
approach for quantifying the burden of disease attributable to risk
factors and the potential impact of interventions.

### Core Concepts: PAF and PIF

The **Population Attributable Fraction (PAF)** answers the question:
*“What fraction of disease cases in a population would be prevented if
we completely eliminated a risk factor?”* It represents the maximum
possible reduction achievable and is often interpreted as the burden of
disease attributed to the exposure.

The **Potential Impact Fraction (PIF)** is a more general measure. It
answers: *“What fraction of disease cases would be prevented if we
changed exposure from its current distribution to a specific
counterfactual scenario?”*

PAF is a specific type of PIF where the counterfactual scenario is the
**theoretical minimum risk exposure level (TMREL)**. We remark that the
TMREL is not always zero. For example, the TMREL for sodium intake is a
specific healthy range (e.g., ~1.6g/day), as both too much and too
little sodium are harmful. For other exposures, such as smoking, the
TMREL can indeed be zero exposure.

> **Note**: The statistical methods underlying the packge assume that
> **the relative risk and exposure prevalence estimates are
> independent** (i.e., derived from different studies or populations).

### Key assumption: Independent (summary) data sources.

The `deltapif` package is designed for a specific, common scenario in
public health:

- The estimate of the log-relative risk (`beta`) comes from one source
  (e.g., a published meta-analysis).

- The estimate of the exposure prevalence (`p`) comes from a separate,
  independent source (e.g., a national survey).

The delta method implementation here relies on this independence. If you
have individual-level data for exposure the [`pifpaf`
package](https://github.com/INSP-RH/pifpaf) is more appropriate as it
leverages the individual-level variability. If individual-level exposure
and outcome data is available from the same source the [`graphPAF`
package](https://CRAN.R-project.org/package=graphPAF) is ideal.

## Usage

### Population Attributable Fraction (PAF)

[Lee et al. (2022)](https://doi.org/10.1001/jamanetworkopen.2022.19672)
estimated the fraction of dementia cases attributable to smoking in the
US. They reported:

- A relative risk of 1.59 (95% CI: 1.15, 2.20)

- A smoking prevalence of 8.5%

The point estimate of the PAF can be calculated using [Levin’s
formula](https://doi.org/10.1016/j.gloepi.2021.100062):

``` r

library(deltapif)

paf(p = 0.085, beta = log(1.59), quiet = TRUE)
#> 
#> ── Population Attributable Fraction: [deltapif-140004351672175] ──
#> 
#> PAF = 4.776% [95% CI: 4.776% to 4.776%]
#> standard_deviation(paf %) = 0.000
```

#### Incorporating Uncertainty

To calculate confidence intervals, we need the variance of the
log-relative risk. The variance can be derived from the confidence
interval following the [Cochrane
Handbook](https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-06#section-6-3-2):

``` r

var_log_rr <- ((log(2.20) - log(1.15)) / (2 * 1.96))^2
var_log_rr
#> [1] 0.0273848
```

We then provide the log-relative risk (`log(1.59)`) and its variance to
[`paf()`](https://rodrigozepeda.github.io/deltapif/reference/pifpaf.md),
specifying the `rr_link` as `exp` to convert the coefficient to a
relative risk by exponentiating the log. Since the prevalence variance
was not reported, we assume `var_p = 0`.

``` r

paf_dementia <- paf(
  p         = 0.085, 
  beta      = log(1.59), 
  var_beta  = var_log_rr, 
  var_p     = 0
)
paf_dementia
#> 
#> ── Population Attributable Fraction: [deltapif-025531705484526] ──
#> 
#> PAF = 4.776% [95% CI: 0.717% to 8.669%]
#> standard_deviation(paf %) = 2.028
```

Lee et al. report **PAF = 4.9% (95% CI: 1.3–9.3)**. Their interval was
obtained by plugging the confidence limits of the relative risk into the
formula (the substitution method), which ignores the uncertainty of the
prevalence and, when several fractions are combined, implicitly assumes
their relative risks are perfectly correlated. The delta-method interval
above is built on the log-complement scale and is therefore slightly
different.

### Potential Impact Fraction (PIF)

[Lee et al. (2022)](https://doi.org/10.1001/jamanetworkopen.2022.19672)
also considered a scenario reducing smoking prevalence by 15% (from 8.5%
to 7.225%). The PIF for this intervention is:

``` r

lee_pif <- pif(
  p        = 0.085, 
  p_cft    = 0.085 * (1 - 0.15), # 15% reduction
  beta     = log(1.59), 
  var_beta = var_log_rr, 
  var_p    = 0
)
lee_pif
#> 
#> ── Potential Impact Fraction: [deltapif-243726361121953] ──
#> 
#> PIF = 0.716% [95% CI: 0.118% to 1.311%]
#> standard_deviation(pif %) = 0.304
```

This result is close to the estimate reported by Lee et al. with the
substitution method: **PIF = 0.7% (95% CI: 0.2–1.4)**.

#### Fixed targets and relative reductions

A numeric `p_cft` is treated as a **fixed target** (e.g. “reduce smoking
prevalence to 7.225%”). Scenarios such as Lee et al.’s are instead
**relative** to the current prevalence (“reduce smoking by 15%”), so the
counterfactual inherits the uncertainty of `p`. In that case, pass
`p_cft` as a function of `p`:

``` r

pif_relative <- pif(
  p        = 0.085,
  p_cft    = function(p) 0.85 * p, # 15% reduction of the current prevalence
  beta     = log(1.59),
  var_beta = var_log_rr,
  var_p    = 0.01^2
)
pif_relative
#> 
#> ── Potential Impact Fraction: [deltapif-000557128674616073] ──
#> 
#> PIF = 0.716% [95% CI: 0.098% to 1.331%]
#> standard_deviation(pif %) = 0.315
```

When `var_p = 0` both forms give the same result. When the prevalence is
uncertain, treating a relative reduction as a fixed target overstates
the variance:

``` r

pif_fixed <- pif(p = 0.085, p_cft = 0.85 * 0.085, beta = log(1.59),
                 var_beta = var_log_rr, var_p = 0.01^2)
pif_fixed
#> 
#> ── Potential Impact Fraction: [deltapif-0621552721415214] ──
#> 
#> PIF = 0.716% [95% CI: -0.537% to 1.954%]
#> standard_deviation(pif %) = 0.635
```

The `variance_components` function shows how much of the variance comes
from the uncertainty in the prevalence (`p`) and in the relative risk
(`beta`):

``` r

variance_components(pif_relative)
#>            p         beta      weights 
#> 6.439931e-07 9.253772e-06 0.000000e+00
variance_components(pif_fixed)
#>            p         beta      weights 
#> 3.111408e-05 9.253772e-06 0.000000e+00
```

### Attributable and averted cases

Attributable and averted cases can be calculated with the
`attributable_cases` function. For example [Dhana et al.
(2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10593099/#SD2) estimate
the number of people with Alzheimer’s Disease in New York, USA 426.5
(400.2, 452.7) thousand. This implies a variance of
`((452.7 - 400.2) / (2*qnorm(0.975)))^2 = 179.3752`.

The number of cases (in thousands) that would be averted if we reduced
smoking by 15% assuming the prevalence of smoking is identical to the
rest of the US is given by:

``` r

averted_cases(426.5, lee_pif, variance = 179.3752)
#> 
#> ── Averted cases: [deltapif-243726361121953] ──
#> 
#> Averted cases = 3.055 [95% CI: 0.504 to 5.606]
#> standard_deviation(averted cases) = 1.302
```

Attributable cases can likewise be estimated using the previous `paf`
as:

``` r

attributable_cases(426.5, paf_dementia, variance = 179.3752)
#> 
#> ── Attributable cases: [deltapif-025531705484526] ──
#> 
#> Attributable cases = 20.368 [95% CI: 3.360 to 37.375]
#> standard_deviation(attributable cases) = 8.677
```

### Combining fractions from subpopulations

Multiple fractions can be combined into totals and ensembles. For
example the fraction among men and women can be combined into an overall
fraction by specifying the distribution of the subgroups in the
population:

``` r

paf_men   <- paf(p = 0.41, beta = 0.31, var_p = 0.001,
                 var_beta = 0.14,
                 label = "Men")
paf_women <- paf(p = 0.37, beta = 0.35, var_p = 0.001, 
                 var_beta = 0.16,
                 label = "Women")
```

Assuming the distribution is 51% women and 49% men:

``` r

paf_total(paf_men, paf_women, weights = c(0.49, 0.51))
#> 
#> ── Population Attributable Fraction: [deltapif-16235488834436] ──
#> 
#> PAF = 13.201% [95% CI: -11.745% to 32.577%]
#> standard_deviation(paf %) = 11.187
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 12.968% (sd %: 15.867) --- [Men]
#> • 13.424% (sd %: 15.773) --- [Women]
#> ────────────────────────────────────────────────────────────────────────────────
```

This is equivalent to calculating:

\\ \textrm{PAF}\_{\text{All}} = 0.49 \cdot \text{PAF}\_{\text{Men}} +
0.51 \cdot \text{PAF}\_{\text{Women}} \\

#### Fractions that share an estimate

In the previous example men and women had different relative risks, so
their fractions are independent. Often, however, subgroups use the
**same** estimate, for example a single relative risk from a
meta-analysis, which makes their fractions correlated. To declare it,
give the parameter as a named list: the name identifies the estimate and
fractions with the same name share it.

``` r

paf_men_shared   <- paf(p = 0.41, beta = list("rr_meta" = 0.31), var_p = 0.001,
                        var_beta = 0.14, label = "Men")
paf_women_shared <- paf(p = 0.37, beta = list("rr_meta" = 0.31), var_p = 0.001,
                        var_beta = 0.14, label = "Women")

paf_total(paf_men_shared, paf_women_shared, weights = c(0.49, 0.51))
#> 
#> ── Population Attributable Fraction: [deltapif-021335574971003] ──
#> 
#> PAF = 12.399% [95% CI: -23.237% to 37.731%]
#> standard_deviation(paf %) = 15.255
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 12.968% (sd %: 15.867) --- [Men]
#> • 11.853% (sd %: 14.693) --- [Women]
#> ────────────────────────────────────────────────────────────────────────────────
```

Without the identifier the same numbers are treated as two independent
estimates, and the interval is too narrow:

``` r

paf_total(
  paf(p = 0.41, beta = 0.31, var_p = 0.001, var_beta = 0.14, label = "Men"),
  paf(p = 0.37, beta = 0.31, var_p = 0.001, var_beta = 0.14, label = "Women"),
  weights = c(0.49, 0.51)
)
#> 
#> ── Population Attributable Fraction: [deltapif-111053489258185] ──
#> 
#> PAF = 12.399% [95% CI: -11.540% to 31.201%]
#> standard_deviation(paf %) = 10.798
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 12.968% (sd %: 15.867) --- [Men]
#> • 11.853% (sd %: 14.693) --- [Women]
#> ────────────────────────────────────────────────────────────────────────────────
```

The same applies to `p` (e.g. `p = list("survey_2020" = c(0.12, 0.30))`
for a categorical exposure). Using one identifier with different values
or variances is an error.

When different estimates are correlated (but not identical) their
covariance can be given to `pif_total`, `pif_ensemble` (and their `paf`
versions) with a matrix whose row and column names are the identifiers:

``` r

cov_beta <- matrix(c(0.14, 0.05,
                     0.05, 0.16), ncol = 2,
                   dimnames = list(c("rr_men", "rr_women"), c("rr_men", "rr_women")))

paf_total(
  paf(p = 0.41, beta = list("rr_men" = 0.31), var_p = 0.001, var_beta = 0.14, label = "Men"),
  paf(p = 0.37, beta = list("rr_women" = 0.35), var_p = 0.001, var_beta = 0.16, label = "Women"),
  weights  = c(0.49, 0.51),
  var_beta = cov_beta
)
#> 
#> ── Population Attributable Fraction: [deltapif-139941046001177] ──
#> 
#> PAF = 13.201% [95% CI: -16.190% to 35.157%]
#> standard_deviation(paf %) = 12.915
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 12.968% (sd %: 15.867) --- [Men]
#> • 13.424% (sd %: 15.773) --- [Women]
#> ────────────────────────────────────────────────────────────────────────────────
```

### Combining fractions from multiple risks

Fractions from disjointed risks can be calculated as an ensemble. For
example the fraction of exposure to lead and the fraction of exposure to
asbestus:

``` r

paf_lead  <- paf(p = 0.41, beta = 0.31, var_p = 0.001,
                 var_beta = 0.014,
                 label = "Lead")
paf_absts <- paf(p = 0.61, beta = 0.15, var_p = 0.001, 
                 var_beta = 0.001,
                 label = "Asbestus")
```

A fraction of **environmental** exposure considering both can be
calculated by multiplying the complements of the weighted fractions.
Here we use (illustrative) weights of `c(0.9, 0.8)`, which correspond to
commonalities of `0.1` and `0.2` (see below):

``` r

paf_ensemble(paf_lead, paf_absts, weights = c(0.9, 0.8))
#> 
#> ── Population Attributable Fraction: [deltapif-0286895218789768] ──
#> 
#> PAF = 18.020% [95% CI: 9.287% to 26.753%]
#> standard_deviation(paf %) = 4.456
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 12.968% (sd %: 5.085) --- [Lead]
#> • 8.985% (sd %: 1.904) --- [Asbestus]
#> ────────────────────────────────────────────────────────────────────────────────
```

where this quantity estimates:

\\ \textrm{PAF}\_{\text{Ensemble}} = 1 - (1 - 0.9 \cdot
\textrm{PAF}\_{\text{Lead}}) \cdot (1 - 0.8 \cdot
\textrm{PAF}\_{\text{Asbestus}}) \\

### Adjusting fractions for commonality

Adjuting for commonality is usually performed when different risks can
be concurrent. In the previous example, exposure to lead and to asbestus
can happen at the same time. [Mukadam et al.
(2019)](https://doi.org/10.1016/S2214-109X(19)30074-9) propose the
individual weighted (adjusted) fractions based on commonality weights.
These weights represent the proportion of the variance shared among risk
factors. To calculate the adjusted fractions one needs to estimate:

\\ \textrm{PIF}\_k^{\text{Adjusted}} = \dfrac{\text{PIF}\_k}{\sum_k
\text{PIF}\_k} \cdot \text{PIF}\_{\text{Overall}} \\ where

\\ \textrm{PIF}^{\text{Overall}} = 1 - \prod\limits_k (1 - w_k
\text{PIF}\_k) \\ with

\\ w_k = 1 - \text{commonality}\_k \\

The adjusted fractions can be calculated with the `weighted_adjusted`
as:

``` r

commonality <- c(0.1, 0.2)
weighted_adjusted_paf(paf_lead, paf_absts, weights = 1 - commonality)
#> $Lead
#> 
#> ── Population Attributable Fraction: [Lead_adj] ──
#> 
#> PAF = 10.645% [95% CI: 4.895% to 23.150%]
#> standard_deviation(paf %) = 4.219
#> 
#> $Asbestus
#> 
#> ── Population Attributable Fraction: [Asbestus_adj] ──
#> 
#> PAF = 7.375% [95% CI: 4.850% to 11.215%]
#> standard_deviation(paf %) = 1.577
```

which returns a named list of the adjusted fractions.

## Additional information

Read the
[examples](https://rodrigozepeda.github.io/deltapif/articles/Examples.html)
vignette

## References

Dhana, Klodian, Todd Beck, Pankaja Desai, Robert S Wilson, Denis A
Evans, and Kumar B Rajan. 2023. “Prevalence of Alzheimer’s Dementia in
the 50 US States and 3142 Counties: A Population Estimate Using the 2020
Bridged-Race Postcensal from the National Center for Health Statistics.”
*Alzheimer’s & Dementia* 19: e074430.

Lee, Mark, Eric Whitsel, Christy Avery, et al. 2022. “Variation in
Population Attributable Fraction of Dementia Associated with Potentially
Modifiable Risk Factors by Race and Ethnicity in the US.” *JAMA Network
Open* 5 (7): e2219672–72.

Mukadam, Naaheed, Andrew Sommerlad, Jonathan Huntley, and Gill
Livingston. 2019. “Population Attributable Fractions for Risk Factors
for Dementia in Low-Income and Middle-Income Countries: An Analysis
Using Cross-Sectional Survey Data.” *The Lancet Global Health* 7 (5):
e596–603.
