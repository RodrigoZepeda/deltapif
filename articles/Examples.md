# Examples

``` r
library(deltapif)
```

## Examples

### Example 1: Estimating the PAF for the proportion of dementia attributable to smoking

The article by [M. Lee et al.
(2022)](https://doi.org/10.1001/jamanetworkopen.2022.19672) estimates
the PAF of dementia associated with 12 risk factors (`logrr`) in US
adults. The package includes this data as well as the proportion of
exposed individuals in the respective race column (Hispanic, Asian,
Black and White):

``` r
data(dementiarisk)
```

    #>            risk_factor     logrr       sdlog total hispanic asian black white
    #> 1       Less education 0.4637340 0.119140710  10.7     27.1   6.4  10.6   5.5
    #> 2         Hearing loss 0.6626880 0.174038430  10.8     13.1   6.9   6.5  10.6
    #> 3                  TBI 0.6097656 0.090990178  17.1     10.3   6.0   9.2  20.1
    #> 4         Hypertension 0.4762342 0.167874478  42.2     38.5  38.5  61.0  39.8
    #> 5    Excessive alcohol 0.1655144 0.054020949   3.6      2.0   0.7   2.7   4.2
    #> 6              Obesity 0.4700036 0.091750556  44.0     48.3  14.6  54.3  43.5
    #> 7              Smoking 0.4637340 0.165486566   8.5      6.9   4.9  11.7   8.4
    #> 8           Depression 0.6418539 0.103984905   7.4     10.7   4.3   6.6   7.2
    #> 9     Social isolation 0.4510756 0.086112272  11.9     24.0   8.0  12.1  10.8
    #> 10 Physical inactivity 0.3293037 0.092961816  62.8     68.6  56.6  73.2  61.3
    #> 11            Diabetes 0.4317824 0.075776055  28.6     41.0  44.1  37.2  25.4
    #> 12       Air pollution 0.0861777 0.009362766  22.8     44.4  55.2  41.3  17.2

In this example, we show how to calculate the population attributable
fraction for **Smoking** among the 4 race groups and then how to
calculate the **Total Population Attributable Fraction** of the overall
population.

Notice that **Smoking** has a log relative risk of:

- `beta = 0.4637340`

with variance of

- `var_beta = 0.0273858` (`= 0.165486566^2`).

Among hispanic individuals 6.9% (`p = 0.069`) smoke. Hence their PAF is:

``` r
paf_hispanic <- paf(p = 0.069, beta = 0.4637340, var_beta = 0.0273858, var_p = 0,
                    rr_link = exp, label = "Hispanic")
paf_hispanic
#> 
#> ── Population Attributable Fraction: [Hispanic] ──
#> 
#> PAF = 3.912% [95% CI: 0.569% to 7.142%]
#> standard_deviation(paf %) = 1.676
```

The other fractions can be computed in the same way noting that 4.9% of
non-hispanic asians, 11.7% of non-hispanic blacks and 8.4% of
non-hispanic whites smoke:

``` r
paf_asian <- paf(p = 0.049, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "Asian")
paf_black <- paf(p = 0.117, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "Black")
paf_white <- paf(p = 0.084, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "White")
```

### Example 2: Combining subpopulations into a total

We can calculate the **total population attributable fraction** by
combining the fractions of the subpopulations weighted by the proportion
of the population. According to Wikipedia the distribution in 2020 of
the US population by race was as follows:

``` r
weights <- c("white" = 0.5784, "hispanic" = 0.1873, 
             "black" = 0.1205, "asian" = 0.0592)

#Normalized weights to sum to 1
weights <- weights / sum(weights)
```

The `paf_total` function computes the aggregated `paf` of the whole
population:

``` r
paf_population <- paf_total(paf_white, paf_hispanic, paf_black, 
                            paf_asian, weights = weights, var_weights = 0,
                            label = "All")
paf_population
#> 
#> ── Population Attributable Fraction: [All] ──
#> 
#> PAF = 4.663% [95% CI: 4.486% to 4.839%]
#> standard_deviation(paf %) = 1.979
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 4.722% (sd %: 2.006) --- [White]
#> • 3.912% (sd %: 1.676) --- [Hispanic]
#> • 6.457% (sd %: 2.694) --- [Black]
#> • 2.810% (sd %: 1.218) --- [Asian]
#> ────────────────────────────────────────────────────────────────────────────────
```

where we set `var_weights = 0` as no covariance matrix was known for the
race distribution data (it does exist its just not on Wikipedia).

The `paf` total implements the following formula:

\\ \textrm{PAF}\_{\text{Total}} = w_1 \cdot \text{PAF}\_{\text{White}} +
w_2 \cdot \text{PAF}\_{\text{Hispanic}} + w_3 \cdot
\text{PAF}\_{\text{Black}} + w_4 \cdot \text{PAF}\_{\text{Asian}} \\

where the weights stand for the proportion of the population in each
subgroup.

The `coef`, `confint`, `summary`, and `as.data.frame` functions have
been extended to allow the user to extract information from a `paf` (or
`pif`):

``` r
as.data.frame(paf_population)
#>        value standard_deviation     ci_low      ci_up confidence type label
#> 1 0.04662895         0.01979259 0.04486196 0.04839268       0.95  PAF   All
```

Attributable cases can be calculated for example using the estimates
from [Dhana et al.
(2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10593099/#SD2). They
estimate the number of people with Alzheimer’s Disease in New York, USA
426.5 (400.2, 452.7) thousand. This implies a variance of
`((452.7 - 400.2) / 2*qnorm(0.975))^2 = 2647.005`.

``` r
attributable_cases(426.5, paf_population, variance = 2647)
#> 
#> ── Attributable cases: [All] ──
#> 
#> Attributable cases = 19.887 [95% CI: 2.572 to 37.203]
#> standard_deviation(attributable cases) = 883.469
```

### Example 3: Estimating the PIF for the proportion of dementia reduced by a 15% reduction in smoking.

[M. Lee et al.
(2022)](https://doi.org/10.1001/jamanetworkopen.2022.19672) also
estimates the proportion of dementia cases that would be reduced by a
15% decrease in smoking prevalence. This 15% reduction corresponds to
multiplying the current smoking prevalence by `1 - 0.15 = 0.85`
specified with the prevalence under the counterfactual (`p_cft`)
variable:

``` r
pif_hispanic <- pif(p = 0.069, p_cft = 0.069*0.85, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "Hispanic")
pif_asian    <- pif(p = 0.049, p_cft = 0.049*0.85, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "Asian")
pif_black    <- pif(p = 0.117, p_cft = 0.117*0.85, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "Black")
pif_white    <- pif(p = 0.084, p_cft = 0.084*0.85, beta = 0.4637340, var_beta = 0.0273858, 
                    var_p = 0, rr_link = exp, label = "White")
```

The total fraction of the overall population can be aggregated with
`pif_total`:

``` r
pif_population <- pif_total(pif_white, pif_hispanic, pif_black, 
                            pif_asian, weights = weights, var_weights = 0,
                            label = "All")
pif_population
#> 
#> ── Potential Impact Fraction: [All] ──
#> 
#> PIF = 0.699% [95% CI: 0.695% to 0.703%]
#> standard_deviation(pif %) = 0.297
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 0.708% (sd %: 0.301) --- [White]
#> • 0.587% (sd %: 0.251) --- [Hispanic]
#> • 0.969% (sd %: 0.404) --- [Black]
#> • 0.421% (sd %: 0.183) --- [Asian]
#> ────────────────────────────────────────────────────────────────────────────────
```

We can use the `as.data.frame` function to create a `data.frame` with
all of the results:

``` r
df <- as.data.frame(pif_population, pif_white, pif_black, pif_hispanic, pif_asian)
df
#>         value standard_deviation       ci_low       ci_up confidence type
#> 1 0.006994343        0.002968888 0.0069537857 0.007034899       0.95  PIF
#> 2 0.007082968        0.003009649 0.0011666075 0.012964284       0.95  PIF
#> 3 0.009685883        0.004040705 0.0017344947 0.017573937       0.95  PIF
#> 4 0.005867629        0.002514437 0.0009271875 0.010783639       0.95  PIF
#> 5 0.004214654        0.001826806 0.0006277357 0.007788699       0.95  PIF
#>      label
#> 1      All
#> 2    White
#> 3    Black
#> 4 Hispanic
#> 5    Asian
```

Averted cases can be calculated with the same number as before with the
`averted_cases` function:

``` r
averted_cases(426.5, pif_population, variance = 2647)
#> 
#> ── Averted cases: [All] ──
#> 
#> Averted cases = 2.983 [95% CI: 0.386 to 5.580]
#> standard_deviation(averted cases) = 132.520
```

### Example 4: A categorical population attributable fraction

[Pandeya et al. (2015)](https://doi.org/10.1111/1753-6405.12456)
estimate the population attributable fraction of different cancers given
distinct levels of alcohol consumption for the Australian population.

The `alcohol` dataframe included in the package contains an estimate of
the intake of alcohol by sex as well as the proportion of individuals in
each of the intake-categories:

``` r
data(alcohol)
```

    #>     sex alcohol_g median_intake age_18_plus
    #> 1  MALE      None           0.0        47.2
    #> 2  MALE   >0-<1 g           0.8         1.4
    #> 3  MALE     1-<2g           1.6         2.2
    #> 4  MALE     2-<5g           3.2         8.1
    #> 5  MALE    5-<10g           7.1         8.1
    #> 6  MALE   10-<15g          12.6         6.9
    #> 7  MALE   15-<20g          17.4         4.6
    #> 8  MALE   20-<30g          24.5         7.3
    #> 9  MALE   30-<40g          34.7         4.2
    #> 10 MALE   40-<50g          44.2         3.2
    #> 11 MALE   50-<60g          54.4         2.0
    #> 12 MALE      ≥60g          85.2         4.8

The relative risk of alcohol varies by consumption level (`alcohol_g`).
For example in the case of rectal cancer the estimated relative risk is
1.10 (1.07–1.12) per 10g/day. We will need to compute the relative risk
for each of the consumption level groups using their median intake as
follows:

``` r
#Get alcohol intake for men divided by 10 cause RR is per 10g/day
alcohol_men_intake <- c(0, 0.8, 1.6, 3.2, 7.1, 12.6, 17.4, 24.5, 34.7, 
                        44.2, 54.4, 85.2) / 10 

#Divided by 10 cause risk is per 10 g/day
relative_risk <- exp(log(1.10)*alcohol_men_intake)

#Divide by 100 to get the prevalence in decimals
prevalence    <- c(0.472, 0.014, 0.022, 0.081, 0.081, 0.069, 0.046, 
                   0.073, 0.042, 0.032, 0.02, 0.048)

#Calculate the paf
paf(prevalence, beta = relative_risk, var_p = 0, var_b = 0)
#> 
#> ── Population Attributable Fraction: [deltapif-0102908670119104] ──
#> 
#> PAF = 70.114% [95% CI: 70.114% to 70.114%]
#> standard_deviation(paf %) = 0.000
```

If we want to add uncertainty we can approximate the covariance of the
relative risk as follows:

``` r
#Calculate covariance of the betas given by variance*intake_i*intake_j
beta_sd <- (1.12 - 1.07) / (2*qnorm(0.975))
beta_var <- matrix(NA, ncol = length(prevalence), nrow = length(prevalence))
for (k in 1:length(prevalence)){
  for (j in 1:length(prevalence)){
    beta_var[k,j] <- beta_sd^2*(alcohol_men_intake[j])*(alcohol_men_intake[k])
  }
}

#Calculate the paf
paf_males <- paf(prevalence, beta = relative_risk, var_p = 0, var_b = beta_var, 
               label = "Males")
paf_males
#> 
#> ── Population Attributable Fraction: [Males] ──
#> 
#> PAF = 70.114% [95% CI: 68.480% to 71.663%]
#> standard_deviation(paf %) = 0.812
```

We can repeat the same process for females:

``` r
#Get alcohol intake for women divided by 10 cause RR is per 10g/day
alcohol_female_intake <- c(0, 0.8, 1.6, 3.9, 7.1, 12.6, 17.4, 
                           24.5, 34.7, 44.2, 54.4, 78.1) / 10

#Divided by 10 cause risk is per 10 g/day
relative_risk <- exp(log(1.10)*alcohol_female_intake)

#Divide by 100 to get the prevalence in decimals
prevalence    <- c(0.603, 0.02, 0.054, 0.091, 0.081, 0.061, 0.026, 
                   0.034, 0.015, 0.007, 0.003, 0.005)

#Calculate covariance of the betas given by variance*intake_i*intake_j
beta_sd <- (1.12 - 1.07) / (2*qnorm(0.975))
beta_var <- matrix(NA, ncol = length(prevalence), nrow = length(prevalence))
for (k in 1:length(prevalence)){
  for (j in 1:length(prevalence)){
    beta_var[k,j] <- beta_sd^2*(alcohol_female_intake[j])*(alcohol_female_intake[k])
  }
}

#Calculate the paf
paf_females <- paf(prevalence, beta = relative_risk, var_p = 0, var_b = beta_var, 
                   label = "Females")

paf_females
#> 
#> ── Population Attributable Fraction: [Females] ──
#> 
#> PAF = 65.264% [95% CI: 64.737% to 65.784%]
#> standard_deviation(paf %) = 0.267
```

The total fraction can be computed with `paf_total` assuming both males
and females are 50% each of the population:

``` r
paf_total(paf_males, paf_females, weights = c(0.5, 0.5))
#> 
#> ── Population Attributable Fraction: [deltapif-0974069593265456] ──
#> 
#> PAF = 67.689% [95% CI: 67.426% to 67.950%]
#> standard_deviation(paf %) = 0.427
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 70.114% (sd %: 0.812) --- [Males]
#> • 65.264% (sd %: 0.267) --- [Females]
#> ────────────────────────────────────────────────────────────────────────────────
```

### Example 5: Combining fractions of different risks into an ensemble

M. Lee et al. (2022) combines multiple population attributable fractions
from different risk exposures can be combined into a single fraction.
For example, we can create an attributable fraction for **behavioral**
risk factors by joining the fractions for: smoking, physical inactivity,
and excesive alcohol consumption into one fraction that represents the
reduction in cases if all of these exposures were taken to the level of
minimum risk.

For that purpose we calculate each of the fractions separately and then
use the `paf_ensemble` function to join them. We first obtain the
log-relative risk, its standard deviation, and the total proportion of
individuals exposed (here we highlight the part of the table):

``` r
data(dementiarisk)
```

    #>            risk_factor     logrr      sdlog total
    #> 5    Excessive alcohol 0.1655144 0.05402095   3.6
    #> 7              Smoking 0.4637340 0.16548657   8.5
    #> 10 Physical inactivity 0.3293037 0.09296182  62.8

We then calculate the population attributable fraction for each:

``` r
paf_alcohol  <- paf(p = 0.036, beta = 0.1655144, var_beta = 0.05402095^2, 
                    var_p = 0, label = "Alcohol")
paf_smoking  <- paf(p = 0.085, beta = 0.4637340, var_beta = 0.16548657^2, 
                    var_p = 0, label = "Smoking")
paf_physical <- paf(p = 0.628, beta = 0.3293037, var_beta = 0.09296182^2, 
                    var_p = 0, label = "Physical Activity")
```

The `paf_ensemble` function allows to calculate the ensemble `paf` by
estimating the product:

\\ \textrm{PAF}\_{\text{Ensemble}} = 1 - (1 - w_1
\cdot\textrm{PAF}\_{\text{Alcohol}})\cdot (1 - w_2
\cdot\textrm{PAF}\_{\text{Smoking}})\cdot (1 - w_3\cdot
\textrm{PAF}\_{\text{Physical}}) \\

where the weights are oftentimes calculated via a communality correction
or are set to 1 for disjointed risks. In the case of this paper the
communality estimates were:

- **Physical inactivity:** 31.1

- **Smoking:** 5.9

- **Excessive alcohol:** 3.8

hence we can use them as the `weights` argument:

``` r
paf_ensemble(paf_alcohol, paf_physical, paf_smoking, 
             weights = c(0.038, 0.311, 0.059),
             label = "Behavioural factors")
#> 
#> ── Population Attributable Fraction: [Behavioural factors] ──
#> 
#> PAF = 6.406% [95% CI: 6.202% to 6.610%]
#> standard_deviation(paf %) = 1.627
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 0.644% (sd %: 0.227) --- [Alcohol]
#> • 19.674% (sd %: 5.236) --- [Physical Activity]
#> • 4.776% (sd %: 2.028) --- [Smoking]
#> ────────────────────────────────────────────────────────────────────────────────
```

Finally, variance among the weights can be implemented with the
`var_weights` argument providing a variance-covariance matrix:

``` r
#Variance covariance matrix
weight_variance <- matrix(
  c( 1.00, 0.21, -0.02, 
     0.21, 1.00,  0.09,
    -0.02, 0.09,  1.00), byrow = TRUE, ncol = 3
)
colnames(weight_variance) <- c("Alcohol", "Physical inactivity", "Smoking")
rownames(weight_variance) <- c("Alcohol", "Physical inactivity", "Smoking")

paf_ensemble(paf_alcohol, paf_physical, paf_smoking, 
             weights     = c(0.038, 0.311, 0.059),
             var_weights = weight_variance,
             label       = "Behavioural factors")
#> 
#> ── Population Attributable Fraction: [Behavioural factors] ──
#> 
#> PAF = 6.406% [95% CI: 3.807% to 9.005%]
#> standard_deviation(paf %) = 20.699
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 0.644% (sd %: 0.227) --- [Alcohol]
#> • 19.674% (sd %: 5.236) --- [Physical Activity]
#> • 4.776% (sd %: 2.028) --- [Smoking]
#> ────────────────────────────────────────────────────────────────────────────────
```

### Example 6: Adjusted fractions

Adjusting for communality is a way to discount the contribution of
correlated risk factors. M. Lee et al. (2022) estimates the adjusted
fractions by computing the ensemble fraction as above (lee calls it
Overall):

\\ \textrm{PAF}\_{\text{Ensemble}} = 1 - (1 - w_1
\cdot\textrm{PAF}\_{\text{Alcohol}})\cdot (1 - w_2
\cdot\textrm{PAF}\_{\text{Smoking}})\cdot (1 - w_3\cdot
\textrm{PAF}\_{\text{Physical}}) \\

And then normalizing the proportion of each fraction as part of the
ensemble. For example with alcohol:

\\ \textrm{PAF}\_{\text{Alcohol}}^{\text{Adj}} =
\dfrac{\textrm{PAF}\_{\text{Alcohol}}}{\textrm{PAF}\_{\text{Alcohol}} +
\textrm{PAF}\_{\text{Smoking}} + \textrm{PAF}\_{\text{Physical}}} \cdot
\textrm{PAF}\_{\text{Ensemble}} \\ The same process is repeated for
adjusting smoking and physical by just changing the numerator of the
fraction. This can be computed in the package as:

``` r
#Calculate each of the fractions
paf_alcohol  <- paf(p = 0.036, beta = 0.1655144, var_beta = 0.05402095^2, 
                    var_p = 0, label = "Alcohol")
paf_smoking  <- paf(p = 0.085, beta = 0.4637340, var_beta = 0.16548657^2, 
                    var_p = 0, label = "Smoking")
paf_physical <- paf(p = 0.628, beta = 0.3293037, var_beta = 0.09296182^2, 
                    var_p = 0, label = "Physical Activity")

#Get the communality weights
cweights <- c(0.038, 0.311, 0.059)

#Calculate the adjusted versions
weighted_adjusted_paf(paf_alcohol, paf_smoking, paf_physical, weights = cweights)
#> $Alcohol
#> 
#> ── Population Attributable Fraction: [Alcohol_adj] ──
#> 
#> PAF = 0.068% [95% CI: 0.016% to 0.120%]
#> standard_deviation(paf %) = 0.027
#> 
#> $Smoking
#> 
#> ── Population Attributable Fraction: [Smoking_adj] ──
#> 
#> PAF = 0.505% [95% CI: 0.099% to 0.911%]
#> standard_deviation(paf %) = 0.207
#> 
#> $`Physical Activity`
#> 
#> ── Population Attributable Fraction: [Physical Activity_adj] ──
#> 
#> PAF = 2.080% [95% CI: 1.359% to 2.800%]
#> standard_deviation(paf %) = 0.368
```

### Example 7: Ensuring non-negative fractions

S. Lee et al. (2024) show an example with a relative risk of `rr = 7.0`
(log relative risk `1.94591`) and a log-variance of `0.35.` that results
in a negative confidence interval when applying directly the delta
method:

``` r
paf(p = 0.05, beta = 1.94, var_beta = 0.591^2, var_p = 0)
#> 
#> ── Population Attributable Fraction: [deltapif-0173820105906731] ──
#> 
#> PAF = 22.955% [95% CI: -5.100% to 43.520%]
#> standard_deviation(paf %) = 12.206
```

The `pif` and `paf` functions contain the `logit` link option for
guaranteeing positivity in the fraction’s intervals:

``` r
paf(p = 0.05, beta = 1.94, var_beta = 0.591^2, var_p = 0, link = "logit")
#> 
#> ── Population Attributable Fraction: [deltapif-0178812018276085] ──
#> 
#> PAF = 22.955% [95% CI: 7.152% to 53.541%]
#> standard_deviation(paf %) = 12.206
```

We remark that this option is not turned on by default as the
**positivity** of the fraction depends on:

1.  The relative risk being strictly greater than 1.

2.  Epidemiological rationale on why the fraction should be positive.

### Example 8: Computing averted and attributable cases

Attributable and averted cases can be calculated with the
`attributable_cases` (resp. `averted_cases`) function. For example
[Dhana et al.
(2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10593099/#SD2)
estimates the number of people with Alzheimer’s Disease in New York, USA
426.5 (400.2, 452.7) thousand. This implies a variance of
`((452.7 - 400.2) / 2*qnorm(0.975))^2 = 2647.005`.

To estimate cases we need to:

1.  Estimate the corresponding potential impact fraction (averted) or
    population attributable fraction (attributable)
2.  Use the `averted_cases` (resp. `attributable_cases` function)

We start by calcualting the same fraction as before from [M. Lee et al.
(2022)](https://doi.org/10.1001/jamanetworkopen.2022.19672). This
considers the impact on dementia of reducing smoking prevalence by 15%
(from 8.5% to 7.225%). The PIF for this intervention is:

``` r
pif_smoking <- pif(
  p        = 0.085, 
  p_cft    = 0.085 * (1 - 0.15), # 15% reduction
  beta     = log(1.59), 
  var_beta = ((log(2.20) - log(1.15)) / (2 * 1.96))^2, 
  var_p    = 0
)
pif_smoking
#> 
#> ── Potential Impact Fraction: [deltapif-0698429449342544] ──
#> 
#> PIF = 0.716% [95% CI: 0.118% to 1.311%]
#> standard_deviation(pif %) = 0.304
```

``` r
averted_cases(426.5, pif_smoking, variance = 2647.005)
#> 
#> ── Averted cases: [deltapif-0698429449342544] ──
#> 
#> Averted cases = 3.055 [95% CI: 0.394 to 5.716]
#> standard_deviation(averted cases) = 135.779
```

### Example 9: Computing averted cases (strictly positive)

In some cases it might be of interest to guarantee the positivity of the
averted cases. Following the previous example, a higher (hypothetical)
variance might result in a confidence interval that includes negative
values:

``` r
averted_cases(426.5, pif_smoking, variance = 100000)
#> 
#> ── Averted cases: [deltapif-0698429449342544] ──
#> 
#> Averted cases = 3.055 [95% CI: -2.398 to 8.508]
#> standard_deviation(averted cases) = 278.207
```

When there is an epidemiological rationale on the interval being
strictly positive (for example, the relative risk is known to be
strictly greater than 1) we can change the link function to `logit` to
guarantee the interval is always positive:

``` r
averted_cases(426.5, pif_smoking, variance = 100000, link = "log")
#> 
#> ── Averted cases: [deltapif-0698429449342544] ──
#> 
#> Averted cases = 3.055 [95% CI: 0.513 to 18.203]
#> standard_deviation(averted cases) = 278.207
```

### Example 10: Using Hazard Ratios (HR) and Odds Ratios (OR) instead of Relative Risks (RR)

Sometimes Hazard Rates or Odds Ratios are reported in the literature. As
the potential impact fraction and the population attributable fraction
require relative risks, one needs to convert the former. One way is via
the `EValues` package which uses the method described in [Maya B. Mathur
and Tyler J. VanderWeele
(2019)](https://doi.org/10.1080/01621459.2018.1529598).

In the following example we used the Hazard Ratio of shingles
vaccination on dementia by Wu et al:

``` r
library(EValue)

upper_rr <- HR(0.72, rare = FALSE) |> toRR() |> as.numeric()
lower_rr <- HR(0.67, rare = FALSE) |> toRR() |> as.numeric()
rr_value <- HR(0.69, rare = FALSE) |> toRR() |> as.numeric()
```

which results in the risk of 0.7735267 with interval 0.76, 0.8. From
here one can calculate the potential impact fraction of doubling
vaccination coverage from the baseline of `p = 0.438` to
`p_cft = 0.876`:

``` r
pif(
  p     = 0.438,
  beta  = log(rr_value),
  p_cft = 0.876,
  var_beta = ((log(upper_rr) - log(lower_rr)) / (2 * 1.96))^2,
  quiet = TRUE
)
#> 
#> ── Potential Impact Fraction: [deltapif-0960449159030135] ──
#> 
#> PIF = 11.012% [95% CI: 9.969% to 12.043%]
#> standard_deviation(pif %) = 0.529
```

### References

Dhana, Klodian, Todd Beck, Pankaja Desai, Robert S Wilson, Denis A
Evans, and Kumar B Rajan. 2023. “Prevalence of Alzheimer’s Dementia in
the 50 US States and 3142 Counties: A Population Estimate Using the 2020
Bridged-Race Postcensal from the National Center for Health Statistics.”
*Alzheimer’s & Dementia* 19: e074430.

Lee, Mark, Eric Whitsel, Christy Avery, Timothy M Hughes, Michael E
Griswold, Sanaz Sedaghat, Rebecca F Gottesman, Thomas H Mosley, Gerardo
Heiss, and Pamela L Lutsey. 2022. “Variation in Population Attributable
Fraction of Dementia Associated with Potentially Modifiable Risk Factors
by Race and Ethnicity in the US.” *JAMA Network Open* 5 (7):
e2219672–72.

Lee, Sangjun, Sungji Moon, Kyungsik Kim, Soseul Sung, Youjin Hong,
Woojin Lim, and Sue K Park. 2024. “A Comparison of Green, Delta, and
Monte Carlo Methods to Select an Optimal Approach for Calculating the
95% Confidence Interval of the Population-Attributable Fraction:
Guidance for Epidemiological Research.” *Journal of Preventive Medicine
and Public Health* 57 (5): 499.

Maya B. Mathur, and Tyler J. VanderWeele. 2019. “Sensitivity Analysis
for Unmeasured Confounding in Meta-Analyses.” *Journal of the American
Statistical Association\>*.

Pandeya, Nirmala, Louise F Wilson, Penelope M Webb, Rachel E Neale,
Christopher J Bain, and David C Whiteman. 2015. “Cancers in Australia in
2010 Attributable to the Consumption of Alcohol.” *Australian and New
Zealand Journal of Public Health* 39 (5): 408–13.
