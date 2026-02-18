# Convert to list

Converts the object into a list.

## Arguments

- x:

  Either a `covariance_structure`, a `pif_atomic_class` or
  `pif_global_ensemble_class`

- ...:

  Additional parameters (ignored)

## Value

The object as a list

## Examples

``` r
#FOR POTENTIAL IMPACT FRACTIONS
#------------------------------------------------------------------------
#Potential impact fraction for women
paf_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001,
    label = "Lead women")
paf_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001,
    label = "Radiation women")
paf_women      <- paf_ensemble(paf_lead_women, paf_rad_women,
    label = "Women")

as.list(paf_women)
#> $Women
#> 
#> ── Population Attributable Fraction: [Women] ──
#> 
#> PAF = 75.299% [95% CI: 71.706% to 78.892%]
#> standard_deviation(paf %) = 2.435
#> ────────────────────────────────── Components: ─────────────────────────────────
#> • 68.422% (sd %: 2.531) --- [Lead women]
#> • 21.778% (sd %: 4.489) --- [Radiation women]
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> $`Lead women`
#> 
#> ── Population Attributable Fraction: [Lead women] ──
#> 
#> PAF = 68.422% [95% CI: 63.051% to 73.012%]
#> standard_deviation(paf %) = 2.531
#> 
#> $`Radiation women`
#> 
#> ── Population Attributable Fraction: [Radiation women] ──
#> 
#> PAF = 21.778% [95% CI: 12.466% to 30.100%]
#> standard_deviation(paf %) = 4.489
#> 

#FOR COVARIANCE STRUCTURES
#------------------------------------------------------------------------
cov_str <- default_parameter_covariance_structure(paf_women)
as.list(cov_str)
#> $Women
#> $Women$Women
#> [1] 0
#> 
#> $Women$`Lead women`
#> [1] 0
#> 
#> $Women$`Radiation women`
#> [1] 0
#> 
#> 
#> $`Lead women`
#> $`Lead women`$Women
#> [1] 0
#> 
#> $`Lead women`$`Lead women`
#> [1] 0.001
#> 
#> $`Lead women`$`Radiation women`
#> [1] 0
#> 
#> 
#> $`Radiation women`
#> $`Radiation women`$Women
#> [1] 0
#> 
#> $`Radiation women`$`Lead women`
#> [1] 0
#> 
#> $`Radiation women`$`Radiation women`
#> [1] 0.001
#> 
#> 
```
