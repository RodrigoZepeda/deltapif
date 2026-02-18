# Names of a PIF's components

Returns a character vector of the names of all of the fractions that
make up a `pif` object. This is particularly useful for `pif_total` and
`pif_ensemble` to obtain the names that build them up.

## Usage

``` r
flatten_names(pif)
```

## Arguments

- pif:

  A potential impact fraction of either `pif_atomic_class` or
  `pif_global_ensemble_class`

## Value

A character vector with the names of all the fractions that make up the
`pif`.

## Examples

``` r
paf_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001,
        label = "Women lead")
paf_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001,
        label = "Women radiation")
paf_women      <- paf_ensemble(paf_lead_women, paf_rad_women,
        label = "Women")
paf_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001,
        label = "Men lead")
paf_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001,
        label = "Men radiation")
paf_men        <- paf_ensemble(paf_lead_men, paf_rad_men,
        label = "Men")
paf_tot        <- paf_total(paf_men, paf_women, weights = c(0.49, 0.51),
        label = "Population")

#For a single PIF return the names
flatten_names(paf_lead_women)
#> [1] "Women lead"

#For an ensemble return the ones that make them up
flatten_names(paf_women)
#> [1] "Women"           "Women lead"      "Women radiation"

#For totals return the ones that make them up
flatten_names(paf_tot)
#> [1] "Population"      "Men"             "Men lead"        "Men radiation"  
#> [5] "Women"           "Women lead"      "Women radiation"
```
