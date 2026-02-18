# Default covariance structures

These are multidimensional arrays of lists where the entry
`list[[i]][[j]]` exists represents the covariance between elements of
the i-th potential impact fraction and the j-th potential impact
fraction. This is particularly useful when handling ensembles and
totals.

## Usage

``` r
covariance_structure(pif, is_variance = FALSE, warning = TRUE)

covariance_structure2(pif1, pif2, warning = TRUE)

default_weight_covariance_structure(pif, is_variance = FALSE, warning = TRUE)

default_weight_covariance_structure2(pif1, pif2, warning = TRUE)

default_parameter_covariance_structure(
  pif,
  parameter = "p",
  is_variance = FALSE,
  warning = TRUE
)

default_parameter_covariance_structure2(
  pif1,
  pif2,
  parameter = "p",
  warning = TRUE
)

default_pif_covariance_structure(pif, is_variance = FALSE, warning = TRUE)

default_pif_covariance_structure2(pif1, pif2, warning = TRUE)

default_weight_pif_covariance_structure(
  pif,
  is_variance = FALSE,
  warning = TRUE
)

default_weight_pif_covariance_structure2(pif1, pif2, warning = TRUE)
```

## Arguments

- pif:

  A potential impact fraction

- is_variance:

  Whether the covariance structure corresponds to a `variance` (i.e.
  when `pif1` and `pif2` are identical)

- warning:

  Whether to throw a warning if `pif1` or `pif2` have common labels

- pif1:

  A potential impact fraction to obtain a covariance structure with
  `pif2`.

- pif2:

  A potential impact fraction to obtain a covariance structure with
  `pif1`,

- parameter:

  Either `beta` or `p`. Indicating which parameter we are calculating
  covariance for.

## Value

A nested list of lists with the entry `[[i]][[j]]` representing the
covariance between elements `i` and `j`.

## Note

The `covariance_structure`s ending in `2` are meant to obtain the
default covariance structure between two fractions `pif1` and `pif2`
while the ones that don't end in `2` are meant to obtain the covariance
structure of a fraction with itself.

## Examples

``` r
pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
                      label = "Women lead")
pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
                      label = "Women radiation")
pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women, label = "Women",
                               weights = c(0.8, 0.72),
                               var_weights = matrix(c(0.3, 0.1, 0.1, 0.4), ncol = 2))

pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
                       label = "Men lead")
pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
                        label = "Men radiation")
pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men, label = "Men",
                        weights = c(0.65, 0.68),
                        var_weights = matrix(c(0.1, -0.2, -0.2, 0.5), ncol = 2))
pif_tot        <- pif_total(pif_men, pif_women,
                        weights = c(0.49, 0.51), label = "Population",
                        var_weights = matrix(c(0.22, 0.4, 0.4, 0.8), ncol = 2))

#This is the default constructor of a covariance. Use it for custom covariances
covariance_structure(pif_lead_women)
#>            Women lead
#> Women lead .         
covariance_structure2(pif_lead_women, pif_lead_men)
#>            Women lead Men lead
#> Women lead .          .       
#> Men lead   .          .       
default_weight_covariance_structure2(pif_men, pif_men)
#>               Men Men lead Men radiation
#> Men           2x2 .        .            
#> Men lead      .   .        .            
#> Men radiation .   .        .            
default_weight_covariance_structure(pif_tot)
#>                 Population Men Men lead Men radiation Women Women lead
#> Population      2x2        .   .        .             .     .         
#> Men             .          2x2 .        .             .     .         
#> Men lead        .          .   .        .             .     .         
#> Men radiation   .          .   .        .             .     .         
#> Women           .          .   .        .             2x2   .         
#> Women lead      .          .   .        .             .     .         
#> Women radiation .          .   .        .             .     .         
#>                 Women radiation
#> Population      .              
#> Men             .              
#> Men lead        .              
#> Men radiation   .              
#> Women           .              
#> Women lead      .              
#> Women radiation .              
default_weight_covariance_structure2(pif_men, pif_women)
#>                 Men Men lead Men radiation Women Women lead Women radiation
#> Men             2x2 .        .             .     .          .              
#> Men lead        .   .        .             .     .          .              
#> Men radiation   .   .        .             .     .          .              
#> Women           .   .        .             2x2   .          .              
#> Women lead      .   .        .             .     .          .              
#> Women radiation .   .        .             .     .          .              
default_parameter_covariance_structure(pif_tot, parameter = "beta")
#>                 Population Men Men lead Men radiation Women Women lead
#> Population      .          .   .        .             .     .         
#> Men             .          .   .        .             .     .         
#> Men lead        .          .   0.015    .             .     0.015     
#> Men radiation   .          .   .        0.022         .     .         
#> Women           .          .   .        .             .     .         
#> Women lead      .          .   0.015    .             .     0.015     
#> Women radiation .          .   .        0.022         .     .         
#>                 Women radiation
#> Population      .              
#> Men             .              
#> Men lead        .              
#> Men radiation   0.022          
#> Women           .              
#> Women lead      .              
#> Women radiation 0.022          
default_parameter_covariance_structure2(pif_lead_women, pif_lead_men, parameter = "beta")
#>            Women lead Men lead
#> Women lead 0.015      0.015   
#> Men lead   0.015      0.015   
```
