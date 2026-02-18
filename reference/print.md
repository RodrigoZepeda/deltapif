# Print or show a potential impact fraction

Function to print or show a potential impact fraction object

Function to print or show a `covariance_structure_class`

## Arguments

- accuracy:

  The accuracy of the printed value

- x:

  A `covariance_structure_class`

- ...:

  Additional arguments to pass to `print`

## Examples

``` r
my_pif <- pif(p = 0.2, beta = 1.3, var_beta = 0.1)
#> ! Assuming parameters `p` have no variance Use `var_p` to input their link_variances and/or covariance
print(my_pif)
#> 
#> ── Potential Impact Fraction: [deltapif-124126307530705] ──
#> 
#> PIF = 34.805% [95% CI: 12.300% to 51.535%]
#> standard_deviation(pif %) = 9.864

# Change the ammount of digits to show just 1
print(my_pif, accuracy = 0.1)
#> 
#> ── Potential Impact Fraction: [deltapif-124126307530705] ──
#> 
#> PIF = 34.8% [95% CI: 12.3% to 51.5%]
#> standard_deviation(pif %) = 9.9
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

print(covariance_structure(pif_lead_women))
#>            Women lead
#> Women lead .         
print(covariance_structure2(pif_lead_women, pif_lead_men))
#>            Women lead Men lead
#> Women lead .          .       
#> Men lead   .          .       
print(default_weight_covariance_structure2(pif_men, pif_women))
#>                 Men Men lead Men radiation Women Women lead Women radiation
#> Men             2x2 .        .             .     .          .              
#> Men lead        .   .        .             .     .          .              
#> Men radiation   .   .        .             .     .          .              
#> Women           .   .        .             2x2   .          .              
#> Women lead      .   .        .             .     .          .              
#> Women radiation .   .        .             .     .          .              
print(default_parameter_covariance_structure(pif_tot, parameter = "beta"))
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
```
