# Get the values of the derivative of link

Obtain the values of the derivative of the link function of a potential
impact fraction (PIF) or a population attributable fraction (PAF) at
`pif`

## Usage

``` r
link_deriv_vals(x)
```

## Arguments

- x:

  A `pif_class` object

## Value

A number indicating the derivative of `link()` evaluated at `pif`.

## Examples

``` r
if (FALSE) { # \dontrun{
#Create a pif object
pif_obj <- pif(p = 0.499, beta = log(3.6), p_cft = 0.499/2, var_p = 0.001,
  var_beta = 0.1, link = "logit", quiet = TRUE)

#Obtain the value of the link at the derivative
link_deriv_vals(pif_obj)

#This is the same as the derivative of logit evaluated at 0.2823
deriv_logit(coef(pif_obj))
} # }

```
