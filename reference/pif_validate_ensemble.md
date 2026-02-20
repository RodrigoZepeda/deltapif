# Prepare a global_ensemble

Helper function that helps validate the inputs on any global ensemble

## Usage

``` r
pif_validate_ensemble(
  pif1,
  ...,
  weights,
  var_weights,
  var_pif_weights,
  conf_level,
  link,
  link_inv,
  link_deriv,
  quiet,
  is_paf = FALSE,
  weights_sum_to_1 = FALSE,
  label
)
```

## Arguments

- pif1:

  A potential impact fraction (class `pif_class`)

- ...:

  The remaining potential impact fractions or (respectively) population
  attributable fractions.

- weights:

  A vector containing the proportion of the population for each of the
  categories (for each of the pifs given).

- var_weights:

  link_covariance structure for the `weights`. Can be `0` (default) if
  the weights are not random, a vector if only the link_variances of the
  weights are available or a link_covariance matrix.

- var_pif_weights:

  covariance matrix with row `i` and column `j` representing the
  covariance between the `i`-th potential impact fraction of the list
  and the `j`-th weight

- conf_level:

  Confidence level for the confidence interval (default 0.95).

- link:

  Link function such that the `pif` confidence intervals stays within
  the expected bounds.

- link_inv:

  (Optional). If `link` is a function then yhe inverse of `link`. For
  example if `link` is `logit` this should be `inv_logit`.

- link_deriv:

  Derivative of the `link` function. The function tries to build it
  automatically from `link` using
  [`Deriv::Deriv()`](https://rdrr.io/pkg/Deriv/man/Deriv.html).

- quiet:

  Whether to show messages.

- is_paf:

  Whether the computed quantity is a population attributable fraction or
  not

- weights_sum_to_1:

  Boolean flag indicating if the weights sum to 1 (normalized weights)
  or if they are not (unnormalized).

- label:

  Character identifier for the impact fraction. This is for

## Value

A list of validated values for use in `pif_total` or `pif_ensemble`.
