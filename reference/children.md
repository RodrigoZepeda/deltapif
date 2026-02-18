# Get the children labels from a `pif_class`

Gets the labels of the pif elements that make up a
`pif_global_ensemble_class`.

## Usage

``` r
children(x, ...)
```

## Arguments

- x:

  A `pif_global_ensemble_class`

- ...:

  Additional arguments (currently ignored)

## Value

A character vector with the names of the fractions that make up `x`

## Note

The result is `NULL` if `x` is a `pif_atomic_class`
