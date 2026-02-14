# Score function for masked-cause DFR series systems

Returns a score (gradient) function computed via numerical
differentiation of the log-likelihood using
[`grad`](https://rdrr.io/pkg/numDeriv/man/grad.html).

## Usage

``` r
# S3 method for class 'dfr_series_md'
score(model, ...)
```

## Arguments

- model:

  A
  [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

A function with signature `function(df, par, ...)` returning the
gradient vector.
