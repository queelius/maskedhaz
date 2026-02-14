# Hessian of log-likelihood for masked-cause DFR series systems

Returns a Hessian function computed via numerical differentiation of the
log-likelihood using
[`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html).

## Usage

``` r
# S3 method for class 'dfr_series_md'
hess_loglik(model, ...)
```

## Arguments

- model:

  A
  [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

A function with signature `function(df, par, ...)` returning the Hessian
matrix.
