# MLE fitting for masked-cause DFR series systems

Returns a solver function that finds the maximum likelihood estimates
for component parameters given masked series system data.

## Usage

``` r
# S3 method for class 'dfr_series_md'
fit(object, ...)
```

## Arguments

- object:

  A
  [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

A solver function with signature
`function(df, par, method = "Nelder-Mead", ..., control = list())` that
returns a
[`fisher_mle`](https://queelius.github.io/likelihood.model/reference/fisher_mle.html)
object.

## Details

Uses [`optim`](https://rdrr.io/r/stats/optim.html) to maximize the
log-likelihood. The score function (gradient) is provided for
gradient-based methods. The Hessian at the MLE is computed for
variance-covariance estimation.
