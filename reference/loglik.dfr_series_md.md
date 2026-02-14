# Log-likelihood for masked-cause DFR series systems

Returns a log-likelihood function for a series system with masked
component cause of failure. Supports four observation types: exact
failures, right-censored, left-censored, and interval-censored.

## Usage

``` r
# S3 method for class 'dfr_series_md'
loglik(model, ...)
```

## Arguments

- model:

  A
  [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

A function with signature `function(df, par, ...)` that computes the
log-likelihood.

## Details

Log-likelihood contributions by observation type:

- Exact (\\\omega = \\ "exact"):

  \\\log L_i = \log(\sum\_{j \in C_i} h_j(t_i)) - H\_{sys}(t_i)\\

- Right-censored (\\\omega = \\ "right"):

  \\\log L_i = -H\_{sys}(t_i)\\

- Left-censored (\\\omega = \\ "left"):

  \\\log L_i = \log \int_0^{t_i} \[\sum\_{j \in C_i} h_j(u)\]
  S\_{sys}(u) \\ du\\

- Interval-censored (\\\omega = \\ "interval"):

  \\\log L_i = \log \int\_{t_i}^{t\_{upper,i}} \[\sum\_{j \in C_i}
  h_j(u)\] S\_{sys}(u) \\ du\\

The exact and right-censored paths use direct hazard/cumulative hazard
calls. Left and interval censoring require numerical integration via
[`integrate`](https://rdrr.io/r/stats/integrate.html).
