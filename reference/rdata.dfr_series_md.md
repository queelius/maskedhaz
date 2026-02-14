# Random data generation for masked-cause DFR series systems

Returns a function that generates random masked series system data from
the model's data-generating process (DGP). Uses
[`sample_components`](https://queelius.github.io/dfr.dist.series/reference/sample_components.html)
for component lifetimes and applies right-censoring and masking
satisfying C1-C2-C3.

## Usage

``` r
# S3 method for class 'dfr_series_md'
rdata(model, ...)
```

## Arguments

- model:

  A
  [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

A function with signature `function(theta, n, tau = Inf, p = 0, ...)`
that returns a data frame with columns for lifetime, observation type,
and candidate sets.
