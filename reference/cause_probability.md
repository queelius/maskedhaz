# Marginal cause-of-failure probability

Returns a closure computing \\P(K=j \| \theta)\\ for all components,
marginalized over the system failure time T. By Theorem 5, this equals
\\E_T\[P(K=j \| T, \theta)\]\\.

## Usage

``` r
cause_probability(model, ...)

# S3 method for class 'dfr_series_md'
cause_probability(model, ...)
```

## Arguments

- model:

  A likelihood model object.

- ...:

  Additional arguments passed to the returned closure.

## Value

A function with signature `function(par, ...)` returning an m-vector
where element j gives P(K=j \| theta).

## Methods (by class)

- `cause_probability(dfr_series_md)`: Method for masked-cause DFR series
  systems using Monte Carlo integration.
