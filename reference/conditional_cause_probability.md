# Conditional cause-of-failure probability

Returns a closure computing \\P(K=j \| T=t, \theta)\\ for all
components, conditional on a specific failure time t. By Theorem 6 of
the foundational paper, this equals \\h_j(t; \theta) / \sum_l h_l(t;
\theta)\\.

## Usage

``` r
conditional_cause_probability(model, ...)

# S3 method for class 'dfr_series_md'
conditional_cause_probability(model, ...)
```

## Arguments

- model:

  A likelihood model object.

- ...:

  Additional arguments passed to the returned closure.

## Value

A function with signature `function(t, par, ...)` returning an n x m
matrix where column j gives P(K=j \| T=t, theta).

## Methods (by class)

- `conditional_cause_probability(dfr_series_md)`: Method for
  masked-cause DFR series systems using component hazard ratios.
