# Generate masked series system data

Creates masked data from pre-generated component lifetimes. Applies
system lifetime calculation (minimum of components), right-censoring at
tau, and candidate set generation satisfying C1-C2-C3.

## Usage

``` r
generate_masked_series_data(
  comp_lifetimes,
  n,
  m,
  tau,
  p,
  default_lifetime,
  default_omega,
  default_candset
)
```

## Arguments

- comp_lifetimes:

  n x m matrix of component lifetimes

- n:

  number of observations

- m:

  number of components

- tau:

  right-censoring time

- p:

  masking probability for non-failed components

- default_lifetime:

  column name for system lifetime

- default_omega:

  column name for observation type

- default_candset:

  column prefix for candidate sets

## Value

data frame with system lifetime, observation type, and candidate sets
