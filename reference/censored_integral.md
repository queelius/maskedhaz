# Compute integral for left/interval censored contributions

Integrates `[sum_{j in C_i} h_j(t)] * S_sys(t)` over `(lower, upper)`.

## Usage

``` r
censored_integral(h_fns, H_fns, layout, par, C_i, lower, upper, m)
```

## Arguments

- h_fns:

  list of hazard closures

- H_fns:

  list of cumulative hazard closures

- layout:

  parameter layout

- par:

  parameter vector

- C_i:

  logical vector of candidate components

- lower:

  lower integration bound

- upper:

  upper integration bound

- m:

  number of components

## Value

numeric integral value
