# Masked-Cause Likelihood Model for DFR Series Systems

Constructs a likelihood model for series systems with masked component
cause of failure, where components are arbitrary
[`dfr_dist`](https://queelius.github.io/flexhaz/reference/dfr_dist.html)
distributions. Supports exact, right-censored, left-censored, and
interval-censored observations with candidate sets satisfying C1-C2-C3.

## Usage

``` r
dfr_series_md(
  series = NULL,
  components = NULL,
  par = NULL,
  n_par = NULL,
  lifetime = "t",
  lifetime_upper = "t_upper",
  omega = "omega",
  candset = "x"
)
```

## Arguments

- series:

  A
  [`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.html)
  object. Ignored if `components` is provided.

- components:

  A list of
  [`dfr_dist`](https://queelius.github.io/flexhaz/reference/dfr_dist.html)
  objects. If provided, a `dfr_dist_series` is built from these.

- par:

  Optional concatenated parameter vector.

- n_par:

  Optional integer vector of parameter counts per component.

- lifetime:

  Column name for system lifetime (default `"t"`).

- lifetime_upper:

  Column name for interval upper bound (default `"t_upper"`).

- omega:

  Column name for observation type (default `"omega"`).

- candset:

  Column prefix for candidate set indicators (default `"x"`).

## Value

An object of class
`c("dfr_series_md", "series_md", "likelihood_model")`.

## Details

The model computes the masked-cause log-likelihood for series systems
where the system lifetime is the minimum of independent component
lifetimes, and the causing component is partially observed through
candidate sets.

**Observation types** (stored in the `omega` column):

- `"exact"`:

  Failed at time t, cause masked among candidates

- `"right"`:

  Right-censored: survived past time t

- `"left"`:

  Left-censored: failed before time t

- `"interval"`:

  Failed in interval (t, t_upper)

**Masking conditions**:

- C1:

  Failed component is in candidate set with probability 1

- C2:

  Uniform probability for candidate sets given component cause

- C3:

  Masking probabilities independent of system parameters

## See also

[`is_dfr_series_md`](https://queelius.github.io/maskedhaz/reference/is_dfr_series_md.md)
for the type predicate,
[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.html)
for the series distribution,
[`loglik`](https://queelius.github.io/likelihood.model/reference/loglik.html)
for the likelihood interface

## Examples

``` r
# \donttest{
library(flexhaz)
library(serieshaz)

# From components
model <- dfr_series_md(components = list(
    dfr_exponential(0.1),
    dfr_exponential(0.2),
    dfr_exponential(0.3)
))

# From pre-built series
sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
))
model2 <- dfr_series_md(series = sys)
# }
```
