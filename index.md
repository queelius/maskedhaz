# maskedhaz

Masked-Cause Likelihood Models for DFR Series Systems

**maskedhaz** provides likelihood-based inference for series systems
with masked component cause of failure. A series system fails when any
component fails, but the causing component may be unknown (masked).
Given candidate sets satisfying conditions C1-C2-C3, this package
computes log-likelihood, score, Hessian, and MLE for the component
parameters.

Unlike closed-form implementations, this package works with
**arbitrary** component hazard functions from `flexhaz` — Weibull,
Gompertz, log-logistic, or any custom `dfr_dist`.

## Installation

Install from [r-universe](https://queelius.r-universe.dev):

``` r
install.packages("maskedhaz", repos = "https://queelius.r-universe.dev")
```

## Quick Start

``` r
library(maskedhaz)

# Three-component exponential series system
model <- dfr_series_md(components = list(
    dfr_exponential(0.1),
    dfr_exponential(0.2),
    dfr_exponential(0.3)
))

# Generate masked data with right-censoring
rdata_fn <- rdata(model)
set.seed(42)
df <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 500, tau = 10, p = 0.5)

# Evaluate log-likelihood
ll_fn <- loglik(model)
ll_fn(df, par = c(0.1, 0.2, 0.3))

# Fit via MLE
solver <- fit(model)
result <- solver(df, par = c(0.5, 0.5, 0.5))
coef(result)    # parameter estimates
vcov(result)    # variance-covariance matrix
confint(result) # confidence intervals
```

## Mixed Component Types

The package supports arbitrary `dfr_dist` component distributions:

``` r
# Weibull wear-out + exponential random failure
model <- dfr_series_md(components = list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
))

rdata_fn <- rdata(model)
set.seed(42)
df <- rdata_fn(theta = c(2, 100, 0.05), n = 300, tau = Inf, p = 0.3)

solver <- fit(model)
result <- solver(df, par = c(1.5, 120, 0.03))
coef(result)
```

## Observation Types

The data frame uses an `omega` column to indicate observation type:

| `omega` value | Meaning                                         |
|---------------|-------------------------------------------------|
| `"exact"`     | Failed at time t, cause masked among candidates |
| `"right"`     | Right-censored: survived past time t            |
| `"left"`      | Left-censored: failed before time t             |
| `"interval"`  | Failed in interval (t, t_upper)                 |

Candidate sets are Boolean columns `x1, x2, ..., xm`.

## Key Features

- **General hazards**: Works with any `dfr_dist` component — no
  closed-form assumptions
- **Full censoring support**: Exact, right, left, and interval censoring
- **Masked cause**: Candidate sets satisfying C1-C2-C3 masking
  conditions
- **MLE fitting**:
  [`fit()`](https://generics.r-lib.org/reference/fit.html) returns
  `fisher_mle` with [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html)
- **Data generation**:
  [`rdata()`](https://queelius.github.io/likelihood.model/reference/rdata.html)
  simulates masked series system data
- **Diagnostics**: Conditional and marginal cause-of-failure
  probabilities
- **Cross-validated**: Matches closed-form results from `maskedcauses`

## Ecosystem

maskedhaz builds on:

- [serieshaz](https://github.com/queelius/serieshaz) — Series system
  distributions
- [flexhaz](https://github.com/queelius/flexhaz) — Dynamic failure rate
  distributions
- [likelihood.model](https://github.com/queelius/likelihood.model) —
  Likelihood model interface
- [algebraic.dist](https://github.com/queelius/algebraic.dist) —
  Distribution generics

Cross-validated against:

- [maskedcauses](https://github.com/queelius/maskedcauses) — Closed-form
  exponential/Weibull series likelihoods
