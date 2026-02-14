# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`maskedhaz` provides likelihood-based inference for series systems with masked component cause of failure, using arbitrary dynamic failure rate (DFR) component distributions. Given series system data where the causing component is uncertain (masked), it computes log-likelihood, score, Hessian, and MLE for the component parameters.

The package decouples from closed-form assumptions: any `dfr_dist` component distribution works, with numerical integration for censored observations and `numDeriv` for derivatives.

## Development Commands

```r
devtools::document()      # Update NAMESPACE and .Rd files from roxygen2
devtools::test()          # Run all tests
devtools::test(filter = "loglik")   # Run specific test file
devtools::check()         # Run R CMD check
devtools::install()       # Install locally
covr::package_coverage()  # Check test coverage (target >= 97%)
```

Always run `devtools::document()` after modifying roxygen2 comments. NAMESPACE is auto-generated — never edit manually.

## Dependencies

**Required (Imports):**
- `serieshaz` — Series system distributions (Phase 1)
- `flexhaz` — DFR component distributions and constructors
- `algebraic.dist` — Distribution interface generics
- `likelihood.model` — Likelihood model interface + `fisher_mle` class
- `generics` — Generic function infrastructure (`fit`)
- `numDeriv` — Numerical score and Hessian
- `stats` — `integrate`, `optim`, `runif`

**Optional (Suggests):**
- `maskedcauses` — Cross-validation against closed-form exponential/Weibull
- `testthat`, `knitr`, `rmarkdown`

## Architecture

### The `dfr_series_md` Object

Constructor in `R/dfr_series_md.R`:
```r
dfr_series_md(series = NULL, components = NULL,
              par = NULL, n_par = NULL,
              lifetime = "t", lifetime_upper = "t_upper",
              omega = "omega", candset = "x")
```

Class hierarchy: `c("dfr_series_md", "series_md", "likelihood_model")`

Stores a `dfr_dist_series` object (from Phase 1) plus column name configuration for data frames.

### Data Convention

Uses `omega` column (string), NOT numeric `delta`:

| `omega` value | Meaning | Columns used |
|---------------|---------|--------------|
| `"exact"` | Failed at time t, cause masked | `t`, `x1..xm` |
| `"right"` | Right-censored at t | `t` |
| `"left"` | Left-censored: failed before t | `t`, `x1..xm` |
| `"interval"` | Failed in (t, t_upper) | `t`, `t_upper`, `x1..xm` |

Candidate sets: Boolean columns `x1, x2, ..., xm` (TRUE = component j is a candidate cause).

### Closure-Returning Pattern

All methods follow the two-step pattern from `flexhaz`:
1. Call method on model → returns closure
2. Call closure with data and parameters

```r
model <- dfr_series_md(components = list(dfr_exponential(), dfr_exponential()))
ll_fn <- loglik(model)         # returns closure
ll_fn(df, par = c(0.1, 0.2))  # evaluates log-likelihood
```

### Log-Likelihood Computation

**Exact + right-censored** (fast path): Direct hazard/cum_haz calls, no integration.
- Exact: `log(sum_{j in C_i} h_j(t_i)) - H_sys(t_i)`
- Right: `-H_sys(t_i)`

**Left + interval-censored** (integration path): Uses `stats::integrate()` per observation.
- Integrates `[sum_{j in C_i} h_j(t)] * S_sys(t)` over the relevant interval.

### Score and Hessian

Always uses `numDeriv::grad()` and `numDeriv::hessian()` on the log-likelihood. No analytical derivatives in this package — that complexity is deferred to component-level `dfr_dist` objects.

### MLE Fitting

`fit.dfr_series_md` uses `stats::optim()` with numerical score. Returns a `fisher_mle` object with standard methods: `coef()`, `vcov()`, `logLik()`, `confint()`, `summary()`.

Auto-switches from Nelder-Mead to BFGS for single-parameter models.

## Source File Organization

- `R/dfr_series_md.R` — Constructor, type predicate, print method
- `R/loglik.R` — `loglik.dfr_series_md` + `censored_integral` helper
- `R/score_hess.R` — `score.dfr_series_md`, `hess_loglik.dfr_series_md` (numDeriv)
- `R/fit.R` — `fit.dfr_series_md` (optim → fisher_mle)
- `R/rdata.R` — `rdata.dfr_series_md` (data generation)
- `R/methods.R` — `assumptions`, `ncomponents`, `component_hazard`, `conditional_cause_probability`, `cause_probability`
- `R/utils.R` — `decode_candidate_matrix`, `extract_md_data`, `generate_masked_series_data`
- `R/reexports.R` — Re-exports from dependency packages

## Testing

Tests in `tests/testthat/` (10 files):
- `helper-fixtures.R` — Shared helpers (`make_exp_model`, `make_weibull_model`, `exp_series_loglik_analytical`)
- `test-constructor.R` — Class hierarchy, validation, print
- `test-loglik.R` — Finite values, analytical match, monotonicity
- `test-score-hess.R` — Dimensions, near-zero at MLE, negative definite Hessian
- `test-fit.R` — Parameter recovery, convergence, fisher_mle methods
- `test-rdata.R` — Data structure, censoring, masking, custom columns
- `test-censoring.R` — Left/interval/right contributions, analytical match
- `test-cross-validate.R` — Against `maskedcauses` (skipped if not installed)
- `test-methods.R` — assumptions, ncomponents, conditional/marginal cause probabilities
- `test-utils.R` — decode_candidate_matrix, extract_md_data, generate_masked_series_data

### Key Identifiability Constraint

Exponential series systems are **not identifiable** for individual rates from system-level masked data alone — only the sum of rates is identifiable. Tests account for this by checking `sum(coef(result))`.

## Masking Conditions (C1-C2-C3)

- **C1**: The failed component is always in the candidate set
- **C2**: Candidate sets are generated uniformly given the failed component
- **C3**: Masking probabilities are independent of system parameters

These conditions come from the foundational paper at `~/github/papers/masked-causes-in-series-systems/paper.tex` (Theorem 7).
