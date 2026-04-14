# maskedhaz 0.1.0

Initial CRAN release.

## Features

* `dfr_series_md()` constructor builds a masked-cause likelihood model for a
  series system whose components are arbitrary `dfr_dist` objects from the
  `flexhaz` and `serieshaz` packages. The resulting model implements the
  `series_md` protocol defined in `maskedcauses`.
* Full support for four observation types via the `omega` column:
  `"exact"`, `"right"`, `"left"`, `"interval"`. Left- and interval-censored
  contributions use `stats::integrate()`; exact and right-censored rows use
  vectorised closed-form expressions.
* `loglik()`, `score()`, `hess_loglik()`, `fit()`, and `rdata()` methods
  that dispatch through the `likelihood.model` generics. `fit()` returns a
  `fisher_mle` object, which realises the `mle_fit` and `algebraic.dist`
  interfaces, so standard MLE diagnostics (`coef`, `vcov`, `confint`, `se`,
  `bias`, `observed_fim`, `as_dist`, `sampler`, `expectation`) all work
  uniformly.
* Methods for the `maskedcauses` domain generics
  `conditional_cause_probability()` and `cause_probability()`.
* Cross-validated against the closed-form exponential-series likelihood
  from `maskedcauses` to confirm the numerical integration path matches
  analytical results to integrator tolerance.

## Vignettes

* `maskedhaz` (overview): protocol, component, and MLE-result stacks;
  quick tour from construction to diagnostics.
* `custom-components`: mixed-distribution series with Weibull, exponential,
  Gompertz, and log-logistic components.
* `censoring-and-masking`: the four observation types with worked examples
  and cross-validation against `maskedcauses`.
