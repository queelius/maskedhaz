# maskedhaz 0.1.0

* Initial release
* `dfr_series_md()` constructor for masked-cause likelihood models
* Log-likelihood supporting exact, right, left, and interval censoring
* Score and Hessian via `numDeriv`
* MLE fitting via `optim` returning `fisher_mle` objects
* Data generation with configurable censoring and masking
* Conditional and marginal cause-of-failure probabilities
* Cross-validated against `maskedcauses` for exponential components
