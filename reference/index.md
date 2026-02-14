# Package index

## Model Constructor

Create masked-cause series system likelihood models

- [`dfr_series_md()`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
  : Masked-Cause Likelihood Model for DFR Series Systems
- [`is_dfr_series_md()`](https://queelius.github.io/maskedhaz/reference/is_dfr_series_md.md)
  : Test whether an object is a dfr_series_md
- [`print(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/print.dfr_series_md.md)
  : Print method for dfr_series_md

## Likelihood Inference

Log-likelihood, score, Hessian, and MLE fitting

- [`loglik(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/loglik.dfr_series_md.md)
  : Log-likelihood for masked-cause DFR series systems
- [`score(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/score.dfr_series_md.md)
  : Score function for masked-cause DFR series systems
- [`hess_loglik(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/hess_loglik.dfr_series_md.md)
  : Hessian of log-likelihood for masked-cause DFR series systems
- [`fit(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/fit.dfr_series_md.md)
  : MLE fitting for masked-cause DFR series systems

## Data Generation

Simulate masked series system failure data

- [`rdata(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/rdata.dfr_series_md.md)
  : Random data generation for masked-cause DFR series systems

## Diagnostics

Cause-of-failure probabilities and model inspection

- [`conditional_cause_probability()`](https://queelius.github.io/maskedhaz/reference/conditional_cause_probability.md)
  : Conditional cause-of-failure probability
- [`cause_probability()`](https://queelius.github.io/maskedhaz/reference/cause_probability.md)
  : Marginal cause-of-failure probability
- [`assumptions(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/assumptions.dfr_series_md.md)
  : Assumptions for masked-cause DFR series systems
- [`ncomponents(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/ncomponents.dfr_series_md.md)
  : Number of components in a masked-cause DFR series system
- [`component_hazard(`*`<dfr_series_md>`*`)`](https://queelius.github.io/maskedhaz/reference/component_hazard.dfr_series_md.md)
  : Component hazard for a masked-cause DFR series system
