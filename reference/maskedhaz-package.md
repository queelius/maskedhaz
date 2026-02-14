# maskedhaz: Masked-Cause Likelihood Models for DFR Series Systems

Likelihood-based inference for series systems with masked component
cause of failure, using arbitrary dynamic failure rate (DFR) component
distributions. Computes log-likelihood, score, Hessian, and MLE for
masked data satisfying conditions C1, C2, C3 under general component
hazard functions.

## Details

The maskedhaz package provides likelihood-based inference for series
systems with masked component cause of failure, using arbitrary dynamic
failure rate (DFR) component distributions from serieshaz.

A series system fails when any component fails, but the causing
component may be unknown (masked). Given candidate sets satisfying
conditions C1, C2, C3, this package computes log-likelihood, score,
Hessian, and MLE for the component parameters.

## Package functions

- [`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md):

  Constructor: create a masked-cause likelihood model

- [`is_dfr_series_md`](https://queelius.github.io/maskedhaz/reference/is_dfr_series_md.md):

  Type predicate

- [`loglik`](https://queelius.github.io/likelihood.model/reference/loglik.html):

  Log-likelihood

- [`score`](https://queelius.github.io/likelihood.model/reference/score.html):

  Score function

- [`hess_loglik`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html):

  Hessian

- [`fit`](https://generics.r-lib.org/reference/fit.html):

  MLE fitting

- [`rdata`](https://queelius.github.io/likelihood.model/reference/rdata.html):

  Data generation

## See also

[`dfr_series_md`](https://queelius.github.io/maskedhaz/reference/dfr_series_md.md)
for the constructor,
[`dfr_dist_series`](https://queelius.github.io/dfr.dist.series/reference/dfr_dist_series.html)
for the series distribution,
[`loglik`](https://queelius.github.io/likelihood.model/reference/loglik.html)
for the likelihood interface

## Author

**Maintainer**: Alexander Towell <lex@metafunctor.com>
