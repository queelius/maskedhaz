#' @keywords internal
#' @details
#' The \pkg{maskedhaz} package provides likelihood-based inference for
#' series systems with masked component cause of failure, using arbitrary
#' dynamic failure rate (DFR) component distributions from
#' \pkg{serieshaz}.
#'
#' A series system fails when any component fails, but the causing component
#' may be unknown (masked). Given candidate sets satisfying conditions C1, C2,
#' C3, this package computes log-likelihood, score, Hessian, and MLE for the
#' component parameters.
#'
#' @section Package functions:
#' \describe{
#'   \item{\code{\link{dfr_series_md}}}{Constructor: create a masked-cause
#'     likelihood model}
#'   \item{\code{\link{is_dfr_series_md}}}{Type predicate}
#'   \item{\code{\link[likelihood.model]{loglik}}}{Log-likelihood}
#'   \item{\code{\link[likelihood.model]{score}}}{Score function}
#'   \item{\code{\link[likelihood.model]{hess_loglik}}}{Hessian}
#'   \item{\code{\link[generics]{fit}}}{MLE fitting}
#'   \item{\code{\link[likelihood.model]{rdata}}}{Data generation}
#' }
#'
#' @seealso
#' \code{\link{dfr_series_md}} for the constructor,
#' \code{\link[serieshaz]{dfr_dist_series}} for the series distribution,
#' \code{\link[likelihood.model]{loglik}} for the likelihood interface
"_PACKAGE"

#' Masked-Cause Likelihood Model for DFR Series Systems
#'
#' Constructs a likelihood model for series systems with masked component cause
#' of failure, where components are arbitrary \code{\link[flexhaz]{dfr_dist}}
#' distributions. Supports exact, right-censored, left-censored, and
#' interval-censored observations with candidate sets satisfying C1-C2-C3.
#'
#' @param series A \code{\link[serieshaz]{dfr_dist_series}} object.
#'   Ignored if \code{components} is provided.
#' @param components A list of \code{\link[flexhaz]{dfr_dist}} objects. If
#'   provided, a \code{dfr_dist_series} is built from these.
#' @param par Optional concatenated parameter vector.
#' @param n_par Optional integer vector of parameter counts per component.
#' @param lifetime Column name for system lifetime (default \code{"t"}).
#' @param lifetime_upper Column name for interval upper bound (default
#'   \code{"t_upper"}).
#' @param omega Column name for observation type (default \code{"omega"}).
#' @param candset Column prefix for candidate set indicators (default
#'   \code{"x"}).
#'
#' @return An object of class
#'   \code{c("dfr_series_md", "series_md", "likelihood_model")}.
#'
#' @details
#' The model computes the masked-cause log-likelihood for series systems where
#' the system lifetime is the minimum of independent component lifetimes, and
#' the causing component is partially observed through candidate sets.
#'
#' \strong{Observation types} (stored in the \code{omega} column):
#' \describe{
#'   \item{\code{"exact"}}{Failed at time t, cause masked among candidates}
#'   \item{\code{"right"}}{Right-censored: survived past time t}
#'   \item{\code{"left"}}{Left-censored: failed before time t}
#'   \item{\code{"interval"}}{Failed in interval (t, t_upper)}
#' }
#'
#' \strong{Masking conditions}:
#' \describe{
#'   \item{C1}{Failed component is in candidate set with probability 1}
#'   \item{C2}{Uniform probability for candidate sets given component cause}
#'   \item{C3}{Masking probabilities independent of system parameters}
#' }
#'
#' @examples
#' \donttest{
#' library(flexhaz)
#' library(serieshaz)
#'
#' # From components
#' model <- dfr_series_md(components = list(
#'     dfr_exponential(0.1),
#'     dfr_exponential(0.2),
#'     dfr_exponential(0.3)
#' ))
#'
#' # From pre-built series
#' sys <- dfr_dist_series(list(
#'     dfr_weibull(shape = 2, scale = 100),
#'     dfr_exponential(0.05)
#' ))
#' model2 <- dfr_series_md(series = sys)
#' }
#'
#' @seealso
#' \code{\link{is_dfr_series_md}} for the type predicate,
#' \code{\link[serieshaz]{dfr_dist_series}} for the series distribution,
#' \code{\link[likelihood.model]{loglik}} for the likelihood interface
#'
#' @importFrom serieshaz dfr_dist_series is_dfr_dist_series
#' @export
dfr_series_md <- function(series = NULL, components = NULL,
                          par = NULL, n_par = NULL,
                          lifetime = "t", lifetime_upper = "t_upper",
                          omega = "omega", candset = "x") {
  if (is.null(series) && !is.null(components)) {
    series <- dfr_dist_series(components, par = par, n_par = n_par)
  }
  if (is.null(series))
    stop("Either 'series' or 'components' must be provided")
  if (!is_dfr_dist_series(series))
    stop("'series' must be a dfr_dist_series object")

  structure(list(
    series = series,
    lifetime = lifetime,
    lifetime_upper = lifetime_upper,
    omega = omega,
    candset = candset
  ), class = c("dfr_series_md", "series_md", "likelihood_model"))
}


#' Test whether an object is a dfr_series_md
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_dfr_series_md <- function(x) {
  inherits(x, "dfr_series_md")
}


#' Print method for dfr_series_md
#'
#' @param x A \code{dfr_series_md} object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.dfr_series_md <- function(x, ...) {
  m <- x$series$m
  cat(sprintf("Masked-cause likelihood model (%d-component series)\n", m))
  for (j in seq_len(m)) {
    np <- x$series$n_par[j]
    par_j <- if (!is.null(x$series$par)) {
      x$series$par[x$series$layout[[j]]]
    }
    par_str <- if (!is.null(par_j)) {
      paste(format(par_j, digits = 4), collapse = ", ")
    } else {
      "unknown"
    }
    cat(sprintf("  Component %d: %d param(s) [%s]\n", j, np, par_str))
  }
  cat("Data columns:", x$lifetime, "(lifetime),",
      x$omega, "(type),", x$candset, "* (candidates)\n")
  invisible(x)
}
