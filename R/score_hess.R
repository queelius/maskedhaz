#' Score function for masked-cause DFR series systems
#'
#' Returns a score (gradient) function computed via numerical differentiation
#' of the log-likelihood using \code{\link[numDeriv]{grad}}.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A function with signature \code{function(df, par, ...)} returning
#'   the gradient vector.
#'
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' set.seed(1)
#' df <- rdata(model)(theta = c(0.1, 0.2), n = 50, tau = 10, p = 0.3)
#' s_fn <- score(model)
#' s_fn(df, par = c(0.1, 0.2))
#' @importFrom likelihood.model score
#' @importFrom numDeriv grad
#' @method score dfr_series_md
#' @export
score.dfr_series_md <- function(model, ...) {
  ll_fn <- loglik(model)
  function(df, par, ...) {
    numDeriv::grad(func = function(p) ll_fn(df, p), x = par)
  }
}


#' Hessian of log-likelihood for masked-cause DFR series systems
#'
#' Returns a Hessian function computed via numerical differentiation of the
#' log-likelihood using \code{\link[numDeriv]{hessian}}.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A function with signature \code{function(df, par, ...)} returning
#'   the Hessian matrix.
#'
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' set.seed(1)
#' df <- rdata(model)(theta = c(0.1, 0.2), n = 50, tau = 10, p = 0.3)
#' H_fn <- hess_loglik(model)
#' H_fn(df, par = c(0.1, 0.2))
#' @importFrom likelihood.model hess_loglik
#' @importFrom numDeriv hessian
#' @method hess_loglik dfr_series_md
#' @export
hess_loglik.dfr_series_md <- function(model, ...) {
  ll_fn <- loglik(model)
  function(df, par, ...) {
    numDeriv::hessian(func = function(p) ll_fn(df, p), x = par)
  }
}
