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
