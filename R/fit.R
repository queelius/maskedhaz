#' MLE fitting for masked-cause DFR series systems
#'
#' Returns a solver function that finds the maximum likelihood estimates for
#' component parameters given masked series system data.
#'
#' @param object A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A solver function with signature
#'   \code{function(df, par, method = "Nelder-Mead", ..., control = list())}
#'   that returns a \code{\link[likelihood.model]{fisher_mle}} object.
#'
#' @details
#' Uses \code{\link[stats]{optim}} to maximize the log-likelihood. The score
#' function (gradient) is provided for gradient-based methods. The Hessian at
#' the MLE is computed for variance-covariance estimation.
#'
#' @importFrom generics fit
#' @importFrom likelihood.model fisher_mle
#' @importFrom stats optim
#' @importFrom utils modifyList
#' @method fit dfr_series_md
#' @export
fit.dfr_series_md <- function(object, ...) {
  ll_fn <- loglik(object)
  s_fn <- score(object)
  H_fn <- hess_loglik(object)

  function(df, par,
           method = c("Nelder-Mead", "BFGS", "SANN", "CG",
                       "L-BFGS-B", "Brent"),
           ..., control = list()) {
    stopifnot(!is.null(par))
    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)
    method <- match.arg(method)

    if (length(par) == 1 && method == "Nelder-Mead")
      method <- "BFGS"

    sol <- optim(
      par = par,
      fn = function(p) ll_fn(df, p),
      gr = if (method == "SANN") NULL else function(p) s_fn(df, p),
      hessian = FALSE,
      method = method,
      control = control
    )

    hessian <- H_fn(df, sol$par)
    score_at_mle <- s_fn(df, sol$par)

    fisher_mle(
      par = sol$par,
      loglik_val = sol$value,
      hessian = hessian,
      score_val = score_at_mle,
      nobs = nrow(df),
      converged = (sol$convergence == 0)
    )
  }
}
