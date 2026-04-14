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
#' function (gradient) is computed from the same \code{\link{loglik}} closure
#' via \code{\link[numDeriv]{grad}}, and the Hessian at the MLE via
#' \code{\link[numDeriv]{hessian}}. One-parameter problems auto-upgrade from
#' Nelder-Mead to BFGS with a warning, because Nelder-Mead is unreliable in
#' one dimension.
#'
#' @examples
#' \donttest{
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' set.seed(1)
#' df <- rdata(model)(theta = c(0.1, 0.2), n = 200, tau = 10, p = 0)
#' solver <- fit(model)
#' result <- solver(df, par = c(0.15, 0.15))
#' coef(result)
#' }
#' @importFrom generics fit
#' @importFrom likelihood.model fisher_mle
#' @importFrom stats optim
#' @importFrom utils modifyList
#' @method fit dfr_series_md
#' @export
fit.dfr_series_md <- function(object, ...) {
  ll_fn <- loglik(object)

  function(df, par,
           method = c("Nelder-Mead", "BFGS", "SANN", "CG",
                       "L-BFGS-B", "Brent"),
           ..., control = list()) {
    stopifnot(!is.null(par))
    user_picked_method <- !missing(method)
    method <- match.arg(method)
    control <- modifyList(list(fnscale = -1), control)

    if (length(par) == 1 && method == "Nelder-Mead") {
      if (user_picked_method)
        warning("Nelder-Mead is unreliable for 1-parameter problems; ",
                "switching to BFGS.")
      method <- "BFGS"
    }

    obj <- function(p) ll_fn(df, p)
    gr  <- function(p) numDeriv::grad(obj, p)

    sol <- optim(
      par = par, fn = obj,
      gr = if (method == "SANN") NULL else gr,
      hessian = FALSE, method = method, control = control
    )

    fisher_mle(
      par = sol$par,
      loglik_val = sol$value,
      hessian = numDeriv::hessian(obj, sol$par),
      score_val = numDeriv::grad(obj, sol$par),
      nobs = nrow(df),
      converged = (sol$convergence == 0)
    )
  }
}
