#' Log-likelihood for masked-cause DFR series systems
#'
#' Returns a log-likelihood function for a series system with masked component
#' cause of failure. Supports four observation types: exact failures,
#' right-censored, left-censored, and interval-censored.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A function with signature \code{function(df, par, ...)} that computes
#'   the log-likelihood.
#'
#' @details
#' Log-likelihood contributions by observation type:
#' \describe{
#'   \item{Exact (\eqn{\omega = } "exact")}{
#'     \eqn{\log L_i = \log(\sum_{j \in C_i} h_j(t_i)) - H_{sys}(t_i)}}
#'   \item{Right-censored (\eqn{\omega = } "right")}{
#'     \eqn{\log L_i = -H_{sys}(t_i)}}
#'   \item{Left-censored (\eqn{\omega = } "left")}{
#'     \eqn{\log L_i = \log \int_0^{t_i} [\sum_{j \in C_i} h_j(u)] S_{sys}(u) \, du}}
#'   \item{Interval-censored (\eqn{\omega = } "interval")}{
#'     \eqn{\log L_i = \log \int_{t_i}^{t_{upper,i}} [\sum_{j \in C_i} h_j(u)] S_{sys}(u) \, du}}
#' }
#'
#' The exact and right-censored paths use direct hazard/cumulative hazard
#' calls. Left and interval censoring require numerical integration via
#' \code{\link[stats]{integrate}}.
#'
#' @importFrom likelihood.model loglik
#' @importFrom algebraic.dist hazard
#' @importFrom dfr.dist cum_haz
#' @importFrom stats integrate
#' @method loglik dfr_series_md
#' @export
loglik.dfr_series_md <- function(model, ...) {
  defaults <- extract_model_defaults(model)
  series <- model$series
  comps <- series$components
  layout <- series$layout
  m <- series$m

  # Pre-build hazard and cum_haz closures for each component
  h_fns <- lapply(comps, hazard)
  H_fns <- lapply(comps, cum_haz)

  function(df, par, ...) {
    if (any(par <= 0)) return(-Inf)

    d <- extract_md_data(df, defaults$lifetime, defaults$omega,
                         defaults$candset, defaults$lifetime_upper)
    ll <- 0

    # Vectorized system cumulative hazard for exact + right observations
    er_idx <- which(d$omega %in% c("exact", "right"))
    if (length(er_idx) > 0) {
      for (i in er_idx) {
        t_i <- d$t[i]
        H_sys <- 0
        for (j in seq_len(m)) {
          H_sys <- H_sys + H_fns[[j]](t_i, par = par[layout[[j]]])
        }
        ll <- ll - H_sys

        if (d$omega[i] == "exact") {
          h_cand <- 0
          for (j in which(d$C[i, ])) {
            h_cand <- h_cand + h_fns[[j]](t_i, par = par[layout[[j]]])
          }
          if (h_cand <= 0) return(-Inf)
          ll <- ll + log(h_cand)
        }
      }
    }

    # Left-censored observations: integrate over (0, t_i)
    for (i in which(d$omega == "left")) {
      val <- censored_integral(h_fns, H_fns, layout, par, d$C[i, ],
                               0, d$t[i], m)
      if (val <= 0) return(-Inf)
      ll <- ll + log(val)
    }

    # Interval-censored observations: integrate over (t_i, t_upper_i)
    for (i in which(d$omega == "interval")) {
      val <- censored_integral(h_fns, H_fns, layout, par, d$C[i, ],
                               d$t[i], d$t_upper[i], m)
      if (val <= 0) return(-Inf)
      ll <- ll + log(val)
    }

    ll
  }
}


#' Compute integral for left/interval censored contributions
#'
#' Integrates `[sum_{j in C_i} h_j(t)] * S_sys(t)` over `(lower, upper)`.
#'
#' @param h_fns list of hazard closures
#' @param H_fns list of cumulative hazard closures
#' @param layout parameter layout
#' @param par parameter vector
#' @param C_i logical vector of candidate components
#' @param lower lower integration bound
#' @param upper upper integration bound
#' @param m number of components
#' @return numeric integral value
#' @keywords internal
censored_integral <- function(h_fns, H_fns, layout, par, C_i,
                              lower, upper, m) {
  integrand <- function(t) {
    h_cand <- 0
    H_sys <- 0
    for (j in seq_len(m)) {
      par_j <- par[layout[[j]]]
      H_sys <- H_sys + H_fns[[j]](t, par = par_j)
      if (C_i[j]) h_cand <- h_cand + h_fns[[j]](t, par = par_j)
    }
    h_cand * exp(-H_sys)
  }
  result <- integrate(integrand, lower, upper, rel.tol = 1e-8,
                      stop.on.error = FALSE)
  result$value
}
