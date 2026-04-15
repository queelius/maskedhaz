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
#' The exact and right-censored paths use vectorized hazard / cumulative hazard
#' calls. Left and interval censoring require per-row numerical integration
#' via \code{\link[stats]{integrate}}.
#'
#' The returned closure caches validated and decoded masked-data extracted
#' from the data frame across repeated calls with the same \code{df}, so that
#' the O(n) validation cost is paid only once per \code{optim}/\code{numDeriv}
#' sweep. The cache is per-closure, kept in the closure's enclosing
#' environment. This is safe for sequential use; if you share the same
#' closure object across forked workers (e.g. \code{parallel::mcparallel}),
#' concurrent writes to the cache are possible but only affect performance,
#' not correctness.
#'
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' set.seed(1)
#' df <- rdata(model)(theta = c(0.1, 0.2), n = 50, tau = 10, p = 0.3)
#' ll_fn <- loglik(model)
#' ll_fn(df, par = c(0.1, 0.2))
#' @importFrom likelihood.model loglik
#' @importFrom algebraic.dist hazard
#' @importFrom flexhaz cum_haz
#' @importFrom stats integrate
#' @method loglik dfr_series_md
#' @export
loglik.dfr_series_md <- function(model, ...) {
  defaults <- extract_model_defaults(model)
  series <- model$series
  layout <- series$layout
  m <- series$m

  h_fns <- lapply(series$components, hazard)
  H_fns <- lapply(series$components, cum_haz)

  # Cache the validated/decoded data when called repeatedly with the same df
  # (e.g. inside optim → numDeriv). identical() is O(1) for same-object input.
  cached_df <- NULL
  cached_d  <- NULL

  function(df, par, ...) {
    if (any(par <= 0)) return(-Inf)

    if (!identical(df, cached_df)) {
      cached_d  <<- extract_md_data(df, defaults$lifetime, defaults$omega,
                                    defaults$candset, defaults$lifetime_upper)
      cached_df <<- df
    }
    d <- cached_d

    # Slice per-component parameter subvectors once per evaluation.
    par_by_comp <- lapply(layout, function(ix) par[ix])

    ll <- 0
    er_idx <- which(d$omega %in% c(OMEGA_EXACT, OMEGA_RIGHT))
    if (length(er_idx) > 0) {
      t_er <- d$t[er_idx]
      H_sys_er <- numeric(length(t_er))
      for (j in seq_len(m)) {
        H_sys_er <- H_sys_er + H_fns[[j]](t_er, par = par_by_comp[[j]])
      }
      ll <- ll - sum(H_sys_er)

      exact_mask <- d$omega[er_idx] == OMEGA_EXACT
      if (any(exact_mask)) {
        t_ex <- t_er[exact_mask]
        C_ex <- d$C[er_idx[exact_mask], , drop = FALSE]
        h_cand <- numeric(length(t_ex))
        for (j in seq_len(m)) {
          if (any(C_ex[, j])) {
            hj <- h_fns[[j]](t_ex, par = par_by_comp[[j]])
            h_cand <- h_cand + hj * C_ex[, j]
          }
        }
        if (any(h_cand <= 0)) return(-Inf)
        ll <- ll + sum(log(h_cand))
      }
    }

    # Left and interval censoring share the same integrand; only bounds differ.
    for (om in c(OMEGA_LEFT, OMEGA_INTERVAL)) {
      idx <- which(d$omega == om)
      for (i in idx) {
        lower <- if (om == OMEGA_LEFT) 0 else d$t[i]
        upper <- if (om == OMEGA_LEFT) d$t[i] else d$t_upper[i]
        val <- censored_integral(h_fns, H_fns, par_by_comp, d$C[i, ],
                                 lower, upper, m)
        if (!is.finite(val) || val <= 0) return(-Inf)
        ll <- ll + log(val)
      }
    }

    if (!is.finite(ll)) return(-Inf)
    ll
  }
}


#' Compute integral for left/interval censored contributions
#'
#' Integrates `[sum_{j in C_i} h_j(t)] * S_sys(t)` over `(lower, upper)`.
#' Returns `NA_real_` when \code{\link[stats]{integrate}} fails to converge,
#' so the caller can steer the optimizer away from pathological regions.
#'
#' @param h_fns list of hazard closures
#' @param H_fns list of cumulative hazard closures
#' @param par_by_comp list of per-component parameter subvectors
#' @param C_i logical vector of candidate components
#' @param lower lower integration bound
#' @param upper upper integration bound
#' @param m number of components
#' @return numeric integral value, or \code{NA_real_} on non-convergence.
#' @keywords internal
censored_integral <- function(h_fns, H_fns, par_by_comp, C_i,
                              lower, upper, m) {
  integrand <- function(t) {
    h_cand <- 0
    H_sys <- 0
    for (j in seq_len(m)) {
      par_j <- par_by_comp[[j]]
      H_sys <- H_sys + H_fns[[j]](t, par = par_j)
      if (C_i[j]) h_cand <- h_cand + h_fns[[j]](t, par = par_j)
    }
    h_cand * exp(-H_sys)
  }
  result <- try(integrate(integrand, lower, upper, rel.tol = 1e-6,
                          stop.on.error = FALSE), silent = TRUE)
  if (inherits(result, "try-error") || result$message != "OK") return(NA_real_)
  result$value
}
