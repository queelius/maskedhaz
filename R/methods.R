#' Assumptions for masked-cause DFR series systems
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (unused).
#' @return Character vector of model assumptions.
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' assumptions(model)
#' @importFrom likelihood.model assumptions
#' @method assumptions dfr_series_md
#' @export
assumptions.dfr_series_md <- function(model, ...) {
  c(
    "iid observations",
    "series system configuration",
    "component independence: component lifetimes are independent",
    "non-negative hazard: h_j(t) >= 0 for all j, t > 0",
    "C1: failed component is in candidate set with probability 1",
    "C2: uniform probability for candidate sets given component cause",
    "C3: masking probabilities independent of system parameters"
  )
}


#' Number of components in a masked-cause DFR series system
#'
#' @param x A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (unused).
#' @return Integer, the number of components.
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2), dfr_exponential(0.3)
#' ))
#' ncomponents(model)   # 3
#' @importFrom serieshaz ncomponents
#' @method ncomponents dfr_series_md
#' @export
ncomponents.dfr_series_md <- function(x, ...) {
  ncomponents(x$series)
}


#' Component hazard for a masked-cause DFR series system
#'
#' @param x A \code{\link{dfr_series_md}} object.
#' @param j Component index.
#' @param ... Additional arguments passed to the closure.
#' @return A closure computing component j's hazard.
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' h1 <- component_hazard(model, 1)
#' h1(t = 5, par = 0.1)   # 0.1 (constant exponential hazard)
#' @importFrom serieshaz component_hazard
#' @method component_hazard dfr_series_md
#' @export
component_hazard.dfr_series_md <- function(x, j, ...) {
  component_hazard(x$series, j, ...)
}


#' Conditional cause-of-failure probability for DFR series systems
#'
#' Method for \code{\link[maskedcauses]{conditional_cause_probability}} that
#' returns a closure computing \eqn{P(K=j \mid T=t, \theta)} for each
#' component. By Theorem 6 of the foundational paper, this equals
#' \eqn{h_j(t; \theta) / \sum_l h_l(t; \theta)}.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments passed to the returned closure.
#' @return A function with signature \code{function(t, par, ...)} returning an
#'   n x m matrix where column j gives P(K=j | T=t, theta).
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2), dfr_exponential(0.3)
#' ))
#' ccp_fn <- conditional_cause_probability(model)
#' ccp_fn(t = c(1, 5, 10), par = c(0.1, 0.2, 0.3))
#' @importFrom maskedcauses conditional_cause_probability
#' @importFrom algebraic.dist hazard
#' @method conditional_cause_probability dfr_series_md
#' @export
conditional_cause_probability.dfr_series_md <- function(model, ...) {
  series <- model$series
  m <- series$m
  h_fns <- lapply(series$components, hazard)
  layout <- series$layout

  function(t, par, ...) {
    H <- matrix(0, nrow = length(t), ncol = m)
    for (j in seq_len(m)) {
      H[, j] <- h_fns[[j]](t, par = par[layout[[j]]])
    }
    H / rowSums(H)
  }
}


#' Marginal cause-of-failure probability for DFR series systems
#'
#' Method for \code{\link[maskedcauses]{cause_probability}} that returns a
#' closure computing \eqn{P(K=j \mid \theta)} for each component, marginalized
#' over the system failure time T via Monte Carlo integration. By Theorem 5,
#' this equals \eqn{E_T[P(K=j \mid T, \theta)]}.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments passed to the returned closure.
#' @return A function with signature \code{function(par, ...)} returning an
#'   m-vector where element j gives P(K=j | theta).
#' @examples
#' \donttest{
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2), dfr_exponential(0.3)
#' ))
#' cp_fn <- cause_probability(model)
#' set.seed(1)
#' cp_fn(par = c(0.1, 0.2, 0.3), n_mc = 2000)
#' }
#' @importFrom maskedcauses cause_probability
#' @method cause_probability dfr_series_md
#' @export
cause_probability.dfr_series_md <- function(model, ...) {
  cond_fn <- conditional_cause_probability(model)
  rdata_fn <- rdata(model)
  defaults <- extract_model_defaults(model)

  function(par, n_mc = 10000, tau = Inf, p = 0, ...) {
    df <- rdata_fn(theta = par, n = n_mc, tau = tau, p = p)
    is_exact <- if (defaults$omega %in% colnames(df)) {
      df[[defaults$omega]] == OMEGA_EXACT
    } else {
      rep(TRUE, nrow(df))
    }
    t_exact <- df[[defaults$lifetime]][is_exact]
    if (length(t_exact) == 0) return(rep(NA_real_, ncomponents(model)))
    colMeans(cond_fn(t_exact, par, ...))
  }
}
