#' Assumptions for masked-cause DFR series systems
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (unused).
#' @return Character vector of model assumptions.
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
#' @importFrom serieshaz component_hazard
#' @method component_hazard dfr_series_md
#' @export
component_hazard.dfr_series_md <- function(x, j, ...) {
  component_hazard(x$series, j, ...)
}


# =========================================================================
# Series MD concept generics
#
# These generics mirror those in maskedcauses but are
# defined here to avoid a hard dependency. dfr_series_md inherits
# "series_md" in its class hierarchy so dispatch works correctly.
# =========================================================================

#' Conditional cause-of-failure probability
#'
#' Returns a closure computing \eqn{P(K=j | T=t, \theta)} for all components,
#' conditional on a specific failure time t. By Theorem 6 of the foundational
#' paper, this equals \eqn{h_j(t; \theta) / \sum_l h_l(t; \theta)}.
#'
#' @param model A likelihood model object.
#' @param ... Additional arguments passed to the returned closure.
#' @return A function with signature \code{function(t, par, ...)} returning an
#'   n x m matrix where column j gives P(K=j | T=t, theta).
#' @export
conditional_cause_probability <- function(model, ...) {
  UseMethod("conditional_cause_probability")
}


#' @describeIn conditional_cause_probability Method for masked-cause DFR
#'   series systems using component hazard ratios.
#' @importFrom algebraic.dist hazard
#' @method conditional_cause_probability dfr_series_md
#' @export
conditional_cause_probability.dfr_series_md <- function(model, ...) {
  series <- model$series
  m <- series$m
  h_fns <- lapply(series$components, hazard)
  layout <- series$layout

  function(t, par, ...) {
    n <- length(t)
    H <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      H[, j] <- h_fns[[j]](t, par = par[layout[[j]]])
    }
    row_sums <- rowSums(H)
    H / row_sums
  }
}


#' Marginal cause-of-failure probability
#'
#' Returns a closure computing \eqn{P(K=j | \theta)} for all components,
#' marginalized over the system failure time T. By Theorem 5, this equals
#' \eqn{E_T[P(K=j | T, \theta)]}.
#'
#' @param model A likelihood model object.
#' @param ... Additional arguments passed to the returned closure.
#' @return A function with signature \code{function(par, ...)} returning an
#'   m-vector where element j gives P(K=j | theta).
#' @export
cause_probability <- function(model, ...) {
  UseMethod("cause_probability")
}


#' @describeIn cause_probability Method for masked-cause DFR series systems
#'   using Monte Carlo integration.
#' @method cause_probability dfr_series_md
#' @export
cause_probability.dfr_series_md <- function(model, ...) {
  cond_fn <- conditional_cause_probability(model)
  rdata_fn <- rdata(model)
  defaults <- extract_model_defaults(model)

  function(par, n_mc = 10000, tau = Inf, p = 0, ...) {
    df <- rdata_fn(theta = par, n = n_mc, tau = tau, p = p)
    omega_col <- defaults$omega
    if (omega_col %in% colnames(df)) {
      is_exact <- df[[omega_col]] == "exact"
    } else {
      is_exact <- rep(TRUE, nrow(df))
    }
    t_exact <- df[[defaults$lifetime]][is_exact]
    if (length(t_exact) == 0) {
      m <- ncomponents(model)
      return(rep(NA_real_, m))
    }
    probs <- cond_fn(t_exact, par, ...)
    colMeans(probs)
  }
}
