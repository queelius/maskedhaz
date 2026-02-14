# Re-exports from likelihood.model
#' @importFrom likelihood.model loglik score hess_loglik assumptions rdata
#' @export
likelihood.model::loglik

#' @export
likelihood.model::score

#' @export
likelihood.model::hess_loglik

#' @export
likelihood.model::assumptions

#' @export
likelihood.model::rdata

# Re-exports from generics
#' @importFrom generics fit
#' @export
generics::fit

# Re-exports from dfr.dist.series
#' @importFrom dfr.dist.series dfr_dist_series is_dfr_dist_series ncomponents
#'   component component_hazard param_layout sample_components
#' @export
dfr.dist.series::dfr_dist_series

#' @export
dfr.dist.series::is_dfr_dist_series

#' @export
dfr.dist.series::ncomponents

#' @export
dfr.dist.series::component

#' @export
dfr.dist.series::component_hazard

#' @export
dfr.dist.series::param_layout

#' @export
dfr.dist.series::sample_components

# Re-exports from dfr.dist
#' @importFrom dfr.dist dfr_exponential dfr_weibull dfr_gompertz
#'   dfr_loglogistic cum_haz
#' @export
dfr.dist::dfr_exponential

#' @export
dfr.dist::dfr_weibull

#' @export
dfr.dist::dfr_gompertz

#' @export
dfr.dist::dfr_loglogistic

#' @export
dfr.dist::cum_haz

# Re-exports from algebraic.dist (via dfr.dist.series)
#' @importFrom algebraic.dist hazard surv cdf sampler params
#' @export
algebraic.dist::hazard

#' @export
algebraic.dist::surv

#' @export
algebraic.dist::sampler

#' @export
algebraic.dist::params
