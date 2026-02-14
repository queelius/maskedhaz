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

# Re-exports from serieshaz
#' @importFrom serieshaz dfr_dist_series is_dfr_dist_series ncomponents
#'   component component_hazard param_layout sample_components
#' @export
serieshaz::dfr_dist_series

#' @export
serieshaz::is_dfr_dist_series

#' @export
serieshaz::ncomponents

#' @export
serieshaz::component

#' @export
serieshaz::component_hazard

#' @export
serieshaz::param_layout

#' @export
serieshaz::sample_components

# Re-exports from flexhaz
#' @importFrom flexhaz dfr_exponential dfr_weibull dfr_gompertz
#'   dfr_loglogistic cum_haz
#' @export
flexhaz::dfr_exponential

#' @export
flexhaz::dfr_weibull

#' @export
flexhaz::dfr_gompertz

#' @export
flexhaz::dfr_loglogistic

#' @export
flexhaz::cum_haz

# Re-exports from algebraic.dist (via serieshaz)
#' @importFrom algebraic.dist hazard surv cdf sampler params
#' @export
algebraic.dist::hazard

#' @export
algebraic.dist::surv

#' @export
algebraic.dist::sampler

#' @export
algebraic.dist::params
