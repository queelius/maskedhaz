#' Random data generation for masked-cause DFR series systems
#'
#' Returns a function that generates random masked series system data from the
#' model's data-generating process (DGP). Uses
#' \code{\link[dfr.dist.series]{sample_components}} for component lifetimes
#' and applies right-censoring and masking satisfying C1-C2-C3.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A function with signature \code{function(theta, n, tau = Inf, p = 0, ...)}
#'   that returns a data frame with columns for lifetime, observation type, and
#'   candidate sets.
#'
#' @importFrom likelihood.model rdata
#' @importFrom dfr.dist.series sample_components
#' @method rdata dfr_series_md
#' @export
rdata.dfr_series_md <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, ...) {
    series <- model$series
    m <- series$m
    if (any(theta <= 0)) stop("All parameters must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    comp_lifetimes <- sample_components(series, n, par = theta)
    generate_masked_series_data(comp_lifetimes, n, m, tau, p,
                                defaults$lifetime, defaults$omega,
                                defaults$candset)
  }
}
