#' Random data generation for masked-cause DFR series systems
#'
#' Returns a function that generates random masked series system data from the
#' model's data-generating process (DGP). Uses
#' \code{\link[serieshaz]{sample_components}} for component lifetimes
#' and applies right-censoring and masking satisfying C1-C2-C3.
#'
#' @param model A \code{\link{dfr_series_md}} object.
#' @param ... Additional arguments (currently unused).
#' @return A function with signature \code{function(theta, n, tau = Inf, p = 0, ...)}
#'   that returns a data frame with columns for lifetime, observation type, and
#'   candidate sets.
#'
#' @examples
#' model <- dfr_series_md(components = list(
#'   dfr_exponential(0.1), dfr_exponential(0.2)
#' ))
#' set.seed(1)
#' df <- rdata(model)(theta = c(0.1, 0.2), n = 20, tau = 10, p = 0.3)
#' head(df)
#' @importFrom likelihood.model rdata
#' @importFrom serieshaz sample_components
#' @method rdata dfr_series_md
#' @export
rdata.dfr_series_md <- function(model, ...) {
  defaults <- extract_model_defaults(model)
  series <- model$series
  m <- series$m

  function(theta, n, tau = Inf, p = 0, ...) {
    if (any(theta <= 0)) stop("All parameters must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    comp_lifetimes <- sample_components(series, n, par = theta)
    generate_masked_series_data(comp_lifetimes, n, m, tau, p,
                                defaults$lifetime, defaults$omega,
                                defaults$candset)
  }
}
