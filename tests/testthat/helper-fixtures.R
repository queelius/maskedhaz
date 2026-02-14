# Shared test helpers for dfr.lik.series.md
library(dfr.dist)
library(dfr.dist.series)

#' Create a 3-component exponential series model
make_exp_model <- function(rates = c(0.1, 0.2, 0.3)) {
  comps <- lapply(rates, dfr_exponential)
  dfr_series_md(components = comps)
}

#' Create a 2-component Weibull series model
make_weibull_model <- function(shapes = c(2, 1.5), scales = c(100, 200)) {
  comps <- mapply(dfr_weibull, shape = shapes, scale = scales,
                  SIMPLIFY = FALSE)
  dfr_series_md(components = comps)
}

#' Create a mixed Weibull + exponential series model
make_mixed_model <- function() {
  comps <- list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
  )
  dfr_series_md(components = comps)
}

#' Generate exact-only (no censoring) masked data for exponential series
make_exp_exact_data <- function(n = 100, rates = c(0.1, 0.2, 0.3),
                                p = 0.5, seed = 42) {
  set.seed(seed)
  model <- make_exp_model(rates)
  rdata_fn <- rdata(model)
  rdata_fn(theta = rates, n = n, tau = Inf, p = p)
}

#' Generate censored masked data for exponential series
make_exp_censored_data <- function(n = 200, rates = c(0.1, 0.2, 0.3),
                                   tau = 5, p = 0.5, seed = 42) {
  set.seed(seed)
  model <- make_exp_model(rates)
  rdata_fn <- rdata(model)
  rdata_fn(theta = rates, n = n, tau = tau, p = p)
}

#' Analytical log-likelihood for exponential series (exact + right only)
#' Used as cross-validation reference
exp_series_loglik_analytical <- function(df, par) {
  m <- length(par)
  lambda_sys <- sum(par)

  cmat_cols <- grep("^x\\d+$", names(df), value = TRUE)
  cmat_cols <- cmat_cols[order(as.integer(sub("x", "", cmat_cols)))]
  cmat <- as.matrix(df[, cmat_cols, drop = FALSE])
  storage.mode(cmat) <- "logical"

  ll <- 0
  for (i in seq_len(nrow(df))) {
    if (df$omega[i] %in% c("exact", "right")) {
      ll <- ll - lambda_sys * df$t[i]
    }
    if (df$omega[i] == "exact") {
      lambda_c <- sum(par[cmat[i, ]])
      ll <- ll + log(lambda_c)
    }
    if (df$omega[i] == "left") {
      lambda_c <- sum(par[cmat[i, ]])
      ll <- ll + log(lambda_c) + log(-expm1(-lambda_sys * df$t[i])) -
        log(lambda_sys)
    }
  }
  ll
}
