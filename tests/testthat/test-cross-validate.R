test_that("loglik matches exp_series_md_c1_c2_c3 on exact data", {
  skip_if_not_installed("maskedcauses")

  par <- c(0.1, 0.2, 0.3)

  # Generate data using the reference model
  ref_model <- maskedcauses::exp_series_md_c1_c2_c3()
  ref_rdata <- likelihood.model::rdata(ref_model)
  set.seed(42)
  df <- ref_rdata(theta = par, n = 200, tau = Inf, p = 0.5)

  # Evaluate loglik using reference (closed-form)
  ref_ll <- likelihood.model::loglik(ref_model)
  ll_ref <- ref_ll(df, par)

  # Evaluate loglik using our package (numerical)
  our_model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L, 1L))
  our_ll <- loglik(our_model)
  ll_ours <- our_ll(df, par)

  expect_equal(ll_ours, ll_ref, tolerance = 1e-6)
})

test_that("loglik matches exp_series_md_c1_c2_c3 on censored data", {
  skip_if_not_installed("maskedcauses")

  par <- c(0.1, 0.2, 0.3)

  ref_model <- maskedcauses::exp_series_md_c1_c2_c3()
  ref_rdata <- likelihood.model::rdata(ref_model)
  set.seed(123)
  df <- ref_rdata(theta = par, n = 300, tau = 5, p = 0.3)

  ref_ll <- likelihood.model::loglik(ref_model)
  ll_ref <- ref_ll(df, par)

  our_model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L, 1L))
  our_ll <- loglik(our_model)
  ll_ours <- our_ll(df, par)

  expect_equal(ll_ours, ll_ref, tolerance = 1e-6)
})

test_that("MLE point estimates match reference on same dataset", {
  skip_if_not_installed("maskedcauses")

  par <- c(0.1, 0.2, 0.3)

  ref_model <- maskedcauses::exp_series_md_c1_c2_c3()
  ref_rdata <- likelihood.model::rdata(ref_model)
  set.seed(42)
  df <- ref_rdata(theta = par, n = 500, tau = Inf, p = 0.5)

  # Fit with reference
  ref_solver <- generics::fit(ref_model)
  ref_result <- suppressWarnings(ref_solver(df, par = c(0.5, 0.5, 0.5)))

  # Fit with ours
  our_model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L, 1L))
  our_solver <- fit(our_model)
  our_result <- suppressWarnings(our_solver(df, par = c(0.5, 0.5, 0.5)))

  # Sum of rates should match (exponential series are identifiable in sum)
  expect_equal(sum(coef(our_result)), sum(coef(ref_result)), tolerance = 0.05)

  # Log-likelihoods at MLE should be very close
  expect_equal(our_result$loglik, ref_result$loglik, tolerance = 0.5)
})

test_that("loglik matches at multiple parameter values", {
  skip_if_not_installed("maskedcauses")

  # Generate shared dataset
  ref_model <- maskedcauses::exp_series_md_c1_c2_c3()
  ref_rdata <- likelihood.model::rdata(ref_model)
  set.seed(77)
  df <- ref_rdata(theta = c(0.2, 0.3, 0.5), n = 100, tau = Inf, p = 0.4)

  ref_ll <- likelihood.model::loglik(ref_model)

  our_model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L, 1L))
  our_ll <- loglik(our_model)

  # Check at several parameter values
  test_pars <- list(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.1, 0.1),
    c(0.5, 0.5, 0.5),
    c(1.0, 2.0, 3.0)
  )

  for (p in test_pars) {
    expect_equal(our_ll(df, p), ref_ll(df, p), tolerance = 1e-6,
                 info = paste("par =", paste(p, collapse = ", ")))
  }
})

test_that("loglik matches 2-component exponential reference", {
  skip_if_not_installed("maskedcauses")

  par <- c(0.5, 1.0)

  ref_model <- maskedcauses::exp_series_md_c1_c2_c3()
  ref_rdata <- likelihood.model::rdata(ref_model)
  set.seed(99)
  df <- ref_rdata(theta = par, n = 150, tau = 3, p = 0.5)

  ref_ll <- likelihood.model::loglik(ref_model)
  ll_ref <- ref_ll(df, par)

  our_model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L))
  our_ll <- loglik(our_model)
  ll_ours <- our_ll(df, par)

  expect_equal(ll_ours, ll_ref, tolerance = 1e-6)
})
