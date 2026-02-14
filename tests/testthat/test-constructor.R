test_that("dfr_series_md builds from components list", {
  model <- dfr_series_md(components = list(
    dfr_exponential(0.1),
    dfr_exponential(0.2),
    dfr_exponential(0.3)
  ))
  expect_s3_class(model, "dfr_series_md")
  expect_s3_class(model, "series_md")
  expect_s3_class(model, "likelihood_model")
  expect_true(is_dfr_series_md(model))
})

test_that("dfr_series_md builds from pre-built dfr_dist_series", {
  sys <- dfr_dist_series(list(
    dfr_weibull(shape = 2, scale = 100),
    dfr_exponential(0.05)
  ))
  model <- dfr_series_md(series = sys)
  expect_s3_class(model, "dfr_series_md")
  expect_true(is_dfr_dist_series(model$series))
})

test_that("dfr_series_md stores column names", {
  model <- make_exp_model()
  expect_equal(model$lifetime, "t")
  expect_equal(model$lifetime_upper, "t_upper")
  expect_equal(model$omega, "omega")
  expect_equal(model$candset, "x")
})

test_that("dfr_series_md custom column names", {
  model <- dfr_series_md(
    components = list(dfr_exponential(0.1), dfr_exponential(0.2)),
    lifetime = "time", omega = "status", candset = "c"
  )
  expect_equal(model$lifetime, "time")
  expect_equal(model$omega, "status")
  expect_equal(model$candset, "c")
})

test_that("dfr_series_md errors on invalid input", {
  expect_error(dfr_series_md(), "Either 'series' or 'components' must be provided")
  expect_error(dfr_series_md(series = "not a series"),
               "'series' must be a dfr_dist_series object")
})

test_that("is_dfr_series_md returns FALSE for non-model objects", {
  expect_false(is_dfr_series_md(42))
  expect_false(is_dfr_series_md(NULL))
  expect_false(is_dfr_series_md(dfr_exponential(0.1)))
})

test_that("print.dfr_series_md works", {
  model <- make_exp_model()
  expect_output(print(model), "Masked-cause likelihood model")
  expect_output(print(model), "3-component series")
})

test_that("series metadata is accessible", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  expect_equal(model$series$m, 3)
  expect_equal(model$series$n_par, c(1L, 1L, 1L))
  expect_equal(model$series$par, c(0.1, 0.2, 0.3))
})

test_that("print works with unknown params", {
  model <- dfr_series_md(components = list(
    dfr_exponential(), dfr_exponential()
  ), n_par = c(1L, 1L))
  expect_output(print(model), "unknown")
})
