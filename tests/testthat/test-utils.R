test_that("decode_candidate_matrix extracts Boolean columns", {
  df <- data.frame(t = 1, omega = "exact", x1 = TRUE, x2 = FALSE, x3 = TRUE)
  mat <- dfr.lik.series.md:::decode_candidate_matrix(df, "x")
  expect_true(is.logical(mat))
  expect_equal(ncol(mat), 3)
  expect_equal(mat[1, ], c(x1 = TRUE, x2 = FALSE, x3 = TRUE))
})

test_that("decode_candidate_matrix orders columns numerically", {
  df <- data.frame(t = 1, x10 = TRUE, x2 = FALSE, x1 = TRUE)
  mat <- dfr.lik.series.md:::decode_candidate_matrix(df, "x")
  expect_equal(colnames(mat), c("x1", "x2", "x10"))
})

test_that("decode_candidate_matrix returns NULL for no matches", {
  df <- data.frame(t = 1, y1 = TRUE, y2 = FALSE)
  expect_null(dfr.lik.series.md:::decode_candidate_matrix(df, "x"))
})

test_that("decode_candidate_matrix handles custom prefix", {
  df <- data.frame(t = 1, c1 = TRUE, c2 = FALSE)
  mat <- dfr.lik.series.md:::decode_candidate_matrix(df, "c")
  expect_equal(ncol(mat), 2)
})

test_that("extract_md_data validates required columns", {
  df <- data.frame(t = 1, omega = "exact", x1 = TRUE)
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "time", "omega", "x"),
    "lifetime variable not in colnames"
  )
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "t", "status", "x"),
    "omega variable"
  )
})

test_that("extract_md_data validates empty data frame", {
  df <- data.frame(t = numeric(0), omega = character(0), x1 = logical(0))
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "t", "omega", "x"),
    "df is empty"
  )
})

test_that("extract_md_data validates omega values", {
  df <- data.frame(t = 1, omega = "invalid", x1 = TRUE)
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "t", "omega", "x"),
    "invalid omega values"
  )
})

test_that("extract_md_data validates C1 for exact observations", {
  df <- data.frame(t = 1, omega = "exact", x1 = FALSE, x2 = FALSE)
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "t", "omega", "x"),
    "C1 violated"
  )
})

test_that("extract_md_data validates interval-censored requirements", {
  # Missing t_upper column
  df <- data.frame(t = 1, omega = "interval", x1 = TRUE)
  expect_error(
    dfr.lik.series.md:::extract_md_data(df, "t", "omega", "x", "t_upper"),
    "interval-censored observations require"
  )

  # t_upper <= t
  df2 <- data.frame(t = 5, t_upper = 3, omega = "interval", x1 = TRUE)
  expect_error(
    dfr.lik.series.md:::extract_md_data(df2, "t", "omega", "x", "t_upper"),
    "t_upper > t"
  )
})

test_that("extract_md_data returns correct structure", {
  df <- data.frame(
    t = c(1, 2, 3),
    omega = c("exact", "right", "exact"),
    x1 = c(TRUE, FALSE, FALSE),
    x2 = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  result <- dfr.lik.series.md:::extract_md_data(df, "t", "omega", "x")
  expect_equal(result$n, 3)
  expect_equal(result$m, 2)
  expect_equal(result$t, c(1, 2, 3))
  expect_equal(result$omega, c("exact", "right", "exact"))
  expect_true(is.logical(result$C))
})

test_that("generate_masked_series_data produces valid output", {
  set.seed(42)
  comp_lifetimes <- matrix(rexp(300, rate = 0.1), nrow = 100, ncol = 3)
  df <- dfr.lik.series.md:::generate_masked_series_data(
    comp_lifetimes, 100, 3, Inf, 0.5, "t", "omega", "x"
  )
  expect_equal(nrow(df), 100)
  expect_true(all(df$omega == "exact")) # tau = Inf, all exact
  # At least some observations should have multiple candidates (p = 0.5)
  cand_sizes <- df$x1 + df$x2 + df$x3
  expect_true(any(cand_sizes > 1))
  # All exact obs have at least one candidate
  expect_true(all(cand_sizes >= 1))
})
