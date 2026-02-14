# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Decode candidate set matrix from Boolean columns
#'
#' Extracts Boolean columns matching the pattern `prefix + digits` from a
#' data frame and assembles them into a logical matrix. This replaces
#' `md.tools::md_decode_matrix()` to avoid the dependency.
#'
#' @param df data frame containing candidate set columns
#' @param prefix column name prefix (default `"x"`)
#' @return logical matrix with one column per component, or NULL if no
#'   matching columns found
#' @keywords internal
decode_candidate_matrix <- function(df, prefix = "x") {
  cols <- grep(paste0("^", prefix, "\\d+$"), names(df), value = TRUE)
  if (length(cols) == 0) return(NULL)
  cols <- cols[order(as.integer(sub(prefix, "", cols)))]
  mat <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(mat) <- "logical"
  mat
}

#' Extract model column name defaults
#'
#' @param model likelihood model object
#' @return list with lifetime, lifetime_upper, omega, candset defaults
#' @keywords internal
extract_model_defaults <- function(model) {
  list(
    lifetime = model$lifetime %||% "t",
    lifetime_upper = model$lifetime_upper %||% "t_upper",
    omega = model$omega %||% "omega",
    candset = model$candset %||% "x"
  )
}

#' Extract and validate masked data from a data frame
#'
#' Shared validation logic for all likelihood model methods. Checks that the
#' data frame is non-empty, required columns exist, decodes the candidate set
#' matrix, and validates observation types.
#'
#' @param df masked data frame
#' @param lifetime column name for system lifetime
#' @param omega column name for observation type
#' @param candset column prefix for candidate set indicators
#' @param lifetime_upper column name for interval upper bound
#' @return list with components: t, t_upper, omega, C, m, n
#' @keywords internal
extract_md_data <- function(df, lifetime, omega, candset,
                            lifetime_upper = NULL) {
  n <- nrow(df)
  if (n == 0) stop("df is empty")
  if (!lifetime %in% colnames(df))
    stop("lifetime variable not in colnames(df)")
  if (!omega %in% colnames(df))
    stop("omega variable '", omega, "' not in colnames(df)")

  cmat <- decode_candidate_matrix(df, candset)
  if (is.null(cmat) || ncol(cmat) == 0)
    stop("no candidate set found for candset prefix '", candset, "'")
  m <- ncol(cmat)

  omega_vals <- as.character(df[[omega]])
  valid_types <- c("exact", "right", "left", "interval")
  invalid <- setdiff(unique(omega_vals), valid_types)
  if (length(invalid) > 0) {
    stop(
      "invalid omega values: ", paste(invalid, collapse = ", "),
      ". Must be one of: ", paste(valid_types, collapse = ", ")
    )
  }

  t_upper <- NULL
  if (!is.null(lifetime_upper) && lifetime_upper %in% colnames(df))
    t_upper <- df[[lifetime_upper]]

  # Validate observations
  for (i in seq_len(n)) {
    has_cand <- any(cmat[i, ])
    if (omega_vals[i] == "exact" && !has_cand)
      stop("C1 violated: exact observation with empty candidate set at row ", i)
    if (omega_vals[i] == "left" && !has_cand)
      stop("left-censored observation must have non-empty candidate set at row ", i)
    if (omega_vals[i] == "interval") {
      if (!has_cand)
        stop("interval-censored observation must have non-empty candidate set at row ", i)
      if (is.null(t_upper))
        stop("interval-censored observations require a '",
             lifetime_upper %||% "t_upper", "' column")
      if (t_upper[i] <= df[[lifetime]][i])
        stop("interval-censored observation requires t_upper > t at row ", i)
    }
  }

  list(t = df[[lifetime]], omega = omega_vals, C = cmat, m = m, n = n,
       t_upper = t_upper)
}

#' Generate masked series system data
#'
#' Creates masked data from pre-generated component lifetimes. Applies system
#' lifetime calculation (minimum of components), right-censoring at tau, and
#' candidate set generation satisfying C1-C2-C3.
#'
#' @param comp_lifetimes n x m matrix of component lifetimes
#' @param n number of observations
#' @param m number of components
#' @param tau right-censoring time
#' @param p masking probability for non-failed components
#' @param default_lifetime column name for system lifetime
#' @param default_omega column name for observation type
#' @param default_candset column prefix for candidate sets
#' @return data frame with system lifetime, observation type, and candidate sets
#' @importFrom stats runif
#' @keywords internal
generate_masked_series_data <- function(comp_lifetimes, n, m, tau, p,
                                        default_lifetime, default_omega,
                                        default_candset) {
  sys_lifetime <- apply(comp_lifetimes, 1, min)
  failed_comp <- apply(comp_lifetimes, 1, which.min)

  is_exact <- sys_lifetime <= tau
  sys_lifetime <- pmin(sys_lifetime, tau)

  omega_vals <- ifelse(is_exact, "exact", "right")

  candset <- matrix(FALSE, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    if (is_exact[i]) {
      candset[i, failed_comp[i]] <- TRUE
      if (p > 0 && m > 1) {
        others <- seq_len(m)[-failed_comp[i]]
        candset[i, others] <- runif(length(others)) < p
      }
    }
  }

  df <- data.frame(sys_lifetime, omega_vals, stringsAsFactors = FALSE)
  names(df) <- c(default_lifetime, default_omega)

  for (j in seq_len(m)) {
    df[[paste0(default_candset, j)]] <- candset[, j]
  }

  df
}
