# Focused regression tests for sim_sens_*() resilience against Fit_*() errors.
#
# Background: each sim_sens_*() iterates over a parameter grid with
# purrr::pmap_dfr and calls Fit_*() inside suppressMessages/suppressWarnings.
# Previously, when Fit_*() errored (e.g. lotri 1.1+ rejecting odd inits with
# "subscript out of bounds"), the error propagated up and aborted the entire
# parameter sweep — which in turn aborted getPPKinits(). The fix wraps the
# Fit_*() call in tryCatch and substitutes rep(NA_real_, 6) for the metrics
# when Fit_*() fails. The downstream selection logic already uses
# min(..., na.rm = TRUE), so NA metrics flow through correctly.

# Minimal observation/dose dataset shared by all tests. Two observations are
# enough for metrics.() to compute non-degenerate values.
make_dat <- function() {
  data.frame(
    ID   = c(1L, 1L, 1L),
    TIME = c(0, 1, 2),
    EVID = c(1L, 0L, 0L),
    AMT  = c(100, 0, 0),
    DV   = c(0, 10, 5),
    CMT  = c(1L, 2L, 2L)
  )
}

# A Fit_*() stand-in that returns a $cp vector of length matching the
# observations in `dat`, so metrics.() succeeds and yields finite numbers.
fake_fit <- function(...) {
  list(cp = c(9, 4.5))
}

err_fit <- function(...) {
  stop("forced error: simulating lotri 'subscript out of bounds' failure")
}

# --- sim_sens_1cmpt_mm -----------------------------------------------------

test_that("sim_sens_1cmpt_mm returns NA metrics when Fit_1cmpt_mm_oral errors", {
  local_mocked_bindings(Fit_1cmpt_mm_oral = err_fit)
  res <- sim_sens_1cmpt_mm(
    dat = make_dat(),
    sim_vmax = list(mode = "manual", values = 1000),
    sim_km   = list(mode = "manual", values = 250),
    sim_vd   = list(mode = "manual", values = 70),
    sim_ka   = list(mode = "manual", values = 1),
    route    = "oral",
    verbose  = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
  expect_equal(res$Vmax, 1000)
  expect_equal(res$Km, 250)
  expect_equal(res$Vd, 70)
  expect_equal(res$Ka, 1)
})

test_that("sim_sens_1cmpt_mm returns NA metrics when Fit_1cmpt_mm_iv errors", {
  local_mocked_bindings(Fit_1cmpt_mm_iv = err_fit)
  res <- sim_sens_1cmpt_mm(
    dat = make_dat(),
    sim_vmax = list(mode = "manual", values = 1000),
    sim_km   = list(mode = "manual", values = 250),
    sim_vd   = list(mode = "manual", values = 70),
    sim_ka   = list(mode = "manual", values = NA),
    route    = "iv",
    verbose  = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
})

test_that("sim_sens_1cmpt_mm returns numeric metrics when Fit_1cmpt_mm_oral succeeds", {
  local_mocked_bindings(Fit_1cmpt_mm_oral = fake_fit)
  res <- sim_sens_1cmpt_mm(
    dat = make_dat(),
    sim_vmax = list(mode = "manual", values = 1000),
    sim_km   = list(mode = "manual", values = 250),
    sim_vd   = list(mode = "manual", values = 70),
    sim_ka   = list(mode = "manual", values = 1),
    route    = "oral",
    verbose  = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.finite(unlist(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")]))))
})

# --- sim_sens_2cmpt --------------------------------------------------------

test_that("sim_sens_2cmpt returns NA metrics when Fit_2cmpt_oral errors", {
  local_mocked_bindings(Fit_2cmpt_oral = err_fit)
  res <- sim_sens_2cmpt(
    dat    = make_dat(),
    sim_ka = list(mode = "manual", values = 1),
    sim_cl = list(mode = "manual", values = 4),
    sim_vc = list(mode = "manual", values = 50),
    sim_vp = list(mode = "manual", values = 50),
    sim_q  = list(mode = "manual", values = 4),
    route  = "oral",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
})

test_that("sim_sens_2cmpt returns NA metrics when Fit_2cmpt_iv errors", {
  local_mocked_bindings(Fit_2cmpt_iv = err_fit)
  res <- sim_sens_2cmpt(
    dat    = make_dat(),
    sim_ka = list(mode = "manual", values = NA),
    sim_cl = list(mode = "manual", values = 4),
    sim_vc = list(mode = "manual", values = 50),
    sim_vp = list(mode = "manual", values = 50),
    sim_q  = list(mode = "manual", values = 4),
    route  = "iv",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
})

test_that("sim_sens_2cmpt returns numeric metrics when Fit_2cmpt_oral succeeds", {
  local_mocked_bindings(Fit_2cmpt_oral = fake_fit)
  res <- sim_sens_2cmpt(
    dat    = make_dat(),
    sim_ka = list(mode = "manual", values = 1),
    sim_cl = list(mode = "manual", values = 4),
    sim_vc = list(mode = "manual", values = 50),
    sim_vp = list(mode = "manual", values = 50),
    sim_q  = list(mode = "manual", values = 4),
    route  = "oral",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.finite(unlist(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")]))))
})

# --- sim_sens_3cmpt --------------------------------------------------------

test_that("sim_sens_3cmpt returns NA metrics when Fit_3cmpt_oral errors", {
  local_mocked_bindings(Fit_3cmpt_oral = err_fit)
  res <- sim_sens_3cmpt(
    dat     = make_dat(),
    sim_vc  = list(mode = "manual", values = 50),
    sim_vp  = list(mode = "manual", values = 50),
    sim_vp2 = list(mode = "manual", values = 50),
    sim_q   = list(mode = "manual", values = 4),
    sim_q2  = list(mode = "manual", values = 4),
    sim_cl  = list(mode = "manual", values = 4),
    sim_ka  = list(mode = "manual", values = 1),
    route   = "oral",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
})

test_that("sim_sens_3cmpt returns NA metrics when Fit_3cmpt_iv errors", {
  local_mocked_bindings(Fit_3cmpt_iv = err_fit)
  res <- sim_sens_3cmpt(
    dat     = make_dat(),
    sim_vc  = list(mode = "manual", values = 50),
    sim_vp  = list(mode = "manual", values = 50),
    sim_vp2 = list(mode = "manual", values = 50),
    sim_q   = list(mode = "manual", values = 4),
    sim_q2  = list(mode = "manual", values = 4),
    sim_cl  = list(mode = "manual", values = 4),
    sim_ka  = list(mode = "manual", values = NA),
    route   = "iv",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.na(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")])))
})

test_that("sim_sens_3cmpt returns numeric metrics when Fit_3cmpt_oral succeeds", {
  local_mocked_bindings(Fit_3cmpt_oral = fake_fit)
  res <- sim_sens_3cmpt(
    dat     = make_dat(),
    sim_vc  = list(mode = "manual", values = 50),
    sim_vp  = list(mode = "manual", values = 50),
    sim_vp2 = list(mode = "manual", values = 50),
    sim_q   = list(mode = "manual", values = 4),
    sim_q2  = list(mode = "manual", values = 4),
    sim_cl  = list(mode = "manual", values = 4),
    sim_ka  = list(mode = "manual", values = 1),
    route   = "oral",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true(all(is.finite(unlist(res[, c("APE", "MAE", "MAPE", "RMSE", "rRMSE1", "rRMSE2")]))))
})
