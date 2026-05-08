## test-r460-diagnosis.R
##
## Purpose: Trace the EXACT failure chain causing:
##   Error in `names(x) <- value`:
##     'names' attribute [12] must be the same length as the vector [11]
##
## Hypothesis being tested:
##   In R 4.6.0, dplyr::case_when() with scalar LHS + vector RHS inside
##   is_ss() causes `doses_required` to become all-NA, which propagates
##   through calculate_cl() → run_single_point() causing cl/vd = NA,
##   which then causes eval_perf_1cmpt() to return rep(NA, 5) instead of
##   rep(NA, 6), making results_df have 11 cols instead of 12,
##   and colnames(base.out) <- c(12 names) crashes.
##
## Each test isolates one link in the chain.
## Run via: testthat::test_file("tests/testthat/test-r460-diagnosis.R")

library(testthat)
library(nlmixr2autoinit)
library(nlmixr2data)

suppress_all <- function(expr) suppressWarnings(suppressMessages(expr))

# =============================================================================
# STEP 0: Session info — always printed, never fails
# =============================================================================
test_that("DIAG-0: print session info for CI logs", {
  cat("\n========== Diagnostic Session Info ==========\n")
  cat("R version      :", R.version$version.string, "\n")
  cat("dplyr version  :", as.character(packageVersion("dplyr")), "\n")
  cat("nlmixr2autoinit:", as.character(packageVersion("nlmixr2autoinit")), "\n")
  cat("=============================================\n")
  expect_true(TRUE)
})

# =============================================================================
# STEP 1: Isolate case_when() scalar-LHS behaviour in this R/dplyr version
#
# Hypothesis: in dplyr >= some version under R 4.6.0,
#   case_when(scalar ~ vector) returns all-NA instead of recycling scalar.
# =============================================================================
test_that("DIAG-1: case_when scalar LHS recycles correctly into a column", {
  ss_method <- "combined"   # length-1 scalar — exactly what is_ss() uses

  df <- data.frame(
    doses_req_fixed = c(5, 5, 5),
    doses_req_hl    = c(3, 4, 6)
  )

  result <- df |>
    dplyr::mutate(
      doses_required = dplyr::case_when(
        ss_method == "combined"       ~ pmin(doses_req_fixed, doses_req_hl, na.rm = TRUE),
        ss_method == "fixed_doses"    ~ doses_req_fixed,
        ss_method == "half_life_based" ~ doses_req_hl
      )
    )

  cat("\nDIAG-1 doses_required values:", result$doses_required, "\n")
  cat("DIAG-1 any NA?              :", any(is.na(result$doses_required)), "\n")

  # If this fails on R 4.6.0, case_when is the root cause
  expect_false(
    any(is.na(result$doses_required)),
    label = paste0(
      "case_when(scalar ~ vector) must not produce NA under dplyr ",
      packageVersion("dplyr"), " / R ", R.version$major, ".", R.version$minor
    )
  )
})

# =============================================================================
# STEP 2: Verify is_ss() doses_required column is never all-NA
#
# If STEP 1 fails, this will also fail and confirm the propagation.
# =============================================================================
test_that("DIAG-2: is_ss() doses_required column contains no all-NA rows for any subject", {
  dat <- suppress_all({
    processData(Bolus_1CPT)$dat
  })

  out <- suppress_all(
    is_ss(df = dat, half_life = 11.26, ssctrl = ss_control())
  )

  cat("\nDIAG-2 doses_required summary:\n")
  print(summary(out$doses_required))
  cat("DIAG-2 NA count:", sum(is.na(out$doses_required)), "/", nrow(out), "\n")

  expect_false(
    all(is.na(out$doses_required)),
    label = "is_ss() doses_required must not be entirely NA — scalar case_when may have failed"
  )
})

# =============================================================================
# STEP 3: Verify is_ss() SteadyState column is well-formed
#
# If doses_required is all-NA, SteadyState will also be wrong.
# =============================================================================
test_that("DIAG-3: is_ss() SteadyState column is logical with no unexpected all-FALSE pattern", {
  dat <- suppress_all({
    processData(Bolus_1CPT)$dat
  })

  out <- suppress_all(
    is_ss(df = dat, half_life = 11.26, ssctrl = ss_control())
  )

  cat("\nDIAG-3 SteadyState table:\n")
  print(table(out$SteadyState, useNA = "always"))

  expect_true(
    is.logical(out$SteadyState),
    label = "SteadyState column must be logical"
  )
  expect_false(
    all(is.na(out$SteadyState)),
    label = "SteadyState must not be entirely NA"
  )
})

# =============================================================================
# STEP 4: Verify eval_perf_1cmpt() always returns length-6 vector
#
# The function has 3 early-return paths that return rep(NA, 5) — wrong length.
# This test checks ALL return paths.
# =============================================================================
test_that("DIAG-4: eval_perf_1cmpt() always returns a length-6 vector", {
  dat <- Bolus_1CPT

  # Path A: normal valid parameters
  res_valid <- suppress_all(
    eval_perf_1cmpt(dat = dat, est.method = "rxSolve",
                    ka = NA, cl = 4, vd = 70, route = "bolus")
  )
  cat("\nDIAG-4A (valid bolus)  length:", length(res_valid),
      " values:", res_valid, "\n")
  expect_equal(length(res_valid), 6,
               label = "eval_perf_1cmpt() valid bolus path must return length 6")

  # Path B: NA cl triggers defensive return
  res_na_cl <- suppress_all(
    eval_perf_1cmpt(dat = dat, est.method = "rxSolve",
                    ka = NA, cl = NA, vd = 70, route = "bolus")
  )
  cat("DIAG-4B (NA cl)        length:", length(res_na_cl),
      " values:", res_na_cl, "\n")
  expect_equal(length(res_na_cl), 6,
               label = "eval_perf_1cmpt() NA-cl defensive path must return length 6 — currently returns 5!")

  # Path C: cl <= 0 triggers defensive return
  res_neg_cl <- suppress_all(
    eval_perf_1cmpt(dat = dat, est.method = "rxSolve",
                    ka = NA, cl = -1, vd = 70, route = "bolus")
  )
  cat("DIAG-4C (cl <= 0)      length:", length(res_neg_cl),
      " values:", res_neg_cl, "\n")
  expect_equal(length(res_neg_cl), 6,
               label = "eval_perf_1cmpt() negative-cl defensive path must return length 6 — currently returns 5!")

  # Path D: oral with NA ka
  res_na_ka <- suppress_all(
    eval_perf_1cmpt(dat = Oral_1CPT, est.method = "rxSolve",
                    ka = NA, cl = 4, vd = 70, route = "oral")
  )
  cat("DIAG-4D (oral NA ka)   length:", length(res_na_ka),
      " values:", res_na_ka, "\n")
  expect_equal(length(res_na_ka), 6,
               label = "eval_perf_1cmpt() oral-NA-ka defensive path must return length 6 — currently returns 5!")
})

# =============================================================================
# STEP 5: Verify hybrid_eval_perf_1cmpt() returns a data.frame with 12 cols
#
# If eval_perf_1cmpt() returns length 5 for any row, do.call(rbind) will
# produce 11 cols instead of 12, and colnames<- will crash.
# =============================================================================
test_that("DIAG-5: hybrid_eval_perf_1cmpt() returns data.frame with exactly 12 columns", {
  result <- suppress_all(
    hybrid_eval_perf_1cmpt(
      route       = "bolus",
      dat         = Bolus_1CPT,
      sp_out_ka   = NA,  sp_out_cl  = 4.0, sp_out_vd  = 70,
      graph_out_ka= NA,  graph_out_cl= 3.8, graph_out_vd= 68,
      nca_fd_ka   = NA,  nca_fd_cl  = 4.2, nca_fd_vd  = 72,
      nca_efd_ka  = NA,  nca_efd_cl = NA,  nca_efd_vd = NA,
      nca_all_ka  = NA,  nca_all_cl = 4.1, nca_all_vd = 71,
      verbose     = FALSE
    )
  )

  cat("\nDIAG-5 hybrid_eval_perf_1cmpt() ncol:", ncol(result), "\n")
  cat("DIAG-5 colnames:", paste(names(result), collapse = ", "), "\n")

  expect_equal(
    ncol(result), 12,
    label = paste0(
      "hybrid_eval_perf_1cmpt() must return 12 cols ",
      "(ka_source, cl_source, vd_source, ka_value, cl_value, vd_value, ",
      "APE, MAE, MAPE, RMSE, rRMSE1, rRMSE2). ",
      "Got ", ncol(result), " — mismatch caused by eval_perf_1cmpt() ",
      "returning rep(NA,5) instead of rep(NA,6) on defensive paths."
    )
  )
})

# =============================================================================
# STEP 6: Verify run_single_point() returns non-NA cl and vd
#
# If is_ss() is broken, calculate_cl() may produce NA cl/vd,
# which then flows into eval_perf_1cmpt() and triggers the defensive return.
# =============================================================================
test_that("DIAG-6: run_single_point() returns non-NA cl and vd for Bolus_1CPT", {
  dat <- suppress_all(processData(Bolus_1CPT)$dat)

  sp_result <- suppress_all(
    run_single_point(
      dat        = dat,
      route      = "bolus",
      half_life  = 11.26,
      dose_type  = "combined_doses",
      pooled_ctrl= pooled_control(),
      ssctrl     = ss_control()
    )
  )

  cat("\nDIAG-6 single-point cl :", sp_result$singlepoint.results$cl, "\n")
  cat("DIAG-6 single-point vd :", sp_result$singlepoint.results$vd, "\n")

  expect_false(
    is.na(sp_result$singlepoint.results$cl),
    label = "run_single_point() cl must not be NA — if NA, is_ss() case_when is broken"
  )
  expect_false(
    is.na(sp_result$singlepoint.results$vd),
    label = "run_single_point() vd must not be NA — if NA, is_ss() case_when is broken"
  )
})
