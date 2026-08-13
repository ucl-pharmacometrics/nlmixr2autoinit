# ---- Sparse-data guard ----

test_that("getPPKinits stops cleanly when only one observation remains after filtering", {
  # Reproducer from GH-2 follow-up: subject 1 has no observations and is
  # filtered out; subject 2's only "observation" is at TIME = 0 (dose time),
  # leaving a single row that nothing downstream can fit.
  d <- data.frame(
    ID   = c(1, 2, 2),
    EVID = c(1, 1, 0),
    CMT  = c(1, 1, 0),
    AMT  = c(0, 1, 0),
    TIME = 0,
    DV   = c(NA, NA, 1)
  )

  expect_error(
    suppressMessages(getPPKinits(d, verbose = FALSE)),
    "At least 2 observations are required"
  )
})

test_that("getPPKinits stops cleanly when zero observations remain (entire dataset is doses)", {
  d_all_doses <- data.frame(
    ID   = c(1, 1, 2),
    EVID = c(1, 1, 1),
    CMT  = c(1, 1, 1),
    AMT  = c(100, 100, 100),
    TIME = c(0, 24, 0),
    DV   = c(NA, NA, NA)
  )

  # processData removes the all-dose subjects first and errors before
  # getPPKinits sees them.
  expect_error(
    suppressMessages(getPPKinits(d_all_doses, verbose = FALSE)),
    "No subjects with observations"
  )
})

# ---- Defensive guard inside hybrid_eval_perf_1cmpt ----

test_that("hybrid_eval_perf_1cmpt errors clearly when no valid parameter sources exist", {
  # All CL/Vd sources are NA - simulates the deeply-sparse case that
  # previously produced NULL and triggered `colnames<-` on a non-matrix.
  expect_error(
    hybrid_eval_perf_1cmpt(
      route        = "oral",
      dat          = data.frame(),
      sp_out_ka    = NA, sp_out_cl    = NA, sp_out_vd    = NA,
      graph_out_ka = NA, graph_out_cl = NA, graph_out_vd = NA,
      nca_fd_ka    = NA, nca_fd_cl    = NA, nca_fd_vd    = NA,
      nca_efd_ka   = NA, nca_efd_cl   = NA, nca_efd_vd   = NA,
      nca_all_ka   = NA, nca_all_cl   = NA, nca_all_vd   = NA,
      verbose      = FALSE
    ),
    "no valid parameter sources"
  )
})
