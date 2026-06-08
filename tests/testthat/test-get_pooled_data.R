result_bolus  <- processData(Bolus_1CPT, verbose = FALSE)
result_oral   <- processData(Oral_1CPT,  verbose = FALSE)
ctrl          <- initsControl()

dose_type_bolus <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Type"]
dose_type_oral  <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Type"]

pooled_bolus <- get_pooled_data(result_bolus$dat, dose_type_bolus, ctrl$pooled.control)
pooled_oral  <- get_pooled_data(result_oral$dat,  dose_type_oral,  ctrl$pooled.control)

test_that("get_pooled_data returns a list for bolus", {
  expect_type(pooled_bolus, "list")
})

test_that("get_pooled_data contains datpooled_fd element", {
  expect_true("datpooled_fd" %in% names(pooled_bolus))
})

test_that("get_pooled_data datpooled_fd contains binned.df", {
  expect_true("binned.df" %in% names(pooled_bolus$datpooled_fd))
})

test_that("get_pooled_data binned.df is a data frame", {
  expect_s3_class(pooled_bolus$datpooled_fd$binned.df, "data.frame")
})

test_that("get_pooled_data returns a list for oral", {
  expect_type(pooled_oral, "list")
})

# ---- Empty most_common_ii guard (GH-2) ----

# A two-subject dataset with one dose per subject and several observation rows
# per subject. With only one dose per ID, `diff(TIME)` within each subject is
# empty, so `most_common_ii` becomes numeric(0).
d_single_dose <- data.frame(
  ID   = rep(1:2, each = 5),
  EVID = rep(c(1, 0, 0, 0, 0), 2),
  CMT  = 1,
  AMT  = rep(c(100, 0, 0, 0, 0), 2),
  TIME = rep(c(0, 1, 2, 4, 8), 2),
  DV   = rep(c(NA, 0.9, 0.7, 0.5, 0.2), 2)
)
processed_single <- suppressMessages(
  processData(d_single_dose, verbose = FALSE)
)

test_that("get_pooled_data returns NA sentinels when no repeated dosing exists", {
  result <- get_pooled_data(
    processed_single$dat,
    dose_type   = "repeated_doses",
    pooled_ctrl = pooled_control()
  )

  expect_type(result, "list")
  expect_identical(result$datpooled_efd, NA)
})

test_that("get_pooled_data combined_doses path skips repeated section when single-dose only", {
  result <- get_pooled_data(
    processed_single$dat,
    dose_type   = "combined_doses",
    pooled_ctrl = pooled_control()
  )

  expect_type(result, "list")
  expect_identical(result$datpooled_efd, NA)
  expect_identical(result$datpooled_all, NA)
  expect_false(identical(result$datpooled_fd, NA))
})
