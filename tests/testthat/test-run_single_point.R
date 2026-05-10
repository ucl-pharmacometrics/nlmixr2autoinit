result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)
ctrl         <- initsControl()

dose_type_bolus <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Type"]
dose_type_oral  <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Type"]
route_bolus     <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Route"]
route_oral      <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Route"]

pooled_bolus <- get_pooled_data(result_bolus$dat, dose_type_bolus, ctrl$pooled.control)
pooled_oral  <- get_pooled_data(result_oral$dat,  dose_type_oral,  ctrl$pooled.control)

hf_bolus <- get_hf(dat = result_bolus$dat, pooled = pooled_bolus, verbose = FALSE)
hf_oral  <- get_hf(dat = result_oral$dat,  pooled = pooled_oral,  verbose = FALSE)

sp_bolus <- run_single_point(
  dat         = result_bolus$dat,
  route       = route_bolus,
  half_life   = hf_bolus$half_life_median,
  dose_type   = dose_type_bolus,
  pooled_ctrl = ctrl$pooled.control,
  ssctrl      = ctrl$ss.control
)

sp_oral <- run_single_point(
  dat         = result_oral$dat,
  route       = route_oral,
  half_life   = hf_oral$half_life_median,
  dose_type   = dose_type_oral,
  pooled_ctrl = ctrl$pooled.control,
  ssctrl      = ctrl$ss.control
)

test_that("run_single_point returns a list", {
  expect_type(sp_bolus, "list")
})

test_that("run_single_point contains singlepoint.results", {
  expect_true("singlepoint.results" %in% names(sp_bolus))
})

test_that("run_single_point singlepoint.results contains cl and vd", {
  expect_true("cl" %in% names(sp_bolus$singlepoint.results))
  expect_true("vd" %in% names(sp_bolus$singlepoint.results))
})

test_that("run_single_point cl and vd are positive numerics for bolus", {
  expect_gt(sp_bolus$singlepoint.results$cl, 0)
  expect_gt(sp_bolus$singlepoint.results$vd, 0)
})

test_that("run_single_point ka is NA for bolus", {
  expect_true(is.na(sp_bolus$singlepoint.results$ka))
})

test_that("run_single_point includes ka field for oral", {
  expect_true("ka" %in% names(sp_oral$singlepoint.results))
})

test_that("run_single_point contains approx.vc.out", {
  expect_true("approx.vc.out" %in% names(sp_bolus))
})
