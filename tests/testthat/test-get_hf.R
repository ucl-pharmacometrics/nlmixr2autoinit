result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)
ctrl         <- initsControl()

dose_type_bolus <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Type"]
dose_type_oral  <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Type"]

pooled_bolus <- get_pooled_data(result_bolus$dat, dose_type_bolus, ctrl$pooled.control)
pooled_oral  <- get_pooled_data(result_oral$dat,  dose_type_oral,  ctrl$pooled.control)

hf_bolus <- get_hf(dat = result_bolus$dat, pooled = pooled_bolus, verbose = FALSE)
hf_oral  <- get_hf(dat = result_oral$dat,  pooled = pooled_oral,  verbose = FALSE)

test_that("get_hf returns a list", {
  expect_type(hf_bolus, "list")
})

test_that("get_hf list contains half_life_median", {
  expect_true("half_life_median" %in% names(hf_bolus))
})

test_that("get_hf half_life_median is a positive finite numeric for bolus", {
  expect_true(is.numeric(hf_bolus$half_life_median))
  expect_gt(hf_bolus$half_life_median, 0)
  expect_true(is.finite(hf_bolus$half_life_median))
})

test_that("get_hf half_life_median is a positive finite numeric for oral", {
  expect_true(is.numeric(hf_oral$half_life_median))
  expect_gt(hf_oral$half_life_median, 0)
  expect_true(is.finite(hf_oral$half_life_median))
})
