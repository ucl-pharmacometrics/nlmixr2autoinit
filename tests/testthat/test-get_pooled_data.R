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
