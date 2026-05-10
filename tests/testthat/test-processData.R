result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)

test_that("processData returns a list", {
  expect_type(result_bolus, "list")
})

test_that("processData list contains dat and Datainfo", {
  expect_true("dat"      %in% names(result_bolus))
  expect_true("Datainfo" %in% names(result_bolus))
})

test_that("processData dat is a data frame", {
  expect_s3_class(result_bolus$dat, "data.frame")
})

test_that("processData Datainfo contains Dose Type and Dose Route rows", {
  info <- result_bolus$Datainfo
  expect_true("Dose Type"  %in% info$Infometrics)
  expect_true("Dose Route" %in% info$Infometrics)
})

test_that("processData identifies bolus route for Bolus_1CPT", {
  route_val <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Route"]
  expect_equal(route_val, "bolus")
})

test_that("processData identifies oral route for Oral_1CPT", {
  route_val <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics == "Dose Route"]
  expect_equal(route_val, "oral")
})

test_that("processData dat retains required columns", {
  expect_true(all(c("ID", "TIME", "EVID", "DV") %in% colnames(result_bolus$dat)))
})
