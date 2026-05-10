dat_bolus <- processData(Bolus_1CPT, verbose = FALSE)$dat
dat_oral  <- processData(Oral_1CPT,  verbose = FALSE)$dat

sigma_bolus <- getsigma(dat_bolus)
sigma_oral  <- getsigma(dat_oral)

test_that("getsigma returns a list for bolus data", {
  expect_type(sigma_bolus, "list")
})

test_that("getsigma list contains summary element", {
  expect_true("summary" %in% names(sigma_bolus))
})

test_that("getsigma summary contains sigma_additive and sigma_proportional", {
  expect_true("sigma_additive"     %in% names(sigma_bolus$summary))
  expect_true("sigma_proportional" %in% names(sigma_bolus$summary))
})

test_that("getsigma sigma values are numeric", {
  expect_true(is.numeric(sigma_bolus$summary$sigma_additive))
  expect_true(is.numeric(sigma_bolus$summary$sigma_proportional))
})

test_that("getsigma works for oral data", {
  expect_type(sigma_oral, "list")
  expect_true("summary" %in% names(sigma_oral))
})
