
result_ok <- ka_calculation_sd(cl = 3.62, ke = 0.0556, t = 0.5, Ct = 310.6, Dose = 60000)
result_high_ct <- ka_calculation_sd(cl = 3.62, ke = 0.0556, t = 0.5, Ct = 1e9, Dose = 60000)

test_that("ka_calculation_sd returns a list", {
  expect_type(result_ok, "list")
})

test_that("ka_calculation_sd list contains ka, full_solution, message", {
  expect_true("ka"            %in% names(result_ok))
  expect_true("full_solution" %in% names(result_ok))
  expect_true("message"       %in% names(result_ok))
})

test_that("ka_calculation_sd ka is a positive numeric for valid input", {
  expect_true(is.numeric(result_ok$ka))
  expect_gt(result_ok$ka, 0)
})

test_that("ka_calculation_sd ka is greater than ke (no flip-flop)", {
  ke <- 0.0556
  expect_gt(result_ok$ka, ke)
})

test_that("ka_calculation_sd message is 'complete' for valid input", {
  expect_equal(result_ok$message, "complete")
})

test_that("ka_calculation_sd returns NA ka when Ct exceeds theoretical max", {
  expect_true(is.na(result_high_ct$ka))
})

test_that("ka_calculation_sd message contains informative text when Ct too high", {
  expect_match(result_high_ct$message, "exceeds", ignore.case = TRUE)
})

test_that("ka_calculation_sd returns NA when boundary signs are equal", {
  # ke so large that ka range is degenerate
  result_bad <- ka_calculation_sd(cl = 0.001, ke = 999, t = 0.5, Ct = 1, Dose = 100)
  expect_true(is.na(result_bad$ka))
})
