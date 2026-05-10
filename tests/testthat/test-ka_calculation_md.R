
result_ok <- ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 852, Dose = 60000, tau = 24)
result_high_ct <- ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 1e9, Dose = 60000, tau = 24)

test_that("ka_calculation_md returns a list", {
  expect_type(result_ok, "list")
})

test_that("ka_calculation_md list contains ka, full_solution, message", {
  expect_true("ka"            %in% names(result_ok))
  expect_true("full_solution" %in% names(result_ok))
  expect_true("message"       %in% names(result_ok))
})

test_that("ka_calculation_md ka is a positive numeric for valid input", {
  expect_true(is.numeric(result_ok$ka))
  expect_gt(result_ok$ka, 0)
})

test_that("ka_calculation_md ka is greater than ke (no flip-flop)", {
  ke <- 0.057
  expect_gt(result_ok$ka, ke)
})

test_that("ka_calculation_md message is 'complete' for valid input", {
  expect_equal(result_ok$message, "complete")
})

test_that("ka_calculation_md returns NA ka when Ct exceeds theoretical max", {
  expect_true(is.na(result_high_ct$ka))
})

test_that("ka_calculation_md message contains informative text when Ct too high", {
  expect_match(result_high_ct$message, "exceeds", ignore.case = TRUE)
})

test_that("ka_calculation_md returns NA when boundary signs are equal", {
  result_bad <- ka_calculation_md(cl = 0.001, ke = 999, t = 1, Ct = 1, Dose = 100, tau = 24)
  expect_true(is.na(result_bad$ka))
})

test_that("ka_calculation_md full_solution is not NULL for valid input", {
  expect_false(is.null(result_ok$full_solution))
})
