omega_out <- getOmegas()

test_that("getOmegas is callable without arguments", {
  expect_no_error(getOmegas())
})

test_that("getOmegas returns a data frame or list", {
  expect_true(is.data.frame(omega_out) || is.list(omega_out))
})

test_that("getOmegas output is not empty", {
  expect_gt(length(omega_out), 0)
})

test_that("getOmegas has more than zero rows when data frame", {
  if (is.data.frame(omega_out)) {
    expect_gt(nrow(omega_out), 0)
  } else {
    succeed()
  }
})
