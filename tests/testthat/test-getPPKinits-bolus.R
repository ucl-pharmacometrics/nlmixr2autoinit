test_that("getPPKinits runs on example bolus data", {
  skip_on_cran()

  out <- getPPKinits(dat = Bolus_1CPT)

  expect_true(!is.null(out))
})
