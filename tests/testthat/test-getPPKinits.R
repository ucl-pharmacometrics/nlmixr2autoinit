test_that("getPPKinits runs on example bolus data", {
  skip_on_cran()

  data("Bolus_1CPT", package = "nlmixr2data")

  out <- getPPKinits(dat = Bolus_1CPT)

  expect_true(!is.null(out))
})

test_that("getPPKinits runs on example oral data", {
  skip_on_cran()

  data("Oral_1CPT", package = "nlmixr2data")

  out <- getPPKinits(dat = Oral_1CPT)

  expect_true(!is.null(out))
})
