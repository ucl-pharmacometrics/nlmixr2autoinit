test_that("bolus route runs end-to-end", {
  skip_on_cran()
  result <- getPPKinits(Bolus_1CPT)
  expect_s3_class(result, "getPPKinits")
})

test_that("oral route runs end-to-end", {
  skip_on_cran()
  result <- getPPKinits(Oral_1CPT)
  expect_s3_class(result, "getPPKinits")
})
