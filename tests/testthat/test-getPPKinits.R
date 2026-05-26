# Use a small subset so the top-level setup is fast.
# IDs 1-5 map to "a"-"e" in alphabetical order, preserving subject order.
d_subset     <- Oral_1CPT[Oral_1CPT$ID %in% 1:5, ]
d_subset_chr <- d_subset
d_subset_chr$ID <- c("a", "b", "c", "d", "e")[d_subset_chr$ID]
d_subset_fac <- d_subset
d_subset_fac$ID <- as.factor(d_subset_chr$ID)

result_num <- getPPKinits(d_subset,     verbose = FALSE)
result_chr <- getPPKinits(d_subset_chr, verbose = FALSE)
result_fac <- getPPKinits(d_subset_fac, verbose = FALSE)

test_that("getPPKinits returns getPPKinits environment for numeric IDs", {
  expect_s3_class(result_num, "getPPKinits")
})

test_that("getPPKinits accepts character subject IDs (issue 6)", {
  expect_s3_class(result_chr, "getPPKinits")
})

test_that("getPPKinits accepts factor subject IDs (issue 6)", {
  expect_s3_class(result_fac, "getPPKinits")
})

test_that("getPPKinits with character IDs produces same estimates as numeric IDs", {
  expect_equal(
    result_num$Recommended_initial_estimates,
    result_chr$Recommended_initial_estimates
  )
})

test_that("getPPKinits with factor IDs produces same estimates as numeric IDs", {
  expect_equal(
    result_num$Recommended_initial_estimates,
    result_fac$Recommended_initial_estimates
  )
})
