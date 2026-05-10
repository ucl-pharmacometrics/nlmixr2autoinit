# Single-dose subset
dat_sd_raw <- Oral_1CPT[Oral_1CPT$SD == 1, ]
dat_sd     <- processData(dat_sd_raw, verbose = FALSE)$dat
dat_sd     <- is_ss(dat_sd)
result_sd  <- run_ka_solution(df = dat_sd, cl = 4, ke = 4 / 70, Fbio = 1)

# Mixed-dose (full dataset)
dat_all    <- processData(Oral_1CPT, verbose = FALSE)$dat
dat_all    <- is_ss(dat_all)
result_all <- run_ka_solution(df = dat_all, cl = 4, ke = 4 / 70, Fbio = 1)

test_that("run_ka_solution returns a list for single-dose data", {
  expect_type(result_sd, "list")
})

test_that("run_ka_solution list contains expected elements", {
  expect_true("ka_calc_median" %in% names(result_sd))
  expect_true("ka_calc_dat_sd" %in% names(result_sd))
  expect_true("ka_calc_dat_md" %in% names(result_sd))
  expect_true("ka_calc_dat"    %in% names(result_sd))
})

test_that("run_ka_solution ka_calc_median is a numeric scalar for single-dose", {
  expect_true(is.numeric(result_sd$ka_calc_median))
  expect_length(result_sd$ka_calc_median, 1)
})

test_that("run_ka_solution ka_calc_median is positive for single-dose", {
  skip_if(is.na(result_sd$ka_calc_median), "ka_calc_median is NA")
  expect_gt(result_sd$ka_calc_median, 0)
})

test_that("run_ka_solution ka_calc_dat_sd is a data frame for single-dose", {
  expect_s3_class(result_sd$ka_calc_dat_sd, "data.frame")
})

test_that("run_ka_solution ka_calc_dat_sd contains ka_calc column", {
  expect_true("ka_calc" %in% names(result_sd$ka_calc_dat_sd))
})

test_that("run_ka_solution returns a list for mixed-dose data", {
  expect_type(result_all, "list")
})

test_that("run_ka_solution ka_calc_median is positive for mixed-dose", {
  skip_if(is.na(result_all$ka_calc_median), "ka_calc_median is NA")
  expect_gt(result_all$ka_calc_median, 0)
})

test_that("run_ka_solution ka_calc_dat contains ka_calcv column", {
  expect_true("ka_calcv" %in% names(result_all$ka_calc_dat))
})
