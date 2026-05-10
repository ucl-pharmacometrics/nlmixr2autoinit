result_oral <- processData(Oral_1CPT, verbose = FALSE)
ctrl        <- initsControl()

dose_type_oral <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics == "Dose Type"]
route_oral     <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics == "Dose Route"]

pooled_oral <- get_pooled_data(result_oral$dat, dose_type_oral, ctrl$pooled.control)

nca_oral <- run_pooled_nca(
  dat         = result_oral$dat,
  dose_type   = dose_type_oral,
  pooled      = pooled_oral,
  route       = route_oral,
  pooled_ctrl = pooled_control(),
  nca_ctrl    = nca_control()
)

ka_wn <- ka_wanger_nelson(
  pooled_oral$datpooled_fd$binned.df,
  nca_oral$nca.fd.results
)

test_that("ka_wanger_nelson returns a list", {
  expect_type(ka_wn, "list")
})

test_that("ka_wanger_nelson list contains ka", {
  expect_true("ka" %in% names(ka_wn))
})

test_that("ka_wanger_nelson ka is numeric", {
  expect_true(is.numeric(ka_wn$ka))
})

test_that("ka_wanger_nelson ka is positive when NCA results are valid", {
  skip_if(is.na(nca_oral$nca.fd.results$clobs) || is.na(nca_oral$nca.fd.results$vzobs),
          "NCA results not available")
  skip_if(is.na(ka_wn$ka), "ka could not be estimated")
  expect_gt(ka_wn$ka, 0)
})
