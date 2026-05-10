result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)
ctrl         <- initsControl()

dose_type_bolus <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Type"]
dose_type_oral  <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Type"]
route_bolus     <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Route"]
route_oral      <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Route"]

pooled_bolus <- get_pooled_data(result_bolus$dat, dose_type_bolus, ctrl$pooled.control)
pooled_oral  <- get_pooled_data(result_oral$dat,  dose_type_oral,  ctrl$pooled.control)

nca_bolus <- run_pooled_nca(
  dat         = result_bolus$dat,
  dose_type   = dose_type_bolus,
  pooled      = pooled_bolus,
  route       = route_bolus,
  pooled_ctrl = pooled_control(),
  nca_ctrl    = nca_control()
)

nca_oral <- run_pooled_nca(
  dat         = result_oral$dat,
  dose_type   = dose_type_oral,
  pooled      = pooled_oral,
  route       = route_oral,
  pooled_ctrl = pooled_control(),
  nca_ctrl    = nca_control()
)

test_that("run_pooled_nca returns a list for bolus", {
  expect_type(nca_bolus, "list")
})

test_that("run_pooled_nca contains nca.fd.results", {
  expect_true("nca.fd.results" %in% names(nca_bolus))
})

test_that("run_pooled_nca nca.fd.results contains clobs and vzobs", {
  expect_true("clobs" %in% names(nca_bolus$nca.fd.results))
  expect_true("vzobs" %in% names(nca_bolus$nca.fd.results))
})

test_that("run_pooled_nca nca.fd.results clobs is numeric", {
  expect_true(is.numeric(nca_bolus$nca.fd.results$clobs))
})

test_that("run_pooled_nca contains nca.all.results", {
  expect_true("nca.all.results" %in% names(nca_bolus))
})

test_that("run_pooled_nca returns a list for oral", {
  expect_type(nca_oral, "list")
})
