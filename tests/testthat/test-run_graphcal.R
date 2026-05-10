result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)
ctrl         <- initsControl()

dose_type_bolus <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Type"]
dose_type_oral  <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Type"]
route_bolus     <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Route"]
route_oral      <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics  == "Dose Route"]

pooled_bolus <- get_pooled_data(result_bolus$dat, dose_type_bolus, ctrl$pooled.control)
pooled_oral  <- get_pooled_data(result_oral$dat,  dose_type_oral,  ctrl$pooled.control)

graph_bolus <- run_graphcal(
  dat         = result_bolus$dat,
  route       = route_bolus,
  data_type   = dose_type_bolus,
  pooled      = pooled_bolus,
  pooled_ctrl = ctrl$pooled.control,
  nlastpoints = 3
)

graph_oral <- run_graphcal(
  dat         = result_oral$dat,
  route       = route_oral,
  data_type   = dose_type_oral,
  pooled      = pooled_oral,
  pooled_ctrl = ctrl$pooled.control,
  nlastpoints = 3
)

test_that("run_graphcal returns a list for bolus", {
  expect_type(graph_bolus, "list")
})

test_that("run_graphcal contains cl and vd for bolus", {
  expect_true("cl" %in% names(graph_bolus))
  expect_true("vd" %in% names(graph_bolus))
})

test_that("run_graphcal cl and vd are positive numerics for bolus", {
  expect_gt(graph_bolus$cl, 0)
  expect_gt(graph_bolus$vd, 0)
})

test_that("run_graphcal contains ka for oral", {
  expect_true("ka" %in% names(graph_oral))
})

test_that("run_graphcal contains time.spent", {
  expect_true("time.spent" %in% names(graph_bolus))
})
