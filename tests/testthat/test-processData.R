result_bolus <- processData(Bolus_1CPT, verbose = FALSE)
result_oral  <- processData(Oral_1CPT,  verbose = FALSE)

test_that("processData returns a list", {
  expect_type(result_bolus, "list")
})

test_that("processData list contains dat and Datainfo", {
  expect_true("dat"      %in% names(result_bolus))
  expect_true("Datainfo" %in% names(result_bolus))
})

test_that("processData dat is a data frame", {
  expect_s3_class(result_bolus$dat, "data.frame")
})

test_that("processData Datainfo contains Dose Type and Dose Route rows", {
  info <- result_bolus$Datainfo
  expect_true("Dose Type"  %in% info$Infometrics)
  expect_true("Dose Route" %in% info$Infometrics)
})

test_that("processData identifies bolus route for Bolus_1CPT", {
  route_val <- result_bolus$Datainfo$Value[result_bolus$Datainfo$Infometrics == "Dose Route"]
  expect_equal(route_val, "bolus")
})

test_that("processData identifies oral route for Oral_1CPT", {
  route_val <- result_oral$Datainfo$Value[result_oral$Datainfo$Infometrics == "Dose Route"]
  expect_equal(route_val, "oral")
})

test_that("processData dat retains required columns", {
  expect_true(all(c("ID", "TIME", "EVID", "DV") %in% colnames(result_bolus$dat)))
})

# ---- Filter subjects with no observations (GH-2) ----

test_that("processData removes subjects with no observations", {
  d_mixed <- data.frame(
    ID   = c(1, 2, 2),
    EVID = c(1, 1, 0),
    CMT  = c(1, 1, 1),
    AMT  = c(100, 100, 0),
    TIME = c(0, 0, 1),
    DV   = c(NA, NA, 0.5)
  )

  result_mixed <- suppressMessages(processData(d_mixed, verbose = FALSE))

  expect_equal(unique(result_mixed$dat$ID), 2)
  n_subj <- result_mixed$Datainfo$Value[
    result_mixed$Datainfo$Infometrics == "Number of Subjects"
  ]
  expect_equal(n_subj, "1")
})

test_that("processData keeps all subjects when every subject has observations", {
  n_ids_before <- dplyr::n_distinct(Bolus_1CPT$ID)
  n_ids_after  <- dplyr::n_distinct(result_bolus$dat$ID)
  expect_equal(n_ids_before, n_ids_after)
})

test_that("processData errors when no subject has observations", {
  d_no_obs <- data.frame(
    ID   = c(1, 1, 2),
    EVID = c(1, 1, 1),
    CMT  = c(1, 1, 1),
    AMT  = c(100, 100, 100),
    TIME = c(0, 24, 0),
    DV   = c(NA, NA, NA)
  )

  expect_error(
    suppressMessages(processData(d_no_obs, verbose = FALSE)),
    "No subjects with observations"
  )
})
