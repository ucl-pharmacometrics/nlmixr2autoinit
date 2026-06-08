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

test_that("processData converts character IDs to integers", {
  d <- data.frame(
    ID   = c("A", "A", "B", "B"),
    EVID = c(1, 0, 1, 0),
    CMT  = c(1, 2, 1, 2),
    AMT  = c(1, 0, 1, 0),
    TIME = c(0, 1, 0, 1),
    DV   = c(NA, 1.5, NA, 1)
  )
  result <- processData(d, verbose = FALSE)
  expect_type(result$dat$ID, "integer")
  expect_equal(sort(unique(result$dat$ID)), 1:2)
})

test_that("processData converts factor IDs to integers", {
  d <- data.frame(
    ID   = factor(c("A", "A", "B", "B")),
    EVID = c(1, 0, 1, 0),
    CMT  = c(1, 2, 1, 2),
    AMT  = c(1, 0, 1, 0),
    TIME = c(0, 1, 0, 1),
    DV   = c(NA, 1.5, NA, 1)
  )
  result <- processData(d, verbose = FALSE)
  expect_type(result$dat$ID, "integer")
  expect_equal(sort(unique(result$dat$ID)), 1:2)
})

test_that("processData output is identical for numeric, character, and factor IDs", {
  d_num <- Oral_1CPT[Oral_1CPT$ID %in% 1:5, ]
  d_chr <- d_num
  d_chr$ID <- c("a", "b", "c", "d", "e")[d_chr$ID]
  d_fac <- d_num
  d_fac$ID <- as.factor(d_chr$ID)

  r_num <- processData(d_num, verbose = FALSE)
  r_chr <- processData(d_chr, verbose = FALSE)
  r_fac <- processData(d_fac, verbose = FALSE)

  expect_equal(r_chr$dat,      r_num$dat)
  expect_equal(r_fac$dat,      r_num$dat)
  expect_equal(r_chr$Datainfo, r_num$Datainfo)
  expect_equal(r_fac$Datainfo, r_num$Datainfo)
})
