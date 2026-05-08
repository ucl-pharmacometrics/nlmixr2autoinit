test_that("sim_sens_1cmpt_mm runs with IV example settings", {
  skip_on_cran()
  
  data("Bolus_1CPTMM", package = "nlmixr2data")
  
  out <- sim_sens_1cmpt_mm(
    dat = Bolus_1CPTMM[Bolus_1CPTMM$ID < 50, ],
    sim_vmax = list(mode = "auto", est.cl = 4),
    sim_km   = list(mode = "auto"),
    sim_vd   = list(mode = "manual", values = 70),
    sim_ka   = list(mode = "manual", values = NA),
    route = "iv",
    verbose = FALSE
  )
  
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) > 0)
})