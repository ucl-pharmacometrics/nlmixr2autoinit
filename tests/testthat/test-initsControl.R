test_that("fallback_control", {
  expect_equal(
    fallback_control(),
    list(
      enable_ka_fallback = TRUE,
      sigma_method_additive = "model",
      sigma_method_proportional = "model",
      sigma_fallback_fraction = 0.2
    )
  )
})
