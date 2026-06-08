# nlmixr2autoinit (development version)

## Bug fixes

- `getPPKinits()` now supports non-integer subject IDs (character or factor).
  IDs are converted to integers via `as.integer(as.factor(ID))` during data
  processing. (#6)
- `sim_sens_1cmpt_mm()`, `sim_sens_2cmpt()`, and `sim_sens_3cmpt()` no longer
  abort the entire parameter sweep when a single `Fit_*()` call errors
  (e.g. `lotri` rejecting NA/Inf/0 initial values with "subscript out of
  bounds"). Failing rows now produce `NA` metric values, which downstream
  selection logic already handles via `min(..., na.rm = TRUE)`.
