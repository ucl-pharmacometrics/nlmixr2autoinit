# nlmixr2autoinit (development version)

## Bug fixes

- `getPPKinits()` now supports non-integer subject IDs (character or factor).
  IDs are converted to integers via `as.integer(as.factor(ID))` during data
  processing. (#6)
