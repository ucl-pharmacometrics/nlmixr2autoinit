# nlmixr2autoinit 1.0.1.9000

* Subjects with no observations are now removed inside `processData()` before
  route detection and downstream calculations. This fixes a crash in
  `getPPKinits()` when the input data contained subjects with dosing records
  but no observation records (GH-2).
* `get_pooled_data()` now returns gracefully (leaving `datpooled_efd` /
  `datpooled_all` as `NA`) when no subject has enough dosing rows to estimate
  a most-common interdose interval, instead of erroring inside the internal
  `tad_check` computation.
