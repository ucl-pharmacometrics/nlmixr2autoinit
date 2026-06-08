# nlmixr2autoinit 1.0.1.9000

* Subjects with no observations are now removed inside `processData()` before
  route detection and downstream calculations. This fixes a crash in
  `getPPKinits()` when the input data contained subjects with dosing records
  but no observation records (GH-2).
* `get_pooled_data()` now returns gracefully (leaving `datpooled_efd` /
  `datpooled_all` as `NA`) when no subject has enough dosing rows to estimate
  a most-common interdose interval, instead of erroring inside the internal
  `tad_check` computation.
* `getPPKinits()` now stops with a clear error when fewer than 2 observations
  remain after data filtering, instead of failing deep inside the
  one-compartment predictive-performance evaluation with a cryptic
  `colnames<-` "less than two dimensions" message (GH-2 follow-up).
* `hybrid_eval_perf_1cmpt()` errors clearly when called with no valid
  CL/Vd (or Ka for oral) parameter sources, rather than silently returning
  `NULL`.
