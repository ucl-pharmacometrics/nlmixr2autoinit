# Evaluates predictive performance of a one-compartment model

Computes predictive error metrics by comparing simulated and observed
concentration–time data using specified pharmacokinetic parameters and
dosing route.

## Usage

``` r
eval_perf_1cmpt(
  dat,
  est.method = "rxSolve",
  ka = NULL,
  cl = NULL,
  vd = NULL,
  route = c("bolus", "infusion", "oral")
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- est.method:

  Estimation method passed to the fitting function. Defaults to using
  `rxSolve` for model simulation and parameter estimation.

- ka:

  Absorption rate constant.

- cl:

  Clearance value.

- vd:

  Volume of distribution.

- route:

  A character string indicating the route of administration. Must be one
  of `"oral"`, `"infusion"`, or `"bolus"`. Defaults to `"bolus"`.

## Value

A numeric vector containing absolute prediction error, mean absolute
error, mean absolute percentage error, root mean square error, and
relative root mean square error.

## Details

Internally selects the appropriate one-compartment model fitting
function, using
[`Fit_1cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_oral.md)
for oral administration and
[`Fit_1cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_iv.md)
for intravenous administration. Predictive performance is quantified
using the
[`metrics.()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/metrics..md)
function.

## See also

[`Fit_1cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_oral.md),
[`Fit_1cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_iv.md),
[`metrics.`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/metrics..md)

## Examples

``` r
eval_perf_1cmpt(
  dat = Oral_1CPT,
  est.method = "rxSolve",
  ka = 1,
  cl = 4,
  vd = 70,
  route = "oral"
)
#>        APE        MAE       MAPE       RMSE     rRMSE1     rRMSE2 
#> 921381.710    132.630     75.620    230.715     45.693     48.466 
```
