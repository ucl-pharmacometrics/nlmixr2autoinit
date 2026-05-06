# Parameter sweeping for a one-compartment Michaelis-Menten model

Performs parameter sweeping by varying pharmacokinetic parameters in a
one-compartment model with Michaelis-Menten elimination.

## Usage

``` r
sim_sens_1cmpt_mm(
  dat,
  sim_vmax = list(mode = "auto", values = NULL, est.cl = NULL),
  sim_km = list(mode = "auto", values = NULL),
  sim_vd = list(mode = "manual", values = NULL),
  sim_ka = list(mode = "manual", values = NULL),
  route = c("iv", "oral"),
  verbose = TRUE
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- sim_vmax:

  List specifying Vmax:

  - mode: "manual" or "auto"

  - values: numeric vector if mode = "manual"

  - est.cl: required if mode = "auto"

- sim_km:

  List specifying Km:

  - mode: "manual" or "auto"

  - values: numeric vector if mode = "manual"

- sim_vd:

  List specifying Vd:

  - mode: must be "manual"

  - values: numeric vector

- sim_ka:

  List specifying Ka (oral route only):

  - mode: must be "manual"

  - values: numeric vector

- route:

  Dosing route, either "iv" or "oral". Default is "iv".

- verbose:

  Logical (default = TRUE). Controls whether progress information is
  displayed during parameter sweeping. When TRUE, a dynamic progress bar
  is shown using the `progressr` package to indicate simulation status
  and elapsed time. When FALSE, progress output is suppressed and the
  function runs silently.

## Value

A data frame containing parameter combinations and model fit metrics.

## Details

The function generates a parameter grid and performs model fitting for
each combination using `Fit_1cmpt_mm_iv`. Parameters can be specified
manually or automatically derived. Model predictions and fit metrics are
computed for each simulation to assess parameter sensitivity.

## See also

[`Fit_1cmpt_mm_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_iv.md),
[`Fit_1cmpt_mm_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
# Example: IV dosing scenario with automatic Vmax and Km
out <- sim_sens_1cmpt_mm(
  dat = Bolus_1CPTMM[Bolus_1CPTMM$ID<50,],
  sim_vmax = list(mode = "auto", est.cl = 4),
  sim_km   = list(mode = "auto"),
  sim_vd   = list(mode = "manual", values = 70),
  sim_ka   = list(mode = "manual", values = NA),
  route = "iv"
)
head(out[out$rRMSE2==min(out$rRMSE2),])
#> # A tibble: 1 × 11
#>    Vmax    Km    Vd    Ka      APE   MAE  MAPE  RMSE rRMSE1 rRMSE2
#>   <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 1427.  178.    70    NA 1697957.  599.  36.5 1286.   78.9   63.1
#> # ℹ 1 more variable: Cumulative.Time.Sec <dbl>
# }
```
