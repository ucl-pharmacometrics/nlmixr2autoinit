# Parameter sweeping for a three-compartment pharmacokinetic model

Performs parameter sweeping by varying pharmacokinetic parameters in a
three-compartment model under IV or oral dosing. Parameter combinations
include Vc, Vp1, Vp2, Q1, Q2, CL, and Ka (oral only).

## Usage

``` r
sim_sens_3cmpt(
  dat,
  sim_vc = list(mode = "manual", values = NULL),
  sim_vp = list(mode = c("auto", "manual"), values = NULL),
  sim_vp2 = list(mode = c("auto", "manual"), values = NULL),
  sim_q = list(mode = c("auto", "manual"), values = NULL, auto.strategy = c("scaled",
    "fixed")),
  sim_q2 = list(mode = c("auto", "manual"), values = NULL, auto.strategy = c("scaled",
    "fixed")),
  sim_cl = list(mode = "manual", values = NULL),
  sim_ka = list(mode = "manual", values = NULL),
  route = c("iv", "oral"),
  verbose = TRUE
)
```

## Arguments

- dat:

  Pharmacokinetic dataset.

- sim_vc:

  List specifying Vc:

  - mode: must be "manual"

  - values: numeric vector

- sim_vp:

  List specifying Vp1:

  - mode: "manual" or "auto"

  - values: numeric vector if manual

- sim_vp2:

  List specifying Vp2:

  - mode: "manual" or "auto"

  - values: numeric vector if manual

- sim_q:

  List specifying Q1:

  - mode: "manual" or "auto"

  - values: numeric vector if manual

  - auto.strategy: "scaled" or "fixed" when auto

- sim_q2:

  List specifying Q2:

  - mode: "manual" or "auto"

  - values: numeric vector if manual

  - auto.strategy: "scaled" or "fixed" when auto

- sim_cl:

  List specifying CL:

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

A data frame containing parameter combinations with model fit metrics.

## Details

The function generates a parameter grid and evaluates each combination
using `Fit_3cmpt_iv` or `Fit_3cmpt_oral`. Model predictions and fit
metrics are calculated for each simulation to assess parameter influence
and identify optimal regions of the parameter space. Parameters can be
provided manually or derived automatically.

## See also

[`Fit_3cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_iv.md),
[`Fit_3cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
out <- sim_sens_3cmpt(
  dat = Bolus_2CPT,
  sim_cl = list(mode = "manual", values = 4),
  sim_vc = list(mode = "manual", values = 50),
  sim_vp = list(mode = "auto"),
  sim_vp2 = list(mode = "auto"),
  sim_q  = list(mode = "auto", auto.strategy = "scaled"),
  sim_q2 = list(mode = "auto", auto.strategy = "scaled"),
  route = "iv",verbose=FALSE
)
head(out[out$rRMSE2==min(out$rRMSE2),])
#> # A tibble: 1 × 14
#>      Vc   Vp1   Vp2    Q1    Q2    CL    Ka      APE   MAE  MAPE  RMSE rRMSE1
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1    50    10    50     8     8     4    NA 1118670.  161.  38.5  281.   46.9
#> # ℹ 2 more variables: rRMSE2 <dbl>, Cumulative.Time.Sec <dbl>
# }
```
