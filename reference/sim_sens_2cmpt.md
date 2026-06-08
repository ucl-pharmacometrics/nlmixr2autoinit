# Parameter sweeping for a two-compartment pharmacokinetic model

Performs parameter sweeping by varying pharmacokinetic parameters in a
two-compartment model under IV or oral dosing. Model fit is evaluated
across combinations of CL, Vc, Vp, Q, and Ka (oral only).

## Usage

``` r
sim_sens_2cmpt(
  dat,
  sim_ka = list(mode = "manual", values = NULL),
  sim_cl = list(mode = "manual", values = NULL),
  sim_vc = list(mode = "manual", values = NULL),
  sim_vp = list(mode = c("auto", "manual"), values = NULL),
  sim_q = list(mode = c("auto", "manual"), values = NULL, auto.strategy = c("scaled",
    "fixed")),
  route = c("iv", "oral"),
  ncores = 2,
  verbose = TRUE
)
```

## Arguments

- dat:

  Pharmacokinetic dataset.

- sim_ka:

  List specifying Ka (oral route only):

  - mode: must be "manual"

  - values: numeric vector

- sim_cl:

  List specifying clearance (CL):

  - mode: must be "manual"

  - values: numeric vector

- sim_vc:

  List specifying central volume (Vc):

  - mode: must be "manual"

  - values: numeric vector

- sim_vp:

  List specifying peripheral volume (Vp):

  - mode: "manual" or "auto"

  - values: numeric vector if manual

- sim_q:

  List specifying inter-compartmental clearance (Q):

  - mode: "manual" or "auto"

  - values: numeric vector if manual

- route:

  Dosing route, either "iv" or "oral". Default is "iv".

- ncores:

  Number of cores to use for parallelization, passed to `rxControl()`.
  Default is 2.

- verbose:

  Logical (default = TRUE). Controls whether progress information is
  displayed during parameter sweeping. When TRUE, a dynamic progress bar
  is shown using the `progressr` package to indicate simulation status
  and elapsed time. When FALSE, progress output is suppressed and the
  function runs silently.

## Value

A data frame containing parameter combinations with model fit metrics.

## Details

The function generates a parameter grid and performs model fitting for
each combination using `Fit_2cmpt_iv` or `Fit_2cmpt_oral`. Parameters
can be specified manually or automatically derived. Model predictions
and fit metrics are computed for each simulation to assess parameter
sensitivity.

## See also

[`Fit_2cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_iv.md),
[`Fit_2cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
out <- sim_sens_2cmpt(
  dat = Bolus_2CPT[Bolus_2CPT$ID<6,],
  sim_cl = list(mode = "manual", values = 4),
  sim_vc = list(mode = "manual", values = 50),
  sim_vp = list(mode = "auto"),
  sim_q  = list(mode = "auto"),
  sim_ka = list(mode = "manual", values = NA),
  route = "iv",verbose=FALSE
)
head(out[out$rRMSE2==min(out$rRMSE2),])
#> # A tibble: 1 × 12
#>      Vc    Vp     Q    CL Ka       APE   MAE  MAPE  RMSE rRMSE1 rRMSE2
#>   <dbl> <dbl> <dbl> <dbl> <lgl>  <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    50   100     8     4 NA    53184.  183.  63.4  309.   55.5   54.5
#> # ℹ 1 more variable: Cumulative.Time.Sec <dbl>
# }
```
