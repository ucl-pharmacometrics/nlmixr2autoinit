# Create full control list for initial parameter estimation

Aggregates modular control functions into a structured list for use in
population pharmacokinetic parameter initialization.

## Usage

``` r
initsControl(
  ss.control = ss_control(),
  pooled.control = pooled_control(),
  nca.control = nca_control(),
  fallback.control = fallback_control(),
  selmetrics = "rRMSE2",
  hybrid.base = TRUE,
  preferNCA = TRUE
)
```

## Arguments

- ss.control:

  A control list consistent with the structure returned by ss_control().

- pooled.control:

  A control list consistent with the structure returned by
  pooled_control().

- nca.control:

  A control list consistent with the structure returned by
  nca_control().

- fallback.control:

  A control list consistent with the structure returned by
  fallback_control().

- selmetrics:

  A character string or vector specifying model performance metrics to
  evaluate. Must be one or more of "APE", "MAE", "MAPE", "RMSE",
  "rRMSE1", or "rRMSE2". Default is "rRMSE2".

- hybrid.base:

  Logical. If TRUE, enables hybrid evaluation mode in which model
  performance is assessed using mixed parameter combinations across
  methods. If FALSE, each method is evaluated independently. Default is
  TRUE.

- preferNCA:

  Logical. If TRUE and selmetrics equals "rRMSE2", the lowest rRMSE2 is
  selected first. If the best-performing method is not NCA-based, the
  function then checks whether an NCA-based method offers a lower
  rRMSE1. If so, the NCA method is selected. Default is TRUE.

## Value

A named list combining all control modules for parameter estimation.

## See also

[ss_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md),
[pooled_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md),
[nca_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/nca_control.md),
[fallback_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/fallback_control.md)

## Examples

``` r
initsControl(
  pooled.control = pooled_control(nbins = 8),
  fallback.control = fallback_control(
    sigma_method_additive = "fixed_fraction"
  )
)
#> $ss.control
#> $ss.control$ss_method
#> [1] "combined"
#> 
#> $ss.control$no.doses
#> [1] 5
#> 
#> $ss.control$no.half_lives
#> [1] 5
#> 
#> $ss.control$allowed_interval_variation
#> [1] 0.25
#> 
#> $ss.control$allowed_dose_variation
#> [1] 0.2
#> 
#> $ss.control$min_doses_required
#> [1] 3
#> 
#> $ss.control$tad_rounding
#> [1] TRUE
#> 
#> 
#> $pooled.control
#> $pooled.control$nbins
#> [1] 8
#> 
#> $pooled.control$bin_method
#> [1] "quantile"
#> 
#> $pooled.control$tad_rounding
#> [1] TRUE
#> 
#> 
#> $nca.control
#> $nca.control$trapezoidal.rule
#> [1] "linear_up_log_down" "linear"            
#> 
#> $nca.control$duration
#> NULL
#> 
#> $nca.control$nlastpoints
#> [1] 3
#> 
#> $nca.control$slope.method
#> [1] "bestfitforce"
#> 
#> 
#> $fallback.control
#> $fallback.control$enable_ka_fallback
#> [1] TRUE
#> 
#> $fallback.control$sigma_method_additive
#> [1] "fixed_fraction"
#> 
#> $fallback.control$sigma_method_proportional
#> [1] "model"
#> 
#> $fallback.control$sigma_fallback_fraction
#> [1] 0.2
#> 
#> 
#> $selmetrics
#> [1] "rRMSE2"
#> 
#> $hybrid.base
#> [1] TRUE
#> 
#> $preferNCA
#> [1] TRUE
#> 
```
