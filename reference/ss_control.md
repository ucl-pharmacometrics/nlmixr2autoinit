# Internal control builder for steady-state evaluation

Constructs a list of control parameters used by is_ss() to determine
pharmacokinetic steady state.

## Usage

``` r
ss_control(
  ss_method = c("combined", "fixed_doses", "half_life_based"),
  no.doses = 5,
  no.half_lives = 5,
  allowed_interval_variation = 0.25,
  allowed_dose_variation = 0.2,
  min_doses_required = 3,
  tad_rounding = TRUE
)
```

## Arguments

- ss_method:

  Character string specifying the method used to determine steady state.
  One of:

  - "combined" (default): uses the smaller of the dose-based estimate
    (no.doses) and the half-life-based estimate (no.half_lives)

  - "fixed_doses": considers steady state reached after no.doses
    administrations

  - "half_life_based": uses the drug half-life and dosing interval to
    estimate the required number of doses

- no.doses:

  Integer indicating the number of doses assumed necessary to reach
  steady state when using the "fixed_doses" method or as part of the
  "combined" method. Default is 5.

- no.half_lives:

  Integer indicating the number of half-lives required to reach steady
  state when using the "half_life_based" or "combined" method. Default
  is 5.

- allowed_interval_variation:

  Numeric value specifying the acceptable fractional variation in dose
  interval. For example, 0.25 allows plus or minus 25 percent variation.
  Default is 0.25.

- allowed_dose_variation:

  Numeric value specifying the acceptable fractional variation in dose
  amount. For example, 0.20 allows plus or minus 20 percent variation.
  Default is 0.20.

- min_doses_required:

  Integer specifying the minimum number of doses that must be
  administered regardless of method. Default is 3.

- tad_rounding:

  Logical value. If TRUE (default), rounding is applied when comparing
  time after dose (tad) to dosing intervals to allow small numerical
  deviations.

## Value

A named list containing the steady-state control parameters, typically
passed as the ssctrl argument to
[`is_ss()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/is_ss.md).

## See also

[is_ss](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/is_ss.md)

## Examples

``` r
ss_control()
#> $ss_method
#> [1] "combined"
#> 
#> $no.doses
#> [1] 5
#> 
#> $no.half_lives
#> [1] 5
#> 
#> $allowed_interval_variation
#> [1] 0.25
#> 
#> $allowed_dose_variation
#> [1] 0.2
#> 
#> $min_doses_required
#> [1] 3
#> 
#> $tad_rounding
#> [1] TRUE
#> 
ss_control(ss_method = "fixed_doses", no.doses = 4)
#> $ss_method
#> [1] "fixed_doses"
#> 
#> $no.doses
#> [1] 4
#> 
#> $no.half_lives
#> [1] 5
#> 
#> $allowed_interval_variation
#> [1] 0.25
#> 
#> $allowed_dose_variation
#> [1] 0.2
#> 
#> $min_doses_required
#> [1] 3
#> 
#> $tad_rounding
#> [1] TRUE
#> 
```
