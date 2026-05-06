# Control options for non-compartmental analysis

Control options for non-compartmental analysis (NCA)

## Usage

``` r
nca_control(
  trapezoidal.rule = c("linear_up_log_down", "linear"),
  duration = NULL,
  nlastpoints = 3,
  slope.method = "bestfitforce"
)
```

## Arguments

- trapezoidal.rule:

  Character. Method for trapezoidal AUC integration:

  - `"linear"` - Linear trapezoidal rule (default)

  - `"linear_up_log_down"` - Linear-up / log-down rule

- duration:

  Numeric. Optional. Duration of the observation window (same units as
  time). Used to restrict the integration or define the evaluation
  range.

- nlastpoints:

  Integer. Number of terminal points for half-life regression (default =
  3).

- slope.method:

  Character. Method for estimating the terminal slope (lambdaz). Options
  are:

  - "bestfit": Performs automated terminal phase selection based on
    adjusted R-squared. If no acceptable segment is found, lambdaz is
    returned as NA.

  - "bestfitforce": First attempts the "bestfit" method. If no valid
    lambdaz is obtained, the function applies a fallback log-linear
    regression using progressively fewer terminal points to force
    estimation. This is the default.

## Value

A list with NCA control parameters.

## Author

Zhonghui Huang

## Examples

``` r
nca_control()
#> $trapezoidal.rule
#> [1] "linear_up_log_down" "linear"            
#> 
#> $duration
#> NULL
#> 
#> $nlastpoints
#> [1] 3
#> 
#> $slope.method
#> [1] "bestfitforce"
#> 
```
