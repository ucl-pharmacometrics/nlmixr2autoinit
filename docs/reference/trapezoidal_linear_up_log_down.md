# Linear-up and log-down trapezoidal rule

Computes the area under the curve (AUC) or the area under the moment
curve (AUMC) using a hybrid trapezoidal rule. The method uses linear
interpolation for increasing or constant concentration segments, and
logarithmic interpolation for decreasing segments.

## Usage

``` r
trapezoidal_linear_up_log_down(x, y, moment = FALSE)
```

## Arguments

- x:

  A numeric vector representing the time points.

- y:

  A numeric vector representing the corresponding concentration values
  at each time point.

- moment:

  Logical. If `TRUE`, computes AUMC by integrating `t * C(t)` instead of
  `C(t)`.

## Value

A numeric value representing the estimated AUC or AUMC using the
linear-up/log-down trapezoidal method.

## Details

If `moment = TRUE`, the function calculates the area under the moment
curve (AUMC), i.e., it integrates `t * C(t)` over time instead of just
`C(t)`.

## Examples

``` r
x <- c(0, 0.5, 1, 2, 4, 6, 8)
y <- c(0, 2, 8, 5, 3, 2, 1)
trapezoidal_linear_up_log_down(x, y)                # AUC
#> [1] 25.03139
trapezoidal_linear_up_log_down(x, y, moment = TRUE) # AUMC
#> [1] 78.85055
```
