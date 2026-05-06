# Linear trapezoidal rule

Computes the area under the curve (AUC) or the area under the moment
curve (AUMC) using the linear trapezoidal rule. If `moment = TRUE`, the
function estimates AUMC by integrating `time * concentration`.

## Usage

``` r
trapezoidal_linear(x, y, moment = FALSE)
```

## Arguments

- x:

  A numeric vector representing the time points.

- y:

  A numeric vector representing the corresponding concentration values
  at each time point.

- moment:

  Logical. If `TRUE`, computes AUMC by integrating `t * C(t)` instead of
  just `C(t)`.

## Value

A numeric value representing the estimated AUC or AUMC using the linear
trapezoidal rule.

## Examples

``` r
x <- c(0.5, 1, 2, 4, 6, 8)
y <- c(12, 8, 5, 3, 2, 1)
trapezoidal_linear(x, y)                # AUC
#> [1] 27.5
trapezoidal_linear(x, y, moment = TRUE) # AUMC
#> [1] 78.5
```
