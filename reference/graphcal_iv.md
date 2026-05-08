# Graphical calculation of clearance and volume of distribution (IV route)

Estimates clearance (CL), volume of distribution (Vd), terminal slope
(lambdaz), and extrapolated concentration at time zero (C0exp) from
intravenous pharmacokinetic data using graphical methods.

## Usage

``` r
graphcal_iv(dat, dose = 1, ...)
```

## Arguments

- dat:

  A data frame containing TIME (time after dosing) and DV (observed
  concentration).

- dose:

  Administered dose amount. Defaults to 1.

- ...:

  Additional arguments passed to
  [`force_find_lambdaz()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/force_find_lambdaz.md).

## Value

A list containing graphical estimates of CL, Vd, lambda_z, and C0exp.

## Details

Terminal slope (lambdaz) is estimated using
[`force_find_lambdaz()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/force_find_lambdaz.md),
which applies an automated phase selection strategy with fallback
regression when required.

## See also

[`force_find_lambdaz`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/force_find_lambdaz.md)

## Author

Zhonghui Huang

## Examples

``` r
dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
                  DV = c(12, 8, 5, 3, 2, 1.5, 1))
graphcal_iv(dat, dose = 100)
#> $cl
#> [1] 2.951299
#> 
#> $vd
#> [1] 16.47151
#> 
#> $slope
#> [1] -0.1791759
#> 
#> $C0exp
#> [1] 6.071088
#> 
#> $method
#> [1] "find_best_lambdaz"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>       1.8035       -0.1792  
#> 
#> 
#> $time.spent
#> [1] 0.003
#> 
```
