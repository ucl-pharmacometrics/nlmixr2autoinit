# Run and evaluate a one-compartment oral model with Michaelis-Menten kinetics

Fits a one-compartment oral pharmacokinetic model with Michaelis-Menten
elimination using a naive pooled data approach, and evaluates model
performance using prediction error metrics.

## Usage

``` r
run_npd_1cmpt_mm_oral(
  dat,
  est.method = "nls",
  input.ka = exp(1),
  input.vmax = exp(1),
  input.km = exp(1),
  input.cl = exp(1),
  input.vd = exp(1),
  input.add = 1,
  km_threshold = FALSE
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in standard
  nlmixr2 format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.ka:

  Initial estimate for the absorption rate constant (ka). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.vmax:

  Initial estimate for the maximum metabolic rate (Vmax). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.km:

  Initial estimate for the Michaelis constant (Km). Defaults to exp(1),
  corresponding to a log-scale value of 1.

- input.cl:

  Initial estimate for clearance (CL). Defaults to exp(1), corresponding
  to a log-scale value of 1. This value is used to derive initial Vmax
  and Km when `km_threshold = TRUE`.

- input.vd:

  Initial estimate for volume of distribution (Vd). Defaults to exp(1),
  corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

- km_threshold:

  Logical indicating whether initial Vmax and Km should be automatically
  adjusted based on observed maximum concentration and clearance.
  Defaults to FALSE.

## Value

A list containing the fitted parameter estimates and prediction error
metrics.

## Details

The function excludes dosing records (`EVID == 2`) prior to model
fitting. When `km_threshold = TRUE`, initial estimates for Vmax and Km
are derived using the observed maximum concentration and clearance. The
model is then fitted using `Fit_1cmpt_mm_oral`, and prediction-based
metrics are calculated to assess performance.

## See also

[`Fit_1cmpt_mm_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
  run_npd_1cmpt_mm_oral(
    dat = Oral_1CPTMM,
    input.ka = 1,
    input.vmax = 1000,
    input.km = 250,
    input.vd = 70
  )
#> $npd.1cmpt.mm_results
#>     ka vmax  km   vd timespent
#> 1 0.87 1050 352 66.1 2.48 secs
#> 
#> $npd.1cmpt.mm.APE
#> metrics.ape 
#>     3043122 
#> 
#> $npd.1cmpt.mm.MAE
#> metrics.mae 
#>       437.3 
#> 
#> $npd.1cmpt.mm.MAPE
#> metrics.mape 
#>        119.4 
#> 
#> $npd.1cmpt.mm.RMSE
#> metrics.rmse 
#>        889.3 
#> 
#> $npd.1cmpt.mm.rRMSE
#> metrics.rrmse1 
#>           63.1 
#> 
#> $npd.1cmpt.mm.list
#> ── nlmixr² nls with LM algorithm ──
#> 
#>           OBJF        AIC        BIC Log-likelihood Condition#(Cov)
#> Pop 5503432579 5503445379 5503445413    -2751722684        1496.186
#>     Condition#(Cor)
#> Pop        203.2179
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table    other
#> elapsed 0.022884 0.136 2.292116
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>            Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tka     -0.1394  0.000243   0.1743 0.8699 (0.8695, 0.8703)                    
#> lvmax     6.956 0.0001296 0.001863       1050 (1050, 1050)                    
#> lkm       5.864 0.0003416 0.005825      352 (351.8, 352.2)                    
#> tv        4.192 4.439e-05 0.001059    66.15 (66.14, 66.15)                    
#> add.err   889.4                                      889.4                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,959 × 15
#>   ID     TIME    DV IPRED   IRES    IWRES    cp depot centre    ka  vmax    km
#>   <fct> <dbl> <dbl> <dbl>  <dbl>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25  27.4  29.4  -1.99 -0.00223  29.4 8045.  1944. 0.870 1050.  352.
#> 2 1      0.5   32.1  52.7 -20.6  -0.0232   52.7 6473.  3489. 0.870 1050.  352.
#> 3 1      0.75  93.5  71.3  22.2   0.0250   71.3 5208.  4714. 0.870 1050.  352.
#> # ℹ 6,956 more rows
#> # ℹ 3 more variables: v <dbl>, tad <dbl>, dosenum <dbl>
#> 
# }
```
