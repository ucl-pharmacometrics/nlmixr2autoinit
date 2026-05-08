# Run and evaluate a one-compartment oral model

Fits a one-compartment oral pharmacokinetic model using a naive pooled
data approach and evaluates model performance based on prediction error
metrics.

## Usage

``` r
run_npd_1cmpt_oral(
  dat,
  est.method = "nls",
  input.ka = exp(1),
  input.cl = exp(1),
  input.vd = exp(1),
  input.add = 1
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in standard
  nlmixr2 format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.ka:

  Initial estimate for the absorption rate constant (Ka). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.cl:

  Initial estimate for clearance (CL). Defaults to exp(1), corresponding
  to a log-scale value of 1.

- input.vd:

  Initial estimate for volume of distribution (Vd). Defaults to exp(1),
  corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

## Value

A list containing the fitted parameter estimates and prediction error
metrics.

## Details

Rows with `EVID == 2` are excluded before model fitting. The model is
fitted using `Fit_1cmpt_oral`, and prediction-based metrics are
calculated to assess performance.

## See also

[`Fit_1cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_1cmpt_oral(dat = Oral_1CPT, input.ka = 1, input.cl = 4, input.vd = 70)
#> $npd.1cmpt_results
#>     ka   cl   vd   timespent
#> 1 1.02 3.89 71.2 4.2775 secs
#> 
#> $npd.1cmpt.APE
#> metrics.ape 
#>    929309.4 
#> 
#> $npd.1cmpt.MAE
#> metrics.mae 
#>       133.8 
#> 
#> $npd.1cmpt.MAPE
#> metrics.mape 
#>         84.7 
#> 
#> $npd.1cmpt.RMSE
#> metrics.rmse 
#>        230.1 
#> 
#> $npd.1cmpt.rRMSE
#> metrics.rrmse1 
#>           45.6 
#> 
#> $nnpd.1cmpt.list
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 367749941 367762717 367762744     -183881354        26.63715
#>     Condition#(Cor)
#> Pop        5.628584
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.024411 0.162    0.001 4.067589
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>            Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tka     0.01948 0.0002116    1.087      1.02 (1.019, 1.02)                    
#> tcl       1.358   4.8e-05 0.003533       3.89 (3.89, 3.89)                    
#> tv        4.266 9.224e-05 0.002162    71.24 (71.23, 71.25)                    
#> add.err   230.1                                      230.1                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,947 × 15
#>   ID     TIME    DV IPRED  IRES   IWRES    cp  depot centre    ka    cl     v
#>   <fct> <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25  205.  188.  16.6  0.0723  188. 46499. 13406.  1.02  3.89  71.2
#> 2 1      0.5   311.  331. -20.9 -0.0906  331. 36036. 23613.  1.02  3.89  71.2
#> 3 1      0.75  389.  440. -50.8 -0.221   440. 27927. 31344.  1.02  3.89  71.2
#> # ℹ 6,944 more rows
#> # ℹ 3 more variables: k <dbl>, tad <dbl>, dosenum <dbl>
#> 
 # }
```
