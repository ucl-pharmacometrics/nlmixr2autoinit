# Run and evaluate a two-compartment oral model

Fits a two-compartment oral pharmacokinetic model using a naive pooled
data approach and evaluates model performance using prediction error
metrics.

## Usage

``` r
run_npd_2cmpt_oral(
  dat,
  est.method = "nls",
  input.ka = exp(1),
  input.cl = exp(1),
  input.vc2cmpt = exp(1),
  input.vp2cmpt = exp(1),
  input.q2cmpt = exp(1),
  input.add = 1
)
```

## Arguments

- dat:

  A data frame containing time–concentration data in standard nlmixr2
  format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.ka:

  Initial estimate for the absorption rate constant (Ka). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.cl:

  Initial estimate for clearance (CL). Defaults to exp(1) ,
  corresponding to a log-scale value of 1.

- input.vc2cmpt:

  Initial estimate for the central volume of distribution (Vc). Defaults
  to exp(1), corresponding to a log-scale value of 1.

- input.vp2cmpt:

  Initial estimate for the peripheral volume of distribution (Vp).
  Defaults to exp(1), corresponding to a log-scale value of 1.

- input.q2cmpt:

  Initial estimate for intercompartmental clearance (Q). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

## Value

A list containing the fitted parameter estimates and prediction error
metrics.

## Details

Rows with `EVID == 2` are excluded before fitting the model. The model
is fitted using `Fit_2cmpt_oral`, and prediction-based metrics are
computed to evaluate performance.

## See also

[`Fit_2cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_2cmpt_oral(
  dat = Oral_2CPT,
  input.ka = 1,
  input.cl = 4,
  input.vc2cmpt = 35,
  input.vp2cmpt = 35,
  input.q2cmpt = 4
)
#> $npd.2cmpt_results
#>     ka   cl vc vp    q    timespent
#> 1 1.06 3.88 72 51 3.29 10.5266 secs
#> 
#> $npd.2cmpt.APE
#> metrics.ape 
#>      940656 
#> 
#> $npd.2cmpt.MAE
#> metrics.mae 
#>       135.2 
#> 
#> $npd.2cmpt.MAPE
#> metrics.mape 
#>         52.8 
#> 
#> $npd.2cmpt.RMSE
#> metrics.rmse 
#>        231.2 
#> 
#> $npd.2cmpt.rRMSE
#> metrics.rrmse1 
#>           45.9 
#> 
#> $npd.list.2cmpt
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 371973525 371986329 371986370     -185993159        720.4805
#>     Condition#(Cor)
#> Pop        32.34632
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.018555 0.169    0.001 10.28944
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>           Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tka     0.0549 0.0003759   0.6848    1.056 (1.056, 1.057)                    
#> tcl      1.355 5.931e-05 0.004379    3.875 (3.875, 3.875)                    
#> tv1      4.277 0.0001987 0.004645    72.02 (71.99, 72.05)                    
#> tv2      3.931 0.0006644   0.0169    50.96 (50.89, 51.02)                    
#> tq       1.191  0.001148  0.09642     3.29 (3.283, 3.298)                    
#> add.err  231.1                                      231.1                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,960 × 20
#>   ID     TIME    DV IPRED   IRES   IWRES    cp  depot     A1    A2    ka    cl
#>   <fct> <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25  196.  191.   5.21  0.0225  191. 46073. 13748.  81.9  1.06  3.88
#> 2 1      0.5   310.  333. -22.7  -0.0983  333. 35379. 23969. 297.   1.06  3.88
#> 3 1      0.75  640.  437. 202.    0.876   437. 27167. 31494. 609.   1.06  3.88
#> # ℹ 6,957 more rows
#> # ℹ 8 more variables: v1 <dbl>, v2 <dbl>, q <dbl>, k <dbl>, k12 <dbl>,
#> #   k21 <dbl>, tad <dbl>, dosenum <dbl>
#> 
# }
```
