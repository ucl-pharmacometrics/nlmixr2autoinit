# Run and evaluate a three-compartment IV model

Fits a three-compartment intravenous pharmacokinetic model using a naive
pooled data approach and evaluates model performance based on prediction
error metrics.

## Usage

``` r
run_npd_3cmpt_iv(
  dat,
  est.method = "nls",
  input.cl = exp(1),
  input.vc3cmpt = exp(1),
  input.vp3cmpt = exp(1),
  input.vp23cmpt = exp(1),
  input.q3cmpt = exp(1),
  input.q23cmpt = exp(1),
  input.add = 1
)
```

## Arguments

- dat:

  A data frame containing raw intravenous concentration–time data in
  standard nlmixr2 format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.cl:

  Initial estimate for clearance (CL). Defaults to exp(1), corresponding
  to a log-scale value of 1.

- input.vc3cmpt:

  Initial estimate for central volume of distribution (Vc). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.vp3cmpt:

  Initial estimate for first peripheral volume (Vp1). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.vp23cmpt:

  Initial estimate for second peripheral volume (Vp2). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.q3cmpt:

  Initial estimate for intercompartmental clearance between central and
  first peripheral compartments (Q1). Defaults to exp(1), corresponding
  to a log-scale value of 1.

- input.q23cmpt:

  Initial estimate for intercompartmental clearance between central and
  second peripheral compartments (Q2). Defaults to exp(1), corresponding
  to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

## Value

A list containing fitted parameter estimates and model prediction error
metrics.

## Details

Rows with `EVID == 2` are excluded prior to model fitting. The model is
fitted using `Fit_3cmpt_iv`, and prediction-based metrics are calculated
to evaluate performance.

## See also

[`Fit_3cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_iv.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_3cmpt_iv(
  dat = Bolus_2CPT,
  input.cl = 4,
  input.vc3cmpt = 70,
  input.vp3cmpt = 40,
  input.vp23cmpt = 10,
  input.q3cmpt = 4,
  input.q23cmpt = 4
)
#> $npd.3cmpt_results
#>     cl   vc   vp  vp2    q   q2    timespent
#> 1 3.79 62.6 60.2 19.3 1.88 3.51 40.6602 secs
#> 
#> $npd.3cmpt.APE
#> metrics.ape 
#>     1095347 
#> 
#> $npd.3cmpt.MAE
#> metrics.mae 
#>       157.4 
#> 
#> $npd.3cmpt.MAPE
#> metrics.mape 
#>         40.3 
#> 
#> $npd.3cmpt.RMSE
#> metrics.rmse 
#>        271.6 
#> 
#> $npd.3cmpt.rRMSE
#> metrics.rrmse1 
#>           45.4 
#> 
#> $npd.list.3cmpt
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 513582700 513595506 513595554     -256797746        32487.15
#>     Condition#(Cor)
#> Pop        231.5548
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.018967   0.1    0.001 40.51203
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>           Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tcl      1.333  6.25e-05 0.004689    3.791 (3.791, 3.792)                    
#> tv1      4.136 8.753e-05 0.002116    62.55 (62.54, 62.56)                    
#> tv2      4.097  0.001449  0.03536    60.16 (59.99, 60.33)                    
#> tv3      2.958  0.005456   0.1845    19.25 (19.05, 19.46)                    
#> tq      0.6289  0.004597   0.7309    1.876 (1.859, 1.893)                    
#> tq2      1.256  0.002116   0.1685    3.511 (3.496, 3.526)                    
#> add.err  271.6                                      271.6                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,960 × 23
#>   ID     TIME    DV IPRED  IRES  IWRES    cp  centre    A2    A3    cl    v1
#>   <fct> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25 1041. 1850. -809. -2.98  1850. 115718.  880. 1616.  3.79  62.6
#> 2 1      0.5  1629  1785. -156. -0.575 1785. 111667. 1722. 3103.  3.79  62.6
#> 3 1      0.75  878. 1724. -846. -3.11  1724. 107832. 2528. 4470.  3.79  62.6
#> # ℹ 6,957 more rows
#> # ℹ 11 more variables: v2 <dbl>, v3 <dbl>, q <dbl>, q2 <dbl>, k <dbl>,
#> #   k12 <dbl>, k21 <dbl>, k13 <dbl>, k31 <dbl>, tad <dbl>, dosenum <dbl>
#> 
# }
```
