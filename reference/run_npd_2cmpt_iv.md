# Run and evaluate a two-compartment IV model

Fits a two-compartment intravenous pharmacokinetic model using a naive
pooled data approach and evaluates model performance based on prediction
error metrics.

## Usage

``` r
run_npd_2cmpt_iv(
  dat,
  est.method = "nls",
  input.cl = exp(1),
  input.vc2cmpt = exp(1),
  input.vp2cmpt = exp(1),
  input.q2cmpt = exp(1),
  input.add = 1
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in standard
  nlmixr2 format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.cl:

  Initial estimate for clearance (CL). Defaults to exp(1), corresponding
  to a log-scale value of 1.

- input.vc2cmpt:

  Initial estimate for central compartment volume (Vc). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.vp2cmpt:

  Initial estimate for peripheral compartment volume (Vp). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.q2cmpt:

  Initial estimate for intercompartmental clearance (Q). Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

## Value

A list containing parameter estimates and prediction error metrics.

## Details

Rows with `EVID == 2` are excluded prior to model fitting. The model is
fitted using `Fit_2cmpt_iv`, and prediction-based metrics are calculated
to assess performance.

## See also

[`Fit_2cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_iv.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_2cmpt_iv(dat = Bolus_2CPT,
                           input.cl = 4,
                           input.vc2cmpt = 35,
                           input.vp2cmpt = 35,
                           input.q2cmpt = 4)
#> $npd.2cmpt_results
#>     cl   vc vp    q   timespent
#> 1 3.84 63.5 71 4.45 4.2397 secs
#> 
#> $npd.2cmpt.APE
#> metrics.ape 
#>     1095258 
#> 
#> $npd.2cmpt.MAE
#> metrics.mae 
#>       157.4 
#> 
#> $npd.2cmpt.MAPE
#> metrics.mape 
#>         40.7 
#> 
#> $npd.2cmpt.RMSE
#> metrics.rmse 
#>        271.8 
#> 
#> $npd.2cmpt.rRMSE
#> metrics.rrmse1 
#>           45.4 
#> 
#> $npd.list.2cmpt
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 514200464 514213265 514213299     -257106628        141.5929
#>     Condition#(Cor)
#> Pop        3.803517
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.036282 0.102    0.001 4.073718
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>          Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tcl     1.345 5.452e-05 0.004055    3.837 (3.836, 3.837)                    
#> tv1     4.151 6.162e-05 0.001484    63.52 (63.52, 63.53)                    
#> tv2     4.262 0.0004843  0.01136    70.97 (70.91, 71.04)                    
#> tq      1.492 0.0004336  0.02906    4.447 (4.443, 4.451)                    
#> add.err 271.8                                      271.8                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,960 × 18
#>   ID     TIME    DV IPRED  IRES  IWRES    cp  centre    A2    cl    v1    v2
#>   <fct> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25 1041. 1829. -788. -2.90  1829. 116167. 2050.  3.84  63.5  71.0
#> 2 1      0.5  1629  1771. -142. -0.522 1771. 112488. 4003.  3.84  63.5  71.0
#> 3 1      0.75  878. 1715. -837. -3.08  1715. 108956. 5863.  3.84  63.5  71.0
#> # ℹ 6,957 more rows
#> # ℹ 6 more variables: q <dbl>, k <dbl>, k12 <dbl>, k21 <dbl>, tad <dbl>,
#> #   dosenum <dbl>
#> 
# }
```
