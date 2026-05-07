# Run and evaluate a one-compartment IV model

Fits a one-compartment intravenous pharmacokinetic model using a naive
pooled data approach and evaluates model performance based on prediction
error metrics.

## Usage

``` r
run_npd_1cmpt_iv(
  dat,
  est.method = "nls",
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

Rows with `EVID == 2` are excluded prior to model fitting. The model is
fitted using `Fit_1cmpt_iv`, and prediction-based metrics are calculated
to assess performance.

## See also

[`Fit_1cmpt_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_iv.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_1cmpt_iv(dat = Bolus_1CPT, input.cl = 4, input.vd = 70)
#> $npd.1cmpt_results
#>     cl   vd   timespent
#> 1 3.82 63.5 1.7981 secs
#> 
#> $npd.1cmpt.APE
#> metrics.ape 
#>     1096334 
#> 
#> $npd.1cmpt.MAE
#> metrics.mae 
#>       157.7 
#> 
#> $npd.1cmpt.MAPE
#> metrics.mape 
#>         65.5 
#> 
#> $npd.1cmpt.RMSE
#> metrics.rmse 
#>        273.2 
#> 
#> $npd.1cmpt.rRMSE
#> metrics.rrmse1 
#>           44.2 
#> 
#> $nnpd.1cmpt.list
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 518958748 518971529 518971550     -259485762        2.289302
#>     Condition#(Cor)
#> Pop        2.244467
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.020958 0.116    0.001 1.639042
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>          Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tcl      1.34 4.631e-05 0.003457    3.818 (3.817, 3.818)                    
#> tv       4.15 5.055e-05 0.001218    63.45 (63.45, 63.46)                    
#> add.err 273.2                                      273.2                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,951 × 13
#>   ID     TIME    DV IPRED  IRES   IWRES    cp centre    cl     v      k   tad
#>   <fct> <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>
#> 1 1      0.25 1126.  931. 195.   0.712   931. 59104.  3.82  63.5 0.0602  0.25
#> 2 1      0.5   870.  918. -47.7 -0.174   918. 58222.  3.82  63.5 0.0602  0.5 
#> 3 1      0.75  884.  904. -20.3 -0.0741  904. 57353.  3.82  63.5 0.0602  0.75
#> # ℹ 6,948 more rows
#> # ℹ 1 more variable: dosenum <dbl>
#> 
# }
```
