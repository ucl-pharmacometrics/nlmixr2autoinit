# Run and evaluate a one-compartment IV Michaelis-Menten model

Fits a one-compartment intravenous pharmacokinetic model with
Michaelis-Menten elimination using a naive pooled data approach and
evaluates model performance based on prediction error metrics.

## Usage

``` r
run_npd_1cmpt_mm_iv(
  dat,
  est.method = "nls",
  npdmm_inputvmax = exp(1),
  npdmm_inputkm = exp(1),
  npdmm_inputcl = exp(1),
  npdmm_inputvd = exp(1),
  input.add = 1,
  km_threshold = FALSE
)
```

## Arguments

- dat:

  A data frame containing pharmacokinetic data in standard nlmixr2
  format.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- npdmm_inputvmax:

  Initial estimate for Vmax. Defaults to exp(1), corresponding to a
  log-scale value of 1.

- npdmm_inputkm:

  Initial estimate for Km. Defaults to exp(1), corresponding to a
  log-scale value of 1.

- npdmm_inputcl:

  Initial estimate for clearance (CL). Defaults to exp(1) ,
  corresponding to a log-scale value of 1.

- npdmm_inputvd:

  Initial estimate for volume of distribution (Vd). Defaults to exp(1),
  corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

- km_threshold:

  Logical value. If TRUE, initial estimates for Vmax and Km are
  calculated based on the maximum observed concentration.

## Value

A list containing parameter estimates and prediction error metrics.

## Details

Rows where `EVID == 2` are excluded before model fitting. The model is
fitted using `Fit_1cmpt_mm_iv`. When `km_threshold = TRUE`, initial
estimates for Vmax and Km are derived from the dataset to provide a
representative starting point for nonlinear elimination.

## See also

[`Fit_1cmpt_mm_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_iv.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_1cmpt_mm_iv(
  dat = Bolus_1CPT,
  npdmm_inputcl = 4,
  npdmm_inputvd = 70,
  km_threshold = TRUE
)
#> $npd.1cmpt.mm_results
#>       vmax       km   vd timespent
#> 1 8.45e+10 2.21e+10 63.5 3.69 secs
#> 
#> $npd.1cmpt.mm.APE
#> metrics.ape 
#>     1096334 
#> 
#> $npd.1cmpt.mm.MAE
#> metrics.mae 
#>       157.7 
#> 
#> $npd.1cmpt.mm.MAPE
#> metrics.mape 
#>         65.5 
#> 
#> $npd.1cmpt.mm.RMSE
#> metrics.rmse 
#>        273.2 
#> 
#> $npd.1cmpt.mm.rRMSE
#> metrics.rrmse1 
#>           44.2 
#> 
#> $npd.1cmpt.mm.list
#> ── nlmixr² nls with LM algorithm ──
#> 
#>          OBJF       AIC       BIC Log-likelihood Condition#(Cov)
#> Pop 518958749 518971532 518971560     -259485762        7.214892
#>     Condition#(Cor)
#> Pop        6.913446
#> 
#> ── Time (sec $time): ──
#> 
#>           setup table compress   other
#> elapsed 0.02029 0.109    0.001 3.51371
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>          Est.        SE      %RSE          Back-transformed(95%CI) BSV(SD)
#> lvmax   25.16 5.504e-05 0.0002188  8.449e+10 (8.448e+10, 8.45e+10)        
#> lkm     23.82 5.504e-05 0.0002311 2.213e+10 (2.213e+10, 2.213e+10)        
#> tv       4.15 4.105e-05 0.0009891             63.45 (63.45, 63.46)        
#> add.err 273.2                                                273.2        
#>         Shrink(SD)%
#> lvmax              
#> lkm                
#> tv                 
#> add.err            
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,951 × 13
#>   ID     TIME    DV IPRED  IRES   IWRES    cp centre    vmax      km     v   tad
#>   <fct> <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>   <dbl>   <dbl> <dbl> <dbl>
#> 1 1      0.25 1126.  931. 195.   0.712   931. 59104. 8.45e10 2.21e10  63.5  0.25
#> 2 1      0.5   870.  918. -47.7 -0.174   918. 58222. 8.45e10 2.21e10  63.5  0.5 
#> 3 1      0.75  884.  904. -20.3 -0.0741  904. 57353. 8.45e10 2.21e10  63.5  0.75
#> # ℹ 6,948 more rows
#> # ℹ 1 more variable: dosenum <dbl>
#> 
# }
```
