# Run and evaluate a three-compartment oral model

Fits a three-compartment oral pharmacokinetic model using a naive pooled
data approach and evaluates model performance using prediction error
metrics.

## Usage

``` r
run_npd_3cmpt_oral(
  dat,
  est.method = "nls",
  input.ka = exp(1),
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

  A data frame containing pharmacokinetic data in the standard nlmixr2
  format, including required columns such as `ID`, `EVID`, `DV`, and
  `dose`.

- est.method:

  Estimation method used in nlmixr2. Defaults to "nls".

- input.ka:

  Initial estimate for the absorption rate constant (Ka). Defaults to
  exp(1), corresponding to a log-scale value of 1.

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

  Initial estimate for intercompartmental clearance Q1. Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.q23cmpt:

  Initial estimate for intercompartmental clearance Q2. Defaults to
  exp(1), corresponding to a log-scale value of 1.

- input.add:

  Additive error term. Defaults to 1.

## Value

A list containing parameter estimates and prediction error metrics.

## Details

Rows with `EVID == 2` are excluded prior to model fitting. The model is
fitted using `Fit_3cmpt_oral`, and prediction-based metrics are
calculated to assess model performance.

## See also

[`Fit_3cmpt_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_oral.md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
run_npd_3cmpt_oral(
  dat = Oral_2CPT,
  input.cl = 4,
  input.vc3cmpt = 70,
  input.vp3cmpt = 40,
  input.vp23cmpt = 40,
  input.q3cmpt = 4,
  input.q23cmpt = 4
)
#> exit of dop853 at x = 0.0000000000000000e+00, more than nmax = 70000 are needed
#> exit of dop853 at x = 0.0000000000000000e+00, more than nmax = 70000 are needed
#> exit of dop853 at x = 0.0000000000000000e+00, more than nmax = 70000 are needed
#> exit of dop853 at x = 0.0000000000000000e+00, more than nmax = 70000 are needed
#> $npd.3cmpt_results
#>     ka   cl vc  vp vp2        q       q2    timespent
#> 1 1.79 3.87 96 Inf   0 1.03e-15 1.47e+15 42.8055 secs
#> 
#> $npd.3cmpt.APE
#> metrics.ape 
#>           0 
#> 
#> $npd.3cmpt.MAE
#> metrics.mae 
#>         NaN 
#> 
#> $npd.3cmpt.MAPE
#> metrics.mape 
#>          NaN 
#> 
#> $npd.3cmpt.RMSE
#> metrics.rmse 
#>          NaN 
#> 
#> $npd.3cmpt.rRMSE
#> metrics.rrmse1 
#>            NaN 
#> 
#> $npd.list.3cmpt
#> ── nlmixr² nls with LM algorithm ──
#> 
#>     OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> Pop    0 12807.62 12862.41      -6395.812               1               1
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.021461 1.177    0.001 41.55454
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>           Est.        SE      %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tka     0.5822 9.491e+07  1.63e+10           1.79 (0, inf)                    
#> tcl      1.352 9.491e+07 7.019e+09          3.866 (0, inf)                    
#> tv1      4.564 9.491e+07 2.079e+09          95.98 (0, inf)                    
#> tv2       1318 9.491e+07   7.2e+06            inf (0, inf)                    
#> tv3      -1312 9.491e+07 7.234e+06              0 (0, inf)                    
#> tq      -34.51 9.491e+07  2.75e+08      1.032e-15 (0, inf)                    
#> tq2      34.92 9.491e+07 2.718e+08      1.469e+15 (0, inf)                    
#> add.err      0                                           0                    
#>  
#>   Covariance Type ($covMethod): r (LM)
#>   Information about run found ($runInfo):
#>    • Problems solving ipred with liblsoda, dop853, returning results from the first method 
#>    • some ID(s) could not solve the ODEs correctly; These values are replaced with 'NA' 
#>    • non-positive definite individual Hessian at solution(ID=119); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=118); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=117); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=116); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=115); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=114); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=113); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=112); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=111); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=110); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=109); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=108); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=107); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=106); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=105); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=104); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=103); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=102); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=101); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=100); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=99); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=98); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=97); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=96); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=95); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=94); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=93); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=92); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=91); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=90); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=89); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=88); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=87); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=86); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=85); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=84); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=83); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=82); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=81); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=80); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=79); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=78); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=77); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=76); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=75); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=74); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=73); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=72); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=71); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=70); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=69); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=68); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=67); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=66); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=65); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=64); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=63); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=62); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=61); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=60); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=59); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=58); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=57); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=56); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=55); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=54); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=53); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=52); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=51); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=50); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=49); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=48); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=47); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=46); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=45); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=44); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=43); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=42); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=41); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=40); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=39); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=38); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=37); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=36); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=35); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=34); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=33); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=32); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=31); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=30); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=29); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=28); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=27); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=26); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=25); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=24); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=23); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=22); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=21); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=20); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=19); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=18); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=17); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=16); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=15); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=14); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=13); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=12); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=11); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=10); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=9); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=8); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=7); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=6); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=5); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=4); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=3); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=2); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=1); FOCEi objective functions may not be comparable 
#>    • non-positive definite individual Hessian at solution(ID=0); FOCEi objective functions may not be comparable 
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     The cosine of the angle between `fvec' and any column of the Jacobian is at most `gtol' in absolute value. 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 6,960 × 25
#>   ID     TIME    DV IPRED  IRES IWRES    cp  depot    A1    A2    A3    ka    cl
#>   <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 1      0.25  196.    NA    NA   NaN   NaN 38354.   NaN   NaN   NaN  1.79  3.87
#> 2 1      0.5   310.    NA    NA   NaN   NaN 24516.   NaN   NaN   NaN  1.79  3.87
#> 3 1      0.75  640.    NA    NA   NaN   NaN 15671.   NaN   NaN   NaN  1.79  3.87
#> # ℹ 6,957 more rows
#> # ℹ 12 more variables: v1 <dbl>, v2 <dbl>, v3 <dbl>, q1 <dbl>, q2 <dbl>,
#> #   k <dbl>, k12 <dbl>, k21 <dbl>, k13 <dbl>, k31 <dbl>, tad <dbl>,
#> #   dosenum <dbl>
#> 
# }
```
