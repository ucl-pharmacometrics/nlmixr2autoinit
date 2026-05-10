# Fit oral pharmacokinetic data to a one-compartment model with Michaelis-Menten elimination

Fits oral pharmacokinetic data to a one-compartment model with
first-order absorption and Michaelis-Menten (nonlinear) elimination
using the naive pooled data approach. Supports multiple estimation
methods available in nlmixr2, and optionally returns only predicted
concentrations to reduce memory use in simulation workflows.

## Usage

``` r
Fit_1cmpt_mm_oral(
  data,
  est.method,
  input.ka,
  input.vmax,
  input.km,
  input.vd,
  input.add,
  return.pred.only = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame of oral pharmacokinetic data formatted for nlmixr2.

- est.method:

  Estimation method to use in nlmixr2, one of: `"rxSolve"`, `"nls"`,
  `"nlm"`, `"nlminb"`, or `"focei"`.

- input.ka:

  Initial estimate of the absorption rate constant (ka).

- input.vmax:

  Initial estimate of the maximum elimination rate (Vmax).

- input.km:

  Initial estimate of the Michaelis constant (Km).

- input.vd:

  Initial estimate of the volume of distribution (V).

- input.add:

  Initial estimate of the additive residual error.

- return.pred.only:

  Logical; if `TRUE`, returns a data frame with only predicted
  concentrations (`cp`) for all observations in the input data.

- ...:

  Optional arguments passed to `nlmixr2()`, such as a custom
  `control = foceiControl(...)` or other control objects.

## Value

If `return.pred.only = TRUE`, returns a `data.frame` with columns `cp`
(predicted concentration) and `DV` (observed data). Otherwise, returns a
fitted model object produced by nlmixr2.

## Author

Zhonghui Huang

## Examples

``` r
 # \donttest{
dat <- Oral_1CPTMM
# Fit using 'nls'
Fit_1cmpt_mm_oral(
  data = dat,
  est.method = "nls",
  input.ka = 1,
  input.vmax = 1000,
  input.km = 250,
  input.vd = 70,
  input.add = 10
)
#> ── nlmixr² nls with LM algorithm ──
#> 
#>           OBJF        AIC        BIC Log-likelihood Condition#(Cov)
#> Pop 5503432579 5503445379 5503445413    -2751722684        1496.186
#>     Condition#(Cor)
#> Pop        203.2179
#> 
#> ── Time (sec value$time): ──
#> 
#>            setup table compress    other
#> elapsed 0.023565 0.131    0.002 3.767435
#> 
#> ── (value$parFixed or value$parFixedDf): ──
#> 
#>            Est.        SE     %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> tka     -0.1394  0.000243   0.1743 0.8699 (0.8695, 0.8703)                    
#> lvmax     6.956 0.0001296 0.001863       1050 (1050, 1050)                    
#> lkm       5.864 0.0003416 0.005825      352 (351.8, 352.2)                    
#> tv        4.192 4.439e-05 0.001059    66.15 (66.14, 66.15)                    
#> add.err   889.4                                      889.4                    
#>  
#>   Covariance Type (value$covMethod): r (LM)
#>   Censoring (value$censInformation): No censoring
#>   Minimization message (value$message):  
#>     Relative error in the sum of squares is at most `ftol'. 
#> 
#> ── Fit Data (object value is a modified tibble): ──
#> # A tibble: 6,960 × 16
#>   ID     EVID  TIME    DV IPRED   IRES    IWRES    cp depot centre    ka  vmax
#>   <fct> <int> <dbl> <dbl> <dbl>  <dbl>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
#> 1 1         0  0.25  27.4  29.4  -1.99 -0.00223  29.4 8045.  1944. 0.870 1050.
#> 2 1         0  0.5   32.1  52.7 -20.6  -0.0232   52.7 6473.  3489. 0.870 1050.
#> 3 1         0  0.75  93.5  71.3  22.2   0.0250   71.3 5208.  4714. 0.870 1050.
#> # ℹ 6,957 more rows
#> # ℹ 4 more variables: km <dbl>, v <dbl>, tad <dbl>, dosenum <dbl>
# Return only predicted concentrations
Fit_1cmpt_mm_oral(
  data = dat,
  est.method = "rxSolve",
  input.ka = 1,
  input.vmax = 1000,
  input.km = 250,
  input.vd = 70,
  input.add = 10
)
#> ── Solved rxode2 object ──
#> ── Parameters (value$params): ──
#> # A tibble: 120 × 6
#>    id      tka lvmax   lkm    tv add.err
#>    <fct> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1 1         0  6.91  5.52  4.25      10
#>  2 2         0  6.91  5.52  4.25      10
#>  3 3         0  6.91  5.52  4.25      10
#>  4 4         0  6.91  5.52  4.25      10
#>  5 5         0  6.91  5.52  4.25      10
#>  6 6         0  6.91  5.52  4.25      10
#>  7 7         0  6.91  5.52  4.25      10
#>  8 8         0  6.91  5.52  4.25      10
#>  9 9         0  6.91  5.52  4.25      10
#> 10 10        0  6.91  5.52  4.25      10
#> # ℹ 110 more rows
#> ── Initial Conditions (value$inits): ──
#>  depot centre 
#>      0      0 
#> ── First part of data (object): ──
#> # A tibble: 6,960 × 12
#>      id  evid  time    ka  vmax    km     v    cp ipredSim   sim depot centre
#>   <int> <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl> <dbl>  <dbl>
#> 1     1     0  0.25     1 1002.  250.  70.1  31.3     31.3  44.3 7788.  2197.
#> 2     1     0  0.5      1 1002.  250.  70.1  55.4     55.4  51.6 6065.  3882.
#> 3     1     0  0.75     1 1002.  250.  70.1  73.8     73.8  88.8 4724.  5172.
#> 4     1     0  1        1 1002.  250.  70.1  87.8     87.8  82.1 3679.  6156.
#> 5     1     0  1.5      1 1002.  250.  70.1 106.     106.  104.  2231.  7462.
#> 6     1     0  2        1 1002.  250.  70.1 117.     117.  114.  1353.  8185.
#> # ℹ 6,954 more rows
# }
```
