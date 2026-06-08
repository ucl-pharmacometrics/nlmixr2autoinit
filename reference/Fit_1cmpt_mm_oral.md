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
  ncores = 2,
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

- ncores:

  Number of cores to use for parallelization, passed to `rxControl()`.
  Default is 2.

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
# Run simulation
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
#> 1     1     0  0.25     1 1002.  250.  70.1  31.3     31.3  29.0 7788.  2197.
#> 2     1     0  0.5      1 1002.  250.  70.1  55.4     55.4  59.0 6065.  3882.
#> 3     1     0  0.75     1 1002.  250.  70.1  73.8     73.8  87.3 4724.  5172.
#> 4     1     0  1        1 1002.  250.  70.1  87.8     87.8  68.1 3679.  6156.
#> 5     1     0  1.5      1 1002.  250.  70.1 106.     106.  102.  2231.  7462.
#> 6     1     0  2        1 1002.  250.  70.1 117.     117.  117.  1353.  8185.
#> # ℹ 6,954 more rows
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
#> 1     1     0  0.25     1 1002.  250.  70.1  31.3     31.3  37.2 7788.  2197.
#> 2     1     0  0.5      1 1002.  250.  70.1  55.4     55.4  58.7 6065.  3882.
#> 3     1     0  0.75     1 1002.  250.  70.1  73.8     73.8  68.5 4724.  5172.
#> 4     1     0  1        1 1002.  250.  70.1  87.8     87.8  97.2 3679.  6156.
#> 5     1     0  1.5      1 1002.  250.  70.1 106.     106.   93.2 2231.  7462.
#> 6     1     0  2        1 1002.  250.  70.1 117.     117.  124.  1353.  8185.
#> # ℹ 6,954 more rows
# }
```
