# Fit oral pharmacokinetic data to a two-compartment model

Fits oral pharmacokinetic data to a two-compartment model with
first-order absorption and first-order elimination using the naive
pooled data approach. Supports multiple estimation methods available in
nlmixr2, and optionally returns only predicted concentrations to support
simulation workflows.

## Usage

``` r
Fit_2cmpt_oral(
  data,
  est.method,
  input.ka,
  input.cl,
  input.vc2cmpt,
  input.vp2cmpt,
  input.q2cmpt,
  input.add,
  ncores = 2,
  return.pred.only = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame containing oral pharmacokinetic data formatted for
  nlmixr2,

- est.method:

  Estimation method to use in nlmixr2, one of: `"rxSolve"`, `"nls"`,
  `"nlm"`, `"nlminb"`, or `"focei"`.

- input.ka:

  Initial estimate of the absorption rate constant (ka).

- input.cl:

  Initial estimate of clearance (CL).

- input.vc2cmpt:

  Initial estimate of central volume of distribution (V1).

- input.vp2cmpt:

  Initial estimate of peripheral volume of distribution (V2).

- input.q2cmpt:

  Initial estimate of inter-compartmental clearance (Q).

- input.add:

  Initial estimate of the additive residual error.

- ncores:

  Number of cores to use for parallelization, passed to `rxControl()`.
  Default is 2.

- return.pred.only:

  Logical; if `TRUE`, returns a data frame with only predicted
  concentrations (`cp`) for all observations in the input data.

- ...:

  Additional arguments passed to `nlmixr2()`, such as a user-defined
  `control = foceiControl(...)` or other control settings.

## Value

If `return.pred.only = TRUE`, returns a `data.frame` with a single
column `cp` (predicted concentrations). Otherwise, returns a fitted
model object produced by nlmixr2.

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
dat <- Oral_2CPT[Oral_2CPT$ID<11,]
# Run simulation
Fit_2cmpt_oral(
  data = dat,
  est.method = "rxSolve",
  input.ka = 1,
  input.cl = 4,
  input.vc2cmpt = 70,
  input.vp2cmpt = 40,
  input.q2cmpt = 10,
  input.add = 10
)
#> ── Solved rxode2 object ──
#> ── Parameters (value$params): ──
#> # A tibble: 10 × 7
#>    id      tka   tcl   tv1   tv2    tq add.err
#>    <fct> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1 1         0  1.39  4.25  3.69   2.3      10
#>  2 2         0  1.39  4.25  3.69   2.3      10
#>  3 3         0  1.39  4.25  3.69   2.3      10
#>  4 4         0  1.39  4.25  3.69   2.3      10
#>  5 5         0  1.39  4.25  3.69   2.3      10
#>  6 6         0  1.39  4.25  3.69   2.3      10
#>  7 7         0  1.39  4.25  3.69   2.3      10
#>  8 8         0  1.39  4.25  3.69   2.3      10
#>  9 9         0  1.39  4.25  3.69   2.3      10
#> 10 10        0  1.39  4.25  3.69   2.3      10
#> ── Initial Conditions (value$inits): ──
#> depot    A1    A2 
#>     0     0     0 
#> ── First part of data (object): ──
#> # A tibble: 580 × 16
#>      id  time    ka    cl    v1    v2     q      k   k12   k21    cp ipredSim
#>   <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1     1  0.25     1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  185.     185.
#> 2     1  0.5      1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  320.     320.
#> 3     1  0.75     1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  417.     417.
#> 4     1  1        1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  486.     486.
#> 5     1  1.5      1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  564.     564.
#> 6     1  2        1  4.01  70.1  40.0  9.97 0.0573 0.142 0.249  591.     591.
#> # ℹ 574 more rows
#> # ℹ 4 more variables: sim <dbl>, depot <dbl>, A1 <dbl>, A2 <dbl>
# }
```
