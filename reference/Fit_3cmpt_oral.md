# Fit oral pharmacokinetic data to a three-compartment linear elimination model

Fits oral pharmacokinetic data to a three-compartment model with
first-order absorption and first-order elimination using the naive
pooled data approach. Supports multiple estimation methods provided by
nlmixr2 and can optionally return only predicted concentrations to
support efficient simulation workflows.

## Usage

``` r
Fit_3cmpt_oral(
  data,
  est.method,
  input.ka,
  input.cl,
  input.vc3cmpt,
  input.vp3cmpt,
  input.vp23cmpt,
  input.q3cmpt,
  input.q23cmpt,
  input.add,
  ncores = 2,
  return.pred.only = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame containing oral pharmacokinetic data formatted for
  nlmixr2.

- est.method:

  Estimation method to use in nlmixr2. Must be one of: `"rxSolve"`,
  `"nls"`, `"nlm"`, `"nlminb"`, or `"focei"`.

- input.ka:

  Initial estimate of the absorption rate constant (ka).

- input.cl:

  Initial estimate of clearance (CL).

- input.vc3cmpt:

  Initial estimate of central volume of distribution (V1).

- input.vp3cmpt:

  Initial estimate of first peripheral volume of distribution (V2).

- input.vp23cmpt:

  Initial estimate of second peripheral volume of distribution (V3).

- input.q3cmpt:

  Initial estimate of first inter-compartmental clearance (Q1).

- input.q23cmpt:

  Initial estimate of second inter-compartmental clearance (Q2).

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
Fit_3cmpt_oral(
  data = dat,
  est.method = "rxSolve",
  input.ka = 1,
  input.cl = 4,
  input.vc3cmpt = 70,
  input.vp3cmpt = 35,
  input.vp23cmpt = 35,
  input.q3cmpt = 4,
  input.q23cmpt = 4,
  input.add = 10
)
#> ── Solved rxode2 object ──
#> ── Parameters (value$params): ──
#> # A tibble: 10 × 9
#>    id      tka   tcl   tv1   tv2   tv3    tq   tq2 add.err
#>    <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1 1         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  2 2         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  3 3         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  4 4         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  5 5         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  6 6         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  7 7         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  8 8         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#>  9 9         0  1.39  4.25  3.56  3.56  1.39  1.39      10
#> 10 10        0  1.39  4.25  3.56  3.56  1.39  1.39      10
#> ── Initial Conditions (value$inits): ──
#> depot    A1    A2    A3 
#>     0     0     0     0 
#> ── First part of data (object): ──
#> # A tibble: 580 × 21
#>      id  time    ka    cl    v1    v2    v3    q1    q2      k    k12   k21
#>   <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>
#> 1     1  0.25     1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> 2     1  0.5      1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> 3     1  0.75     1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> 4     1  1        1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> 5     1  1.5      1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> 6     1  2        1  4.01  70.1  35.2  35.2  4.01  4.01 0.0573 0.0573 0.114
#> # ℹ 574 more rows
#> # ℹ 9 more variables: k13 <dbl>, k31 <dbl>, cp <dbl>, ipredSim <dbl>,
#> #   sim <dbl>, depot <dbl>, A1 <dbl>, A2 <dbl>, A3 <dbl>
# }
```
