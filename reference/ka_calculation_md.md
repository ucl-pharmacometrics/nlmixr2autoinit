# Calculate absorption rate constant (ka) in a multiple-dose one-compartment model

This estimates the absorption rate constant in a multiple-dose oral
model using first-order pharmacokinetics.

## Usage

``` r
ka_calculation_md(cl, ke, t, Ct, Fbio = 1, Dose, tau)
```

## Arguments

- cl:

  Numeric. Clearance of the drug (in L/hr).

- ke:

  Numeric. Elimination rate constant (in 1/hr).

- t:

  Numeric. Time after the last dose (in hours) at which the
  concentration is measured.

- Ct:

  Numeric. Observed concentration of the drug at time `t` (in mg/L).

- Fbio:

  Numeric. Bioavailability fraction (default = 1, meaning 100%
  bioavailability).

- Dose:

  Numeric. Administered dose of the drug (in mg).

- tau:

  Numeric. Dosing interval (in hours) between successive doses.

## Value

A list containing the following components:

- ka:

  The calculated absorption rate constant.

- full_solution:

  The full solution object returned by the
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) function, which
  includes additional details about the root-finding process.

## Details

The value of ka is obtained numerically using the uniroot unction by
solving the following equation:

\$\$Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left(
\frac{e^{-ke \cdot t}}{1 - e^{-ke \cdot \tau}} - \frac{e^{-ka \cdot
t}}{1 - e^{-ka \cdot \tau}} \right)\$\$

ka is estimated using
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html), which solves for
the root of the residual function (predicted Ct - observed Ct) within a
bounded interval (ka \> ke and ka \<= 1000)

## Author

Zhonghui Huang

## Examples

``` r
# Example from Oral_1CPT dataset (ID = 1, 5th dose, t = 2 h)
ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 852, Dose = 60000, tau = 24)
#> $ka
#> [1] 0.6135834
#> 
#> $full_solution
#> $full_solution$root
#> [1] 0.6135834
#> 
#> $full_solution$f.root
#> [1] 0.001792214
#> 
#> $full_solution$iter
#> [1] 15
#> 
#> $full_solution$init.it
#> [1] NA
#> 
#> $full_solution$estim.prec
#> [1] 6.103516e-05
#> 
#> 
#> $message
#> [1] "complete"
#> 
```
