# Estimate absorption rate constant in a one-compartment oral model

This estimates the absorption rate constant in a single-dose oral model
using first-order pharmacokinetics.

## Usage

``` r
ka_calculation_sd(cl, ke, t, Ct, Fbio = 1, Dose)
```

## Arguments

- cl:

  Numeric. Clearance of the drug.

- ke:

  Numeric. Elimination rate constant.

- t:

  Numeric. Time after administration.

- Ct:

  Numeric. Observed plasma concentration at time t.

- Fbio:

  Numeric. Absolute bioavailability fraction. Default is 1.

- Dose:

  Numeric. Administered oral dose.

## Value

A list containing:

- ka:

  Estimated absorption rate constant.

- full_solution:

  The full result object returned by the root-finding process.

- message:

  A character string indicating the status of the estimation or any
  warnings.

## Details

The model assumes a one-compartment structure with first-order
absorption and first-order elimination.

The concentration-time relationship is: \$\$Ct = \frac{Fbio \cdot Dose
\cdot ka}{Vd \cdot (ka - ke)} \left( e^{-ke \cdot t} - e^{-ka \cdot t}
\right)\$\$ where the volume of distribution is defined as: \$\$Vd =
\frac{cl}{ke}\$\$

ka is estimated using
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html), which solves for
the root of the residual function (predicted Ct - observed Ct) within a
bounded interval (ka \> ke and ka \<= 1000)

## Author

Zhonghui Huang

## Examples

``` r
# Example from Oral_1CPT dataset (ID = 1, 1st dose, t = 0.5 h)
ka_calculation_sd(cl = 3.62, ke = 0.0556, t = 0.5, Ct = 310.6, Dose = 60000)
#> $ka
#> [1] 0.837342
#> 
#> $full_solution
#> $full_solution$root
#> [1] 0.837342
#> 
#> $full_solution$f.root
#> [1] 7.349658e-05
#> 
#> $full_solution$iter
#> [1] 11
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
