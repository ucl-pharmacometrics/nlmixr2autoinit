# Calculate the absorption rate constant using the Wagner-Nelson method

Calculates absorption rate constant using the Wagner–Nelson method for
single-dose extravascular pharmacokinetics.

## Usage

``` r
ka_wanger_nelson(dat, nca.out = NULL)
```

## Arguments

- dat:

  A data frame containing two columns: 'TIME' for sampling time points
  and 'DV' for observed plasma drug concentrations.

- nca.out:

  Optional object containing results from a previous noncompartmental
  analysis. It must include 'auc0_inf' for the area under the
  concentration-time curve extrapolated to infinity and 'lambdaz' for
  the terminal elimination rate constant. If not provided, the function
  calls 'getnca' internally using the input data.

## Value

A list containing:

- ka: Estimated absorption rate constant

- dat_out_wanger_nelson: Input data frame augmented with calculated
  pharmacokinetic variables including cumulative AUC, fraction absorbed,
  and fraction remaining

## Details

The Wagner-Nelson method estimates the fraction of drug absorbed over
time based on the principle of mass balance, where the unabsorbed
fraction is quantified as the proportion of the administered dose that
has not yet entered systemic circulation. A linear regression is applied
to the natural logarithm of the unabsorbed fraction versus time, and the
negative slope of this regression corresponds to the first-order
absorption rate constant 'ka'.

Key assumptions:

- Single-dose oral or extravascular administration

- First-order absorption and first-order elimination

- Linear pharmacokinetics with 'ka' greater than 'ke'

Computational steps:

- AUC is calculated using trapezoidal integration.

- The fraction absorbed is calculated from AUC and the terminal
  elimination rate constant.

- The remaining fraction is transformed using the natural logarithm.

- Linear regression of log(remaining fraction) against time yields 'ka'.

## References

Wagner JG and Nelson E (1963). Percent absorbed time plots derived from
blood level and/or urinary excretion data. Journal of Pharmaceutical
Sciences, 52(6), 610-611.

## Author

Zhonghui Huang

## Examples

``` r
# Simulated one-compartment oral absorption data
Dose <- 100
Fbio <- 1
Vd   <- 70
CL   <- 4
ka   <- 1.2
ke   <- CL / Vd
t  <- seq(0.5, 8, by = 0.5)
Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))

dat <- data.frame(TIME = t, DV = Ct)

ka_wanger_nelson(dat)
#> $ka
#> [1] 1.220523
#> 
#> $dat_out_wanger_nelson
#>    TIME        DV auc_intervals auc_accumulate  frac_abs frac_remained
#> 1   0.5 0.6345319     0.1586330      0.1586330 0.4536953  0.5463047204
#> 2   1.0 0.9648974     0.3998573      0.5584903 0.7024742  0.2975257803
#> 3   1.5 1.1288363     0.5234334      1.0819237 0.8387987  0.1612012752
#> 4   2.0 1.2019277     0.5826910      1.6646147 0.9134128  0.0865872063
#> 5   2.5 1.2256362     0.6068910      2.2715057 0.9541652  0.0458348416
#> 6   3.0 1.2227051     0.6120853      2.8835910 0.9763394  0.0236606336
#> 7   3.5 1.2056028     0.6070770      3.4906680 0.9883231  0.0116769496
#> 8   4.0 1.1811596     0.5966906      4.0873586 0.9947193  0.0052807108
#> 9   4.5 1.1531117     0.5835678      4.6709264 0.9980542  0.0019458416
#> 10  5.0 1.1234978     0.5691524      5.2400788 0.9997139  0.0002861454
#> 11  5.5 1.0934250     0.5542307      5.7943095 1.0004590 -0.0004589994
#> 12  6.0 1.0634895     0.5392286      6.3335381 1.0007069 -0.0007068955
#> 13  6.5 1.0340078     0.5243743      6.8579124 1.0006864 -0.0006864320
#> 14  7.0 1.0051428     0.5097876      7.3677001 1.0005231 -0.0005230981
#> 15  7.5 0.9769735     0.4955291      7.8632291 1.0002856 -0.0002856395
#> 16  8.0 0.9495332     0.4816267      8.3448558 1.0000117 -0.0000116642
#> 
```
