# Generate pooled data for pharmacokinetic analysis

Processes pharmacokinetic data and produces pooled datasets according to
the dosing context. Data can be grouped based on first dose, repeated
dosing, or a combination of both, with control over binning and time
alignment.

## Usage

``` r
get_pooled_data(
  dat,
  dose_type = c("first_dose", "repeated_doses", "combined_doses"),
  pooled_ctrl = pooled_control()
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- dose_type:

  Specifies the dosing context of the pharmacokinetic observations.
  Classified as:

  - first_dose: data include only observations following the initial
    administration

  - repeated_doses: data include only observations during repeated or
    steady-state dosing

  - combined_doses: data include observations from both first-dose and
    repeated-dose intervals

- pooled_ctrl:

  A list of control parameters created by 'pooled_control', including
  settings for binning and time rounding.

## Value

A list containing pooled pharmacokinetic datasets depending on the
specified dose type:

- datpooled_fd: pooled data for first-dose observations

- datpooled_efd: pooled data for repeated dosing

- datpooled_all: pooled data combining first-dose and repeated-dose
  observations

## Details

For repeated-doses and combined-doses classifications, the most common
interdose interval is identified from dosing records and used to
determine whether observations fall within the relevant interval. If
tad_rounding is TRUE, both time after dose and dosing interval are
rounded before comparison.

## See also

[pooled_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md),
[trimmed_geom_mean](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trimmed_geom_mean.md)

## Author

Zhonghui Huang

## Examples

``` r
dat <- processData(Bolus_1CPT)$dat
#> 
#> 
#> Infometrics                               Value          
#> ----------------------------------------  ---------------
#> Dose Route                                bolus          
#> Dose Type                                 combined_doses 
#> Number of Subjects                        120            
#> Number of Observations                    6951           
#> Subjects with First-Dose Interval Data    120            
#> Observations in the First-Dose Interval   2276           
#> Subjects with Multiple-Dose Data          120            
#> Observations after Multiple Doses         4675           
#> ----------------------------------------  ------
get_pooled_data(dat, dose_type = "combined_doses")
#> $datpooled_fd
#> $datpooled_fd$binned.df
#>      Time        Conc
#> 1   0.375 0.014535000
#> 2   0.875 0.013858333
#> 3   1.750 0.013546667
#> 4   2.750 0.012431667
#> 5   5.000 0.010956667
#> 6   8.000 0.008958333
#> 7  14.000 0.006098333
#> 8  22.000 0.003826667
#> 9  42.000 0.001225000
#> 10 60.000 0.000310000
#> 
#> $datpooled_fd$bin_limits.df
#>    Group Lower Upper
#> 1      1  0.25  0.70
#> 2      2  0.70  1.30
#> 3      3  1.30  2.20
#> 4      4  2.20  3.20
#> 5      5  3.20  6.00
#> 6      6  6.00 11.20
#> 7      7 11.20 18.40
#> 8      8 18.40 28.80
#> 9      9 28.80 50.40
#> 10    10 50.40 71.99
#> 
#> $datpooled_fd$breaks
#>  [1]  0.25  0.70  1.30  2.20  3.20  6.00 11.20 18.40 28.80 50.40 71.99
#> 
#> $datpooled_fd$method_used
#> [1] "quantile"
#> 
#> $datpooled_fd$nbins_final
#> [1] 10
#> 
#> 
#> $datpooled_efd
#> $datpooled_efd$binned.df
#>      Time        Conc
#> 1   0.375 0.019555000
#> 2   0.875 0.018667083
#> 3   1.500 0.017518333
#> 4   2.250 0.017312500
#> 5   3.500 0.015998333
#> 6   6.000 0.013527500
#> 7  10.000 0.010872500
#> 8  16.000 0.007835833
#> 9  21.995 0.004849583
#> 10 23.990 0.004757500
#> 
#> $datpooled_efd$bin_limits.df
#>    Group Lower Upper
#> 1      1  0.25  0.65
#> 2      2  0.65  1.10
#> 3      3  1.10  1.90
#> 4      4  1.90  2.70
#> 5      5  2.70  4.00
#> 6      6  4.00  7.20
#> 7      7  7.20 12.80
#> 8      8 12.80 19.20
#> 9      9 19.20 23.99
#> 10    10 23.99 24.00
#> 
#> $datpooled_efd$breaks
#>  [1]  0.25  0.65  1.10  1.90  2.70  4.00  7.20 12.80 19.20 23.99 24.00
#> 
#> $datpooled_efd$method_used
#> [1] "quantile"
#> 
#> $datpooled_efd$nbins_final
#> [1] 10
#> 
#> 
#> $datpooled_all
#> $datpooled_all$binned.df
#>      Time        Conc
#> 1   0.375 0.018370000
#> 2   0.875 0.017265833
#> 3   1.500 0.016350833
#> 4   2.250 0.015772500
#> 5   3.500 0.014385000
#> 6   6.000 0.012228750
#> 7  10.000 0.009681667
#> 8  16.000 0.006780833
#> 9  20.000 0.004700000
#> 10 23.990 0.004450000
#> 
#> $datpooled_all$bin_limits.df
#>    Group Lower Upper
#> 1      1  0.25  0.65
#> 2      2  0.65  1.10
#> 3      3  1.10  1.90
#> 4      4  1.90  2.70
#> 5      5  2.70  4.00
#> 6      6  4.00  7.20
#> 7      7  7.20 12.80
#> 8      8 12.80 19.20
#> 9      9 19.20 23.99
#> 10    10 23.99 24.00
#> 
#> $datpooled_all$breaks
#>  [1]  0.25  0.65  1.10  1.90  2.70  4.00  7.20 12.80 19.20 23.99 24.00
#> 
#> $datpooled_all$method_used
#> [1] "quantile"
#> 
#> $datpooled_all$nbins_final
#> [1] 10
#> 
#> 
```
