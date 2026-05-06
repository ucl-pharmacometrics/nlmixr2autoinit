# Compute overall residual variability from elimination phase

Applies `getsigmas` to each individual and dose group after filtering
observation records (EVID == 0), and calculates trimmed mean estimates
of additive and proportional residual variability.

## Usage

``` r
getsigma(df, nlastpoints = 3, sigma_trim = 0.05)
```

## Arguments

- df:

  Full pharmacokinetic dataset containing at least the columns: EVID,
  ID, TIME, DV, and routeobs.

- nlastpoints:

  Number of terminal points used for elimination phase regression in
  each group (passed to getsigmas).

- sigma_trim:

  Trimming proportion used when calculating trimmed means of residual
  standard deviations. Default is 0.05.

## Value

A list containing:

- summary: Named list with trimmed mean values of additive and
  proportional residual variability

- full: Data frame with residual estimates for each individual-dose
  group

## Details

The function groups the dataset by subject and dose occasion, applies
elimination-phase residual analysis using `getsigmas`, and summarizes
the individual residual standard deviations by their trimmed means. This
provides population-level estimates of additive and proportional
residual unexplained variability (RUV).

## See also

[getsigmas](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getsigmas.md)

## Author

Zhonghui Huang

## Examples

``` r
dat <- Bolus_1CPT
dat <- processData(dat)$dat
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
getsigma(dat)
#> $summary
#> $summary$sigma_additive
#> [1] 11.11136
#> 
#> $summary$sigma_proportional
#> [1] 0.1138021
#> 
#> 
#> $full
#> # A tibble: 960 × 7
#>       ID resetflag dose_number intercept   slope residual_sd_additive
#>    <int>     <int>       <int>     <dbl>   <dbl>                <dbl>
#>  1     1         1           1      6.97 -0.0679                1.70 
#>  2     1         1           2     NA    NA                    NA    
#>  3     1         1           3     NA    NA                    NA    
#>  4     1         1           4     NA    NA                    NA    
#>  5     1         1           5     19.6  -0.0843               72.9  
#>  6     1         1           6     NA    NA                    NA    
#>  7     1         1           7     NA    NA                    NA    
#>  8     1         2           1     15.5  -0.0446                0.608
#>  9     2         1           1      6.53 -0.0732                0.892
#> 10     2         1           2     NA    NA                    NA    
#> # ℹ 950 more rows
#> # ℹ 1 more variable: residual_sd_proportional <dbl>
#> 
```
