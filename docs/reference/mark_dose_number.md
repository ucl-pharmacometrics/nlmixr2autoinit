# Mark dose number

Assigns sequential dose numbers based on dosing events (EVID) within
each subject.

## Usage

``` r
mark_dose_number(dat)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

## Value

A modified data frame with an added column named dose_number, indicating
the sequential dose count within each subject and reset group.

## Author

Zhonghui Huang

## Examples

``` r
mark_dose_number(Bolus_1CPT)
#> # A tibble: 7,920 × 16
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID  DOSE     V    CL    SS    II    SD
#>    <int> <dbl> <dbl> <dbl> <int> <int> <int> <int> <dbl> <dbl> <int> <int> <int>
#>  1     1  0       0   0        1 60000     1 60000  65.1  4.07    99     0     1
#>  2     1  0.25 1126.  7.03     0     0     0 60000  65.1  4.07    99     0     1
#>  3     1  0.5   870.  6.77     0     0     0 60000  65.1  4.07    99     0     1
#>  4     1  0.75  884.  6.78     0     0     0 60000  65.1  4.07    99     0     1
#>  5     1  1    1244   7.13     0     0     0 60000  65.1  4.07    99     0     1
#>  6     1  1.5   995.  6.90     0     0     0 60000  65.1  4.07    99     0     1
#>  7     1  2     946.  6.85     0     0     0 60000  65.1  4.07    99     0     1
#>  8     1  2.5   589.  6.38     0     0     0 60000  65.1  4.07    99     0     1
#>  9     1  3     754.  6.63     0     0     0 60000  65.1  4.07    99     0     1
#> 10     1  4    1061.  6.97     0     0     0 60000  65.1  4.07    99     0     1
#> # ℹ 7,910 more rows
#> # ℹ 3 more variables: CMT <int>, resetflag <dbl>, dose_number <int>
mark_dose_number(Infusion_1CPT)
#> # A tibble: 7,920 × 17
#>       ID  TIME    DV  LNDV   MDV    AMT   RATE  EVID   DOSE     V    CL    SS
#>    <int> <dbl> <dbl> <dbl> <int>  <int>  <int> <int>  <int> <dbl> <dbl> <int>
#>  1     1  0       0   0        1 120000 120000     1 120000  75.2  3.36    99
#>  2     1  0.25  343.  5.84     0      0      0     0 120000  75.2  3.36    99
#>  3     1  0.5   756.  6.63     0      0      0     0 120000  75.2  3.36    99
#>  4     1  0.75  888.  6.79     0      0      0     0 120000  75.2  3.36    99
#>  5     1  1    1620.  7.39     0      0      0     0 120000  75.2  3.36    99
#>  6     1  1.5  1511.  7.32     0      0      0     0 120000  75.2  3.36    99
#>  7     1  2    1690.  7.43     0      0      0     0 120000  75.2  3.36    99
#>  8     1  2.5  1069.  6.97     0      0      0     0 120000  75.2  3.36    99
#>  9     1  3    1114.  7.02     0      0      0     0 120000  75.2  3.36    99
#> 10     1  4    1164.  7.06     0      0      0     0 120000  75.2  3.36    99
#> # ℹ 7,910 more rows
#> # ℹ 5 more variables: II <int>, SD <int>, CMT <int>, resetflag <dbl>,
#> #   dose_number <int>
mark_dose_number(Oral_1CPT)
#> # A tibble: 7,920 × 17
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID  DOSE     V    CL    KA    SS    II
#>    <int> <dbl> <dbl> <dbl> <int> <int> <int> <int> <dbl> <dbl> <dbl> <int> <int>
#>  1     1  0       0   0        1 60000     1 60000  86.5  4.88  1.09    99     0
#>  2     1  0.25  205.  5.32     0     0     0 60000  86.5  4.88  1.09    99     0
#>  3     1  0.5   311.  5.74     0     0     0 60000  86.5  4.88  1.09    99     0
#>  4     1  0.75  389.  5.96     0     0     0 60000  86.5  4.88  1.09    99     0
#>  5     1  1     620.  6.43     0     0     0 60000  86.5  4.88  1.09    99     0
#>  6     1  1.5   400   5.99     0     0     0 60000  86.5  4.88  1.09    99     0
#>  7     1  2     538.  6.29     0     0     0 60000  86.5  4.88  1.09    99     0
#>  8     1  2.5   975.  6.88     0     0     0 60000  86.5  4.88  1.09    99     0
#>  9     1  3     457.  6.13     0     0     0 60000  86.5  4.88  1.09    99     0
#> 10     1  4     580.  6.36     0     0     0 60000  86.5  4.88  1.09    99     0
#> # ℹ 7,910 more rows
#> # ℹ 4 more variables: SD <int>, CMT <int>, resetflag <dbl>, dose_number <int>
```
