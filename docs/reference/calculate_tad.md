# Calculate time after dose for pharmacokinetic data

Calculate time after dose (TAD) for pharmacokinetic observations.

## Usage

``` r
calculate_tad(dat, verbose = FALSE)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- verbose:

  Logical; if TRUE, prints informational messages during processing
  (e.g., when generating dose numbers). Default is FALSE.

## Value

A modified data frame with added columns:

- tad: time after dose, calculated as the observation time minus the
  time of the most recent prior dose; set to NA for dosing records

- iiobs: interdose interval inherited from the most recent dosing record

- rateobs: infusion rate inherited from the most recent dosing record

- routeobs (optional): route of administration inherited from the most
  recent dosing record, included only if route information is present

- dose_number: sequential dose number, generated if not already present

## Details

The procedure identifies dosing events based on the event identifier
(EVID) and assigns each observation the attributes of the most recent
prior dose. The time after dose is then calculated for observation rows.
If `dose_number` column is not present in the input, it is automatically
created for each subject.

## Author

Zhonghui Huang

## Examples

``` r
calculate_tad(Bolus_1CPT)
#> # A tibble: 7,920 × 21
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    SS    II
#>    <int> <dbl> <dbl> <dbl> <int> <int> <int>    <int> <dbl> <dbl> <int> <int>
#>  1     1  0       0   0        1 60000     1    60000  65.1  4.07    99     0
#>  2     1  0.25 1126.  7.03     0     0     0    60000  65.1  4.07    99     0
#>  3     1  0.5   870.  6.77     0     0     0    60000  65.1  4.07    99     0
#>  4     1  0.75  884.  6.78     0     0     0    60000  65.1  4.07    99     0
#>  5     1  1    1244   7.13     0     0     0    60000  65.1  4.07    99     0
#>  6     1  1.5   995.  6.90     0     0     0    60000  65.1  4.07    99     0
#>  7     1  2     946.  6.85     0     0     0    60000  65.1  4.07    99     0
#>  8     1  2.5   589.  6.38     0     0     0    60000  65.1  4.07    99     0
#>  9     1  3     754.  6.63     0     0     0    60000  65.1  4.07    99     0
#> 10     1  4    1061.  6.97     0     0     0    60000  65.1  4.07    99     0
#> # ℹ 7,910 more rows
#> # ℹ 9 more variables: SD <int>, CMT <int>, resetflag <dbl>, dose_number <int>,
#> #   RATE <dbl>, tad <dbl>, dose <int>, iiobs <dbl>, rateobs <dbl>
calculate_tad(Infusion_1CPT)
#> # A tibble: 7,920 × 21
#>       ID  TIME    DV  LNDV   MDV    AMT   RATE  EVID raw_dose     V    CL    SS
#>    <int> <dbl> <dbl> <dbl> <int>  <int>  <int> <int>    <int> <dbl> <dbl> <int>
#>  1     1  0       0   0        1 120000 120000     1   120000  75.2  3.36    99
#>  2     1  0.25  343.  5.84     0      0      0     0   120000  75.2  3.36    99
#>  3     1  0.5   756.  6.63     0      0      0     0   120000  75.2  3.36    99
#>  4     1  0.75  888.  6.79     0      0      0     0   120000  75.2  3.36    99
#>  5     1  1    1620.  7.39     0      0      0     0   120000  75.2  3.36    99
#>  6     1  1.5  1511.  7.32     0      0      0     0   120000  75.2  3.36    99
#>  7     1  2    1690.  7.43     0      0      0     0   120000  75.2  3.36    99
#>  8     1  2.5  1069.  6.97     0      0      0     0   120000  75.2  3.36    99
#>  9     1  3    1114.  7.02     0      0      0     0   120000  75.2  3.36    99
#> 10     1  4    1164.  7.06     0      0      0     0   120000  75.2  3.36    99
#> # ℹ 7,910 more rows
#> # ℹ 9 more variables: II <int>, SD <int>, CMT <int>, resetflag <dbl>,
#> #   dose_number <int>, tad <dbl>, dose <int>, iiobs <dbl>, rateobs <dbl>
calculate_tad(Oral_1CPT)
#> # A tibble: 7,920 × 22
#>       ID  TIME    DV  LNDV   MDV   AMT  EVID raw_dose     V    CL    KA    SS
#>    <int> <dbl> <dbl> <dbl> <int> <int> <int>    <int> <dbl> <dbl> <dbl> <int>
#>  1     1  0       0   0        1 60000     1    60000  86.5  4.88  1.09    99
#>  2     1  0.25  205.  5.32     0     0     0    60000  86.5  4.88  1.09    99
#>  3     1  0.5   311.  5.74     0     0     0    60000  86.5  4.88  1.09    99
#>  4     1  0.75  389.  5.96     0     0     0    60000  86.5  4.88  1.09    99
#>  5     1  1     620.  6.43     0     0     0    60000  86.5  4.88  1.09    99
#>  6     1  1.5   400   5.99     0     0     0    60000  86.5  4.88  1.09    99
#>  7     1  2     538.  6.29     0     0     0    60000  86.5  4.88  1.09    99
#>  8     1  2.5   975.  6.88     0     0     0    60000  86.5  4.88  1.09    99
#>  9     1  3     457.  6.13     0     0     0    60000  86.5  4.88  1.09    99
#> 10     1  4     580.  6.36     0     0     0    60000  86.5  4.88  1.09    99
#> # ℹ 7,910 more rows
#> # ℹ 10 more variables: II <int>, SD <int>, CMT <int>, resetflag <dbl>,
#> #   dose_number <int>, RATE <dbl>, tad <dbl>, dose <int>, iiobs <dbl>,
#> #   rateobs <dbl>
```
