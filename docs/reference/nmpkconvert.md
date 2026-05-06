# Expand additional dosing (ADDL) records for pharmacokinetic analysis

Expands dosing records that contain additional doses (ADDL) using the
specified interdose interval (II). Each additional dose is converted
into an explicit record to provide a complete dosing history suitable
for population pharmacokinetic modeling.

## Usage

``` r
nmpkconvert(dat)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

## Value

A data frame with expanded dosing records. The columns ADDL and II are
reset to zero after expansion.

## Details

Dosing records with ADDL greater than zero are expanded using the
formula: TIME_new = TIME + n × II, where n ranges from 1 up to ADDL.
Observation records (EVID = 0) are not modified.

## Examples

``` r
# Dataset with a single subject and additional dosing
dat <- data.frame(
  ID   = rep(1, 6),
  EVID = c(1, 0, 0, 1, 0, 0),
  ADDL = c(2, 0, 0, 0, 0, 0),
  TIME = c(0, 1, 2, 3, 4, 5),
  II   = c(24, 0, 0, 0, 0, 0),
  AMT  = c(100, 0, 0, 0, 0, 0)
)
nmpkconvert(dat)
#>     ID EVID ADDL TIME II AMT
#> 1    1    1    0    0  0 100
#> 2    1    0    0    1  0   0
#> 3    1    0    0    2  0   0
#> 21   1    1    0    3  0   0
#> 5    1    0    0    4  0   0
#> 6    1    0    0    5  0   0
#> 11   1    1    0   24  0 100
#> 1.1  1    1    0   48  0 100
```
