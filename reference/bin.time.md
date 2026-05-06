# Bin time-concentration data using quantile or algorithmic binning

Bins data by time using either equal-frequency (quantile) binning or
algorithmic binning methods.

## Usage

``` r
bin.time(
  dat,
  nbins = "auto",
  bin.method = c("quantile", "jenks", "kmeans", "pretty", "sd", "equal", "density")
)
```

## Arguments

- dat:

  A data frame containing PK data. Must include:

  - tad: time after dose

  - DVstd: standardized concentration (DV/dose)

  - EVID: optional event ID column used to filter observations (EVID ==
    0)

- nbins:

  Number of bins or "auto". If numeric with `bin.method = "quantile"`,
  specifies equal-frequency bins. If "auto", 10 bins are used for
  quantile; otherwise binning is determined by
  [`vpc::auto_bin`](https://rdrr.io/pkg/vpc/man/auto_bin.html). Numeric
  nbins for non-quantile methods is passed to
  [`vpc::auto_bin`](https://rdrr.io/pkg/vpc/man/auto_bin.html).

- bin.method:

  Binning strategy (default = "quantile"). Available options are:

  - quantile: equal-frequency binning by empirical quantiles

  - jenks: natural breaks minimizing within-bin variance

  - kmeans, pretty, sd, equal, density: alternative binning methods from
    vpc::auto_bin

## Value

A list containing summary results of the time-concentration binning
process.

## Details

Supports quantile-based binning and other data-driven methods (jenks,
kmeans, pretty, sd, equal, density), with optional automatic bin count
selection.

## See also

[`vpc::auto_bin`](https://rdrr.io/pkg/vpc/man/auto_bin.html)

## Author

Zhonghui Huang

## Examples

``` r
dat <- Bolus_1CPT
dat <- nmpkconvert(dat)
dat <- calculate_tad(dat)
dat$DVstd <- dat$DV / dat$dose
bin.time(dat)
#> $binned.df
#>     Time         Conc
#> 1   0.50 0.0180495833
#> 2   1.25 0.0167466667
#> 3   2.25 0.0157725000
#> 4   3.50 0.0143850000
#> 5   7.00 0.0116075000
#> 6  14.00 0.0076733333
#> 7  20.00 0.0047000000
#> 8  23.99 0.0044500000
#> 9  42.00 0.0014350000
#> 10 60.00 0.0003566667
#> 
#> $bin_limits.df
#>    Group  Lower  Upper
#> 1      1  0.250  0.775
#> 2      2  0.775  1.600
#> 3      3  1.600  2.650
#> 4      4  2.650  4.800
#> 5      5  4.800 10.000
#> 6      6 10.000 18.400
#> 7      7 18.400 23.990
#> 8      8 23.990 33.600
#> 9      9 33.600 58.800
#> 10    10 58.800 72.000
#> 
#> $breaks
#>  [1]  0.250  0.775  1.600  2.650  4.800 10.000 18.400 23.990 33.600 58.800
#> [11] 72.000
#> 
#> $method_used
#> [1] "quantile"
#> 
#> $nbins_final
#> [1] 10
#> 
```
