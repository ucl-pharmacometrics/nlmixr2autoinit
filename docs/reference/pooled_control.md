# Control settings for pooled data analysis

Defines control parameters for time binning and preprocessing in pooled
data analysis. These parameters are typically passed to
'get_pooled_data'.

## Usage

``` r
pooled_control(
  nbins = 10,
  bin_method = c("quantile", "jenks", "kmeans", "pretty", "sd", "equal", "density"),
  tad_rounding = TRUE
)
```

## Arguments

- nbins:

  Integer or the character string auto. Number of time bins used to
  group observations. Default is 10.

- bin_method:

  Character string specifying the binning method. Must be one of
  "quantile", "jenks", "kmeans", "pretty", "sd", "equal", or "density".

- tad_rounding:

  Logical value indicating whether tad and the most common dosing
  interval should be rounded to the nearest whole unit before
  comparison. Default is TRUE, allowing small deviations (for example, a
  tad of 24.3 is treated as within a 24-unit interval).

## Value

A named list containing control parameters for pooled pharmacokinetic
analysis.

## See also

[get_pooled_data](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)

## Examples

``` r
pooled_control()
#> $nbins
#> [1] 10
#> 
#> $bin_method
#> [1] "quantile"
#> 
#> $tad_rounding
#> [1] TRUE
#> 
```
