# Computes the trimmed geometric mean

Computes the trimmed geometric mean of a numeric vector

## Usage

``` r
trimmed_geom_mean(x, trim = 0, na.rm = TRUE)
```

## Arguments

- x:

  A numeric vector containing the values for geometric mean calculation.

- trim:

  A numeric value between 0 and 0.5 indicating the proportion of values
  to be trimmed from each end of the vector. Default is 0.

- na.rm:

  Logical value indicating whether missing values should be removed
  before computation. Default is TRUE.

## Value

A numeric value representing the trimmed geometric mean. Returns NA if
no values remain after trimming.

## Examples

``` r
x <- c(1, 2, 3, 4, 5, 100)
trimmed_geom_mean(x, trim = 0.05)
#> [1] 4.784797
```
