# Control settings for fallback rules in parameter estimation

Control settings for fallback rules in parameter estimation

## Usage

``` r
fallback_control(
  enable_ka_fallback = TRUE,
  sigma_method_additive = "model",
  sigma_method_proportional = "model",
  sigma_fallback_fraction = 0.2
)
```

## Arguments

- enable_ka_fallback:

  Logical value indicating whether to apply a fallback to ka = 1 if the
  estimated value is invalid.

- sigma_method_additive:

  Method for additive sigma. Options are "model" or "fixed_fraction".

- sigma_method_proportional:

  Method for proportional sigma. Options are "model" or
  "fixed_fraction".

- sigma_fallback_fraction:

  Numeric value specifying the fallback fraction, for example, 0.2
  corresponds to 20 percent of the mean of observed concentrations.

## Value

A list of fallback control parameters.

## Examples

``` r
fallback_control()
#> $enable_ka_fallback
#> [1] TRUE
#> 
#> $sigma_method_additive
#> [1] "model"
#> 
#> $sigma_method_proportional
#> [1] "model"
#> 
#> $sigma_fallback_fraction
#> [1] 0.2
#> 
```
