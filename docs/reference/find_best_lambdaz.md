# Find the best terminal elimination rate constant (lambdaz)

Identifies the optimal terminal phase for lambdaz estimation using a
systematic log-linear regression approach with adjusted R-squared
optimization criteria.

## Usage

``` r
find_best_lambdaz(
  time,
  conc,
  route = "bolus",
  duration = NULL,
  adj_r_squared_threshold = 0.7,
  nlastpoints = 3,
  tolerance = 1e-04
)
```

## Arguments

- time:

  Numeric vector of observation time points.

- conc:

  Numeric vector of concentration measurements corresponding to time
  points.

- route:

  Administration method specification:

  - "bolus" (default) - Excludes time of maximum concentration (Tmax)
    point

  - "infusion" - Includes Tmax point in terminal phase evaluation

  - "oral" - Starts regression from Tmax point

- duration:

  Numeric (optional). Duration of infusion administration, in the same
  time units as `time`. Required only when `route = "infusion"`. Used to
  determine the first post-infusion observation time for terminal phase
  regression.

- adj_r_squared_threshold:

  Minimum acceptable adjusted R-squared value for valid estimation
  (default = 0.7). Values below this threshold will generate warnings.

- nlastpoints:

  Integer. Minimum number of terminal points (from the end of the
  profile) to include when evaluating candidate regression segments for
  \\\lambda_z\\. Default is `3`.

- tolerance:

  Threshold for considering adjusted R-squared values statistically
  equivalent (default = 1e-4). Used when selecting between fits with
  similar goodness-of-fit.

## Value

A list containing:

- `lambdaz`: Estimated terminal elimination rate constant
  (\\\lambda_z\\), or NA if no valid fit

- `UsedPoints`: Number of data points used in the optimal fit

- `adj.r.squared`: Adjusted R-squared value (\\R^2\\) of the optimal
  regression

- `message`: Character vector containing diagnostic messages or warnings

## Details

The algorithm implements the following decision logic:

1.  Identifies the time of maximum observed concentration (Tmax)

2.  Defines candidate terminal phases starting from the last 3
    measurable concentrations

3.  Iteratively evaluates longer time spans by including preceding data
    points

4.  For each candidate phase:

    - Performs log-concentration vs. time linear regression

    - Requires negative regression slope (positive \\\lambda_z\\)

    - Calculates adjusted R-squared metric

5.  Selects the optimal phase based on:

    - Highest adjusted R-squared value

    - When R-squared differences are \< tolerance, selects the fit with
      more points

6.  Validates final selection against R-squared threshold

## Author

Zhonghui Huang

## Examples

``` r
# Basic usage
time <- c(0.5, 1, 2, 4, 6, 8, 10)
conc <- c(12, 8, 5, 3, 2, 1.5, 1)
find_best_lambdaz(time, conc)
#> $lambdaz
#> [1] 0.1791759
#> 
#> $UsedPoints
#> [1] 4
#> 
#> $adj.r.squared
#> [1] 0.9935461
#> 
#> $message
#> [1] "Selected 4 points (higher Rsquare) Rsquare=0.9935 lambdaz=0.1792"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>       1.8035       -0.1792  
#> 
#> 

# With infusion route specification
find_best_lambdaz(time, conc, route = "bolus",duration=1)
#> $lambdaz
#> [1] 0.1791759
#> 
#> $UsedPoints
#> [1] 4
#> 
#> $adj.r.squared
#> [1] 0.9935461
#> 
#> $message
#> [1] "Selected 4 points (higher Rsquare) Rsquare=0.9935 lambdaz=0.1792"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>       1.8035       -0.1792  
#> 
#> 

# Custom threshold settings
find_best_lambdaz(time, conc, adj_r_squared_threshold = 0.8, tolerance = 0.001)
#> $lambdaz
#> [1] 0.1791759
#> 
#> $UsedPoints
#> [1] 4
#> 
#> $adj.r.squared
#> [1] 0.9935461
#> 
#> $message
#> [1] "Selected 4 points (higher Rsquare) Rsquare=0.9935 lambdaz=0.1792"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>       1.8035       -0.1792  
#> 
#> 
```
