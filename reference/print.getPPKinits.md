# Print method for `getPPKinits` objects

Prints a summary of the results from the initial parameter estimation
pipeline, including recommended initial estimates, ETA variance
estimates, and parameter descriptions. It is the default S3 `print`
method for objects of class `getPPKinits`.

## Usage

``` r
# S3 method for class 'getPPKinits'
print(x, ...)
```

## Arguments

- x:

  An object of class `getPPKinits` containing the initial parameter
  estimation results. Expected components include:

  - `Recommended_initial_estimates`: A data frame with estimated values
    and selection methods.

  - `Parameter.descriptions`: A character vector explaining the meaning
    of each parameter.

  - `time.spent`: Time taken to compute the estimates.

- ...:

  Additional arguments (for compatibility with the generic
  [`print()`](https://rdrr.io/r/base/print.html)).

## Value

Prints a formatted summary to the console.

## Examples

``` r
# \donttest{
## Oral example
inits.out <- getPPKinits(Bolus_1CPT)
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
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 11.26 h
#> Evaluating the predictive performance of calculated one-compartment model parameters....................
#> Error in loadNamespace(x): there is no package called ‘progress’
print(inits.out)
#> Error: object 'inits.out' not found
# }
```
