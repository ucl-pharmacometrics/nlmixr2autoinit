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
#> Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 4, estimated Vd: 66 
#> Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
#> Run parameter sweeping on multi-compartmental PK parameters....................
print(inits.out)
#> ===============Initial Parameter Estimation Summary ===============
#> 
#> Recommended initial estimates :
#>            Parameters               Methods    Values
#> 1                  Ka                    IV        NA
#> 2                  CL Naive pooled NCA (MD)     4.000
#> 3                  Vd Naive pooled NCA (MD)    66.000
#> 4                Vmax    Parameter sweeping 21551.820
#> 5                  Km    Parameter sweeping  5321.437
#> 6           Vc(2CMPT)    Parameter sweeping    66.000
#> 7           Vp(2CMPT)    Parameter sweeping     6.600
#> 8            Q(2CMPT)    Parameter sweeping     1.000
#> 9           Vc(3CMPT)    Parameter sweeping    47.601
#> 10          Vp(3CMPT)    Parameter sweeping     9.520
#> 11         Vp2(3CMPT)    Parameter sweeping     9.520
#> 12           Q(3CMPT)    Parameter sweeping     8.000
#> 13          Q2(3CMPT)    Parameter sweeping     8.000
#> 14     Sigma additive           Model-based    11.111
#> 15 Sigma proportional           Model-based     0.114
#> 
#> Time spent :
#> [1] "39.453s"
#> 
#> ETA variances and derived covariances:
#>         Parameters                  Methods Values
#> 1           eta.ka             fixed_values    0.1
#> 2           eta.cl             fixed_values    0.1
#> 3           eta.vc             fixed_values    0.1
#> 4           eta.vp             fixed_values    0.1
#> 5            eta.q             fixed_values    0.1
#> 6          eta.vp2             fixed_values    0.1
#> 7           eta.q2             fixed_values    0.1
#> 8         eta.vmax             fixed_values    0.1
#> 9           eta.km             fixed_values    0.1
#> 10 cor.eta_vmax_km eta_corr_derived (r=0.1)   0.01
#> 11   cor.eta_cl_vc eta_corr_derived (r=0.1)   0.01
#> 12   cor.eta_cl_vp eta_corr_derived (r=0.1)   0.01
#> 13  cor.eta_cl_vp2 eta_corr_derived (r=0.1)   0.01
#> 14    cor.eta_cl_q eta_corr_derived (r=0.1)   0.01
#> 15   cor.eta_cl_q2 eta_corr_derived (r=0.1)   0.01
#> 16   cor.eta_vc_vp eta_corr_derived (r=0.1)   0.01
#> 17  cor.eta_vc_vp2 eta_corr_derived (r=0.1)   0.01
#> 18    cor.eta_vc_q eta_corr_derived (r=0.1)   0.01
#> 19   cor.eta_vc_q2 eta_corr_derived (r=0.1)   0.01
#> 20  cor.eta_vp_vp2 eta_corr_derived (r=0.1)   0.01
#> 21    cor.eta_vp_q eta_corr_derived (r=0.1)   0.01
#> 22   cor.eta_vp_q2 eta_corr_derived (r=0.1)   0.01
#> 23   cor.eta_vp2_q eta_corr_derived (r=0.1)   0.01
#> 24  cor.eta_vp2_q2 eta_corr_derived (r=0.1)   0.01
#> 25    cor.eta_q_q2 eta_corr_derived (r=0.1)   0.01
#> Note: The ETA variances and covariances listed above are predefined default initialization values automatically assigned by the package.
#> 
#> Parameter descriptions:
#>  [1] "Ka: absorption constant rate"                                                       
#>  [2] "CL: clearance"                                                                      
#>  [3] "Vd: volume of distribution"                                                         
#>  [4] "Vmax: maximum metabolic rate"                                                       
#>  [5] "Km: Michaelis constant"                                                             
#>  [6] "Vc: volume of distribution of the central compartment"                              
#>  [7] "Vp: volume of distribution of the peripheral compartment"                           
#>  [8] "Vp2: volume of distribution of the second peripheral compartment"                   
#>  [9] "Q: inter-compartmental clearance"                                                   
#> [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
#> [11] "Sigma additive: standard deviation of additive residual error"                      
#> [12] "Sigma proportional: standard deviation of proportional residual error"              
#> 
#> =============== End of Summary ===============
# }
```
