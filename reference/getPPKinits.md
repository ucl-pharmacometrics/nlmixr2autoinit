# Automated pipeline for generating initial estimates in population PK models

Provides a unified and fully automated workflow to generate initial
pharmacokinetic and residual variability parameters for population PK
models using concentration–time data from bolus, infusion, or oral
administration.

## Usage

``` r
getPPKinits(dat, control = initsControl(), verbose = TRUE)
```

## Arguments

- dat:

  A data frame containing pharmacokinetic records in standard nlmixr2
  format, including ID, TIME, EVID, and DV.

- control:

  A list created by
  [`initsControl()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/initsControl.md)
  specifying configuration for pooling, non-compartmental analysis,
  steady-state detection, fallback rules, statistical model components,
  and parameter selection metrics.

- verbose:

  Logical (default = TRUE); when TRUE, displays key progress messages
  and stepwise updates during the initialization process. When FALSE,
  the function runs quietly without printing intermediate information.

## Value

An object of class `getPPKinits` containing recommended initial
parameter estimates, intermediate results, and computation diagnostics.

## Details

The pipeline integrates four model-informed analytical components
applied to raw or pooled concentration–time profiles:

1.  Adaptive single-point methods

2.  Naive pooled graphic methods

3.  Naive pooled non-compartmental analysis (NCA) with optional
    Wagner–Nelson Ka calculation for oral dosing

4.  Parameter sweeping across one-, two-, three-compartment and
    Michaelis–Menten models

In addition to structural PK parameters, the framework also initializes
statistical model components:

- Inter-individual variability (IIV): pragmatic fixed \\\omega^2\\
  values are assigned to random effects.

- Residual unexplained variability (RUV): estimated either using a
  data-driven method based on trimmed residual summaries or a
  fixed-fraction approach consistent with NONMEM User Guide
  recommendations.

- Model applicability: the automated and model-informed strategy
  generates robust initial values suitable for both linear and nonlinear
  mixed-effects pharmacokinetic models.

## See also

[`initsControl`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/initsControl.md),
[`run_single_point`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point.md),
[`run_graphcal`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_graphcal.md),
[`run_pooled_nca`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_pooled_nca.md),
[`sim_sens_1cmpt_mm`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_1cmpt_mm.md),
[`sim_sens_2cmpt`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_2cmpt.md),
[`sim_sens_3cmpt`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_3cmpt.md),
[`metrics.`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/metrics..md)

## Author

Zhonghui Huang

## Examples

``` r
# \donttest{
## Bolus example
getPPKinits(Bolus_1CPT,verbose = TRUE)
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
#> [1] "31.277s"
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
## Oral example (run quietly)
getPPKinits(Oral_1CPT,verbose = FALSE)
#> ===============Initial Parameter Estimation Summary ===============
#> 
#> Recommended initial estimates :
#>            Parameters               Methods    Values
#> 1                  Ka Naive pooled NCA (MD)     1.000
#> 2                  CL Naive pooled NCA (MD)     4.120
#> 3                  Vd Naive pooled NCA (MD)    70.600
#> 4                Vmax    Parameter sweeping 17735.186
#> 5                  Km    Parameter sweeping  4251.513
#> 6           Vc(2CMPT)    Parameter sweeping    66.700
#> 7           Vp(2CMPT)    Parameter sweeping     6.670
#> 8            Q(2CMPT)    Parameter sweeping     1.030
#> 9           Vc(3CMPT)    Parameter sweeping    66.700
#> 10          Vp(3CMPT)    Parameter sweeping     6.670
#> 11         Vp2(3CMPT)    Parameter sweeping     6.670
#> 12           Q(3CMPT)    Parameter sweeping     1.030
#> 13          Q2(3CMPT)    Parameter sweeping     1.030
#> 14     Sigma additive           Model-based     9.529
#> 15 Sigma proportional           Model-based     0.105
#> 
#> Time spent :
#> [1] "54.941s"
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
