# Automated pipeline for generating initial estimates in population PK models

Provides a unified and fully automated workflow to generate initial
pharmacokinetic and residual variability parameters for population PK
models using concentration–time data from bolus, infusion, or oral
administration.

## Usage

``` r
getPPKinits(dat, control = initsControl(), ncores = 2, verbose = TRUE)
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

- ncores:

  Number of cores to use for parallelization, passed to `rxControl()`.
  Default is 2.

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
getPPKinits(pheno_sd[pheno_sd$ID<11,])
#> 
#> 
#> Infometrics                               Value          
#> ----------------------------------------  ---------------
#> Dose Route                                bolus          
#> Dose Type                                 combined_doses 
#> Number of Subjects                        10             
#> Number of Observations                    30             
#> Subjects with First-Dose Interval Data    10             
#> Observations in the First-Dose Interval   10             
#> Subjects with Multiple-Dose Data          10             
#> Observations after Multiple Doses         20             
#> ----------------------------------------  ------
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 1.47 h
#> Evaluating the predictive performance of calculated one-compartment model parameters....................
#> Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0112, estimated Vd: 0.351 
#> Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
#> Run parameter sweeping on multi-compartmental PK parameters....................
#> ===============Initial Parameter Estimation Summary ===============
#> 
#> Recommended initial estimates :
#>            Parameters                   Methods   Values
#> 1                  Ka                        IV       NA
#> 2                  CL    Naive pooled NCA (all)   0.0112
#> 3                  Vd    Naive pooled NCA (all)   0.3510
#> 4                Vmax        Parameter sweeping   1.2057
#> 5                  Km        Parameter sweeping 119.9200
#> 6           Vc(2CMPT)        Parameter sweeping   0.3510
#> 7           Vp(2CMPT)        Parameter sweeping   0.3510
#> 8            Q(2CMPT)        Parameter sweeping   0.0224
#> 9           Vc(3CMPT)        Parameter sweeping   0.3510
#> 10          Vp(3CMPT)        Parameter sweeping   0.1755
#> 11         Vp2(3CMPT)        Parameter sweeping   0.3510
#> 12           Q(3CMPT)        Parameter sweeping   0.0224
#> 13          Q2(3CMPT)        Parameter sweeping   0.0224
#> 14     Sigma additive Fallback (fixed fraction)   4.8307
#> 15 Sigma proportional Fallback (fixed fraction)   0.2000
#> 
#> Time spent :
#> [1] "52.186s"
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
