# nlmixr2autoinit

nlmixr2autoinit provides an automated framework for deriving initial
parameter estimates for population pharmacokinetic (PopPK) models via
multiple PK analysis methods.

It aims to streamline model parameter initialization in R, reducing
reliance on manual trial-and-error and facilitating automated model
development workflows.

## Installation

You can install nlmixr2autoinit either from CRAN or from the development
version on GitHub:

From CRAN (stable release)

``` r

install.packages("nlmixr2autoinit", dependencies = TRUE)
```

From GitHub (development version)

``` r

library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
```

For additional information on installing the broader nlmixr2 ecosystem
please see the installation instructions in the main nlmixr2 repository
(<https://github.com/nlmixr2/nlmixr2>). A nlmixr2autoinit tutorial is
available at <https://ucl-pharmacometrics.github.io/nlmixr2autoinit>.

## Example 1 (IV case)

``` r

library(nlmixr2autoinit)
inits.out<-getPPKinits(dat = Bolus_1CPT)

# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        120            
# Number of Observations                    6951           
# Subjects with First-Dose Interval Data    120            
# Observations in the First-Dose Interval   2276           
# Subjects with Multiple-Dose Data          120            
# Observations after Multiple Doses         4675           
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 11.26 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 4, estimated Vd: 66 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# > inits.out                                                                     
# 
# ===============Initial Parameter Estimation Summary ===============
# 
# Recommended initial estimates :
#            Parameters               Methods    Values
# 1                  Ka                    IV        NA
# 2                  CL Naive pooled NCA (MD)     4.000
# 3                  Vd Naive pooled NCA (MD)    66.000
# 4                Vmax    Parameter sweeping 21551.820
# 5                  Km    Parameter sweeping  5321.437
# 6           Vc(2CMPT)    Parameter sweeping    66.000
# 7           Vp(2CMPT)    Parameter sweeping     6.600
# 8            Q(2CMPT)    Parameter sweeping     1.000
# 9           Vc(3CMPT)    Parameter sweeping    47.601
# 10          Vp(3CMPT)    Parameter sweeping     9.520
# 11         Vp2(3CMPT)    Parameter sweeping     9.520
# 12           Q(3CMPT)    Parameter sweeping     8.000
# 13          Q2(3CMPT)    Parameter sweeping     8.000
# 14     Sigma additive           Model-based    11.111
# 15 Sigma proportional           Model-based     0.114
# 
# Time spent :
# [1] "44.972s"
# 
# ETA variances and derived covariances:
#         Parameters                  Methods Values
# 1           eta.ka             fixed_values    0.1
# 2           eta.cl             fixed_values    0.1
# 3           eta.vc             fixed_values    0.1
# 4           eta.vp             fixed_values    0.1
# 5            eta.q             fixed_values    0.1
# 6          eta.vp2             fixed_values    0.1
# 7           eta.q2             fixed_values    0.1
# 8         eta.vmax             fixed_values    0.1
# 9           eta.km             fixed_values    0.1
# 10 cor.eta_vmax_km eta_corr_derived (r=0.1)   0.01
# 11   cor.eta_cl_vc eta_corr_derived (r=0.1)   0.01
# 12   cor.eta_cl_vp eta_corr_derived (r=0.1)   0.01
# 13  cor.eta_cl_vp2 eta_corr_derived (r=0.1)   0.01
# 14    cor.eta_cl_q eta_corr_derived (r=0.1)   0.01
# 15   cor.eta_cl_q2 eta_corr_derived (r=0.1)   0.01
# 16   cor.eta_vc_vp eta_corr_derived (r=0.1)   0.01
# 17  cor.eta_vc_vp2 eta_corr_derived (r=0.1)   0.01
# 18    cor.eta_vc_q eta_corr_derived (r=0.1)   0.01
# 19   cor.eta_vc_q2 eta_corr_derived (r=0.1)   0.01
# 20  cor.eta_vp_vp2 eta_corr_derived (r=0.1)   0.01
# 21    cor.eta_vp_q eta_corr_derived (r=0.1)   0.01
# 22   cor.eta_vp_q2 eta_corr_derived (r=0.1)   0.01
# 23   cor.eta_vp2_q eta_corr_derived (r=0.1)   0.01
# 24  cor.eta_vp2_q2 eta_corr_derived (r=0.1)   0.01
# 25    cor.eta_q_q2 eta_corr_derived (r=0.1)   0.01
# Note: The ETA variances and covariances listed above are predefined default initialization values automatically assigned by the package.
# 
# Parameter descriptions:
#  [1] "Ka: absorption constant rate"                                                       
#  [2] "CL: clearance"                                                                      
#  [3] "Vd: volume of distribution"                                                         
#  [4] "Vmax: maximum metabolic rate"                                                       
#  [5] "Km: Michaelis constant"                                                             
#  [6] "Vc: volume of distribution of the central compartment"                              
#  [7] "Vp: volume of distribution of the peripheral compartment"                           
#  [8] "Vp2: volume of distribution of the second peripheral compartment"                   
#  [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# [11] "Sigma additive: standard deviation of additive residual error"                      
# [12] "Sigma proportional: standard deviation of proportional residual error"              
# 
# =============== End of Summary ===============
```

## Example 2 (Using multiple cores)

By default,
[`getPPKinits()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getPPKinits.md)
runs with 2 cores. You can increase performance by setting `ncores` to
half of your available system cores:

``` r

library(nlmixr2autoinit)
# Use half of the available system cores to speed up parameter sweeping
ncores <- max(1, floor(parallel::detectCores() / 2))
inits.out <- getPPKinits(dat = Bolus_1CPT, ncores = ncores)
inits.out
# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        120            
# Number of Observations                    6951           
# Subjects with First-Dose Interval Data    120            
# Observations in the First-Dose Interval   2276           
# Subjects with Multiple-Dose Data          120            
# Observations after Multiple Doses         4675           
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 11.26 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 4, estimated Vd: 66 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# > inits.out                                                                     
# 
# ===============Initial Parameter Estimation Summary ===============
# 
# Recommended initial estimates :
#            Parameters               Methods    Values
# 1                  Ka                    IV        NA
# 2                  CL Naive pooled NCA (MD)     4.000
# 3                  Vd Naive pooled NCA (MD)    66.000
# 4                Vmax    Parameter sweeping 21551.820
# 5                  Km    Parameter sweeping  5321.437
# 6           Vc(2CMPT)    Parameter sweeping    66.000
# 7           Vp(2CMPT)    Parameter sweeping     6.600
# 8            Q(2CMPT)    Parameter sweeping     1.000
# 9           Vc(3CMPT)    Parameter sweeping    47.601
# 10          Vp(3CMPT)    Parameter sweeping     9.520
# 11         Vp2(3CMPT)    Parameter sweeping     9.520
# 12           Q(3CMPT)    Parameter sweeping     8.000
# 13          Q2(3CMPT)    Parameter sweeping     8.000
# 14     Sigma additive           Model-based    11.111
# 15 Sigma proportional           Model-based     0.114
# 
# Time spent :
# [1] "40.771s"
# 
# ETA variances and derived covariances:
#         Parameters                  Methods Values
# 1           eta.ka             fixed_values    0.1
# 2           eta.cl             fixed_values    0.1
# 3           eta.vc             fixed_values    0.1
# 4           eta.vp             fixed_values    0.1
# 5            eta.q             fixed_values    0.1
# 6          eta.vp2             fixed_values    0.1
# 7           eta.q2             fixed_values    0.1
# 8         eta.vmax             fixed_values    0.1
# 9           eta.km             fixed_values    0.1
# 10 cor.eta_vmax_km eta_corr_derived (r=0.1)   0.01
# 11   cor.eta_cl_vc eta_corr_derived (r=0.1)   0.01
# 12   cor.eta_cl_vp eta_corr_derived (r=0.1)   0.01
# 13  cor.eta_cl_vp2 eta_corr_derived (r=0.1)   0.01
# 14    cor.eta_cl_q eta_corr_derived (r=0.1)   0.01
# 15   cor.eta_cl_q2 eta_corr_derived (r=0.1)   0.01
# 16   cor.eta_vc_vp eta_corr_derived (r=0.1)   0.01
# 17  cor.eta_vc_vp2 eta_corr_derived (r=0.1)   0.01
# 18    cor.eta_vc_q eta_corr_derived (r=0.1)   0.01
# 19   cor.eta_vc_q2 eta_corr_derived (r=0.1)   0.01
# 20  cor.eta_vp_vp2 eta_corr_derived (r=0.1)   0.01
# 21    cor.eta_vp_q eta_corr_derived (r=0.1)   0.01
# 22   cor.eta_vp_q2 eta_corr_derived (r=0.1)   0.01
# 23   cor.eta_vp2_q eta_corr_derived (r=0.1)   0.01
# 24  cor.eta_vp2_q2 eta_corr_derived (r=0.1)   0.01
# 25    cor.eta_q_q2 eta_corr_derived (r=0.1)   0.01
# Note: The ETA variances and covariances listed above are predefined default initialization values automatically assigned by the package.
# 
# Parameter descriptions:
#  [1] "Ka: absorption constant rate"                                                       
#  [2] "CL: clearance"                                                                      
#  [3] "Vd: volume of distribution"                                                         
#  [4] "Vmax: maximum metabolic rate"                                                       
#  [5] "Km: Michaelis constant"                                                             
#  [6] "Vc: volume of distribution of the central compartment"                              
#  [7] "Vp: volume of distribution of the peripheral compartment"                           
#  [8] "Vp2: volume of distribution of the second peripheral compartment"                   
#  [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# [11] "Sigma additive: standard deviation of additive residual error"                      
# [12] "Sigma proportional: standard deviation of proportional residual error"              
```

## Example 3 (Oral case)

``` r

inits.out<-getPPKinits(dat = Oral_1CPT,ncores=10)
inits.out
# ===============Initial Parameter Estimation Summary ===============
# 
# Recommended initial estimates :
#            Parameters               Methods    Values
# 1                  Ka Naive pooled NCA (MD)     1.000
# 2                  CL Naive pooled NCA (MD)     4.120
# 3                  Vd Naive pooled NCA (MD)    70.600
# 4                Vmax    Parameter sweeping 17735.186
# 5                  Km    Parameter sweeping  4251.513
# 6           Vc(2CMPT)    Parameter sweeping    66.700
# 7           Vp(2CMPT)    Parameter sweeping     6.670
# 8            Q(2CMPT)    Parameter sweeping     1.030
# 9           Vc(3CMPT)    Parameter sweeping    66.700
# 10          Vp(3CMPT)    Parameter sweeping     6.670
# 11         Vp2(3CMPT)    Parameter sweeping     6.670
# 12           Q(3CMPT)    Parameter sweeping     1.030
# 13          Q2(3CMPT)    Parameter sweeping     1.030
# 14     Sigma additive           Model-based     9.529
# 15 Sigma proportional           Model-based     0.105
# 
# Time spent :
# [1] "58.415s"
# 
# ETA variances and derived covariances:
#         Parameters                  Methods Values
# 1           eta.ka             fixed_values    0.1
# 2           eta.cl             fixed_values    0.1
# 3           eta.vc             fixed_values    0.1
# 4           eta.vp             fixed_values    0.1
# 5            eta.q             fixed_values    0.1
# 6          eta.vp2             fixed_values    0.1
# 7           eta.q2             fixed_values    0.1
# 8         eta.vmax             fixed_values    0.1
# 9           eta.km             fixed_values    0.1
# 10 cor.eta_vmax_km eta_corr_derived (r=0.1)   0.01
# 11   cor.eta_cl_vc eta_corr_derived (r=0.1)   0.01
# 12   cor.eta_cl_vp eta_corr_derived (r=0.1)   0.01
# 13  cor.eta_cl_vp2 eta_corr_derived (r=0.1)   0.01
# 14    cor.eta_cl_q eta_corr_derived (r=0.1)   0.01
# 15   cor.eta_cl_q2 eta_corr_derived (r=0.1)   0.01
# 16   cor.eta_vc_vp eta_corr_derived (r=0.1)   0.01
# 17  cor.eta_vc_vp2 eta_corr_derived (r=0.1)   0.01
# 18    cor.eta_vc_q eta_corr_derived (r=0.1)   0.01
# 19   cor.eta_vc_q2 eta_corr_derived (r=0.1)   0.01
# 20  cor.eta_vp_vp2 eta_corr_derived (r=0.1)   0.01
# 21    cor.eta_vp_q eta_corr_derived (r=0.1)   0.01
# 22   cor.eta_vp_q2 eta_corr_derived (r=0.1)   0.01
# 23   cor.eta_vp2_q eta_corr_derived (r=0.1)   0.01
# 24  cor.eta_vp2_q2 eta_corr_derived (r=0.1)   0.01
# 25    cor.eta_q_q2 eta_corr_derived (r=0.1)   0.01
# Note: The ETA variances and covariances listed above are predefined default initialization values automatically assigned by the package.
# 
# Parameter descriptions:
#  [1] "Ka: absorption constant rate"                                                       
#  [2] "CL: clearance"                                                                      
#  [3] "Vd: volume of distribution"                                                         
#  [4] "Vmax: maximum metabolic rate"                                                       
#  [5] "Km: Michaelis constant"                                                             
#  [6] "Vc: volume of distribution of the central compartment"                              
#  [7] "Vp: volume of distribution of the peripheral compartment"                           
#  [8] "Vp2: volume of distribution of the second peripheral compartment"                   
#  [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# [11] "Sigma additive: standard deviation of additive residual error"                      
# [12] "Sigma proportional: standard deviation of proportional residual error"              
# 
# =============== End of Summary ===============
```
