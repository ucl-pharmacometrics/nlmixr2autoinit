# nlmixr2autoinit
Automated Generation of Initial Estimates for Population Pharmacokinetic Modelling

## Installation
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
```
## Example 1 (IV case)
``` r
library(nlmixr2autoinit)
inits.out<-getPPKinits(dat = Bolus_1CPT)
# Infometrics                               Value          
# ----------------------------------------  ---------------
#  Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Total Number of Subjects                  120            
# Total Number of Observations              6951           
# Subjects with First-Dose Interval Data    120            
# Observations in the First-Dose Interval   2276           
# Subjects with Multiple-Dose Data          120            
# Observations after Multiple Doses         4675           
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t½ = 11.26 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# (hybrid mode: parameters combined across sources)....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 4, estimated Vd: 66 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters.................... 
# 
# > inits.out
# ===============Initial Parameter Estimation Summary ===============
#   
#   Recommended initial estimates :
#   Parameters               Methods    Values
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
# [1] "53.291s"
# 
# ETA variances and derived covariances:
#   Parameters                  Methods Values
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
# 
# Parameter descriptions:
# [1] "Ka: absorption constant rate"                                                       
# [2] "CL: clearance"                                                                      
# [3] "Vd: volumn of distribution"                                                         
# [4] "Vmax: maximum metobolic rate"                                                       
# [5] "Km: Michaelis constant"                                                             
# [6] "Vc: volume of distribution of the central compartment"                              
# [7] "Vp: volume of distribution of the peripheral compartment"                           
# [8] "Vp: volume of distribution of the second peripheral compartment"                    
# [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# [11] "Sigma additive: stanadard deviation of additive residual additive error"            
# [12] "Sigma proportional: stanadard deviation of proportional residual additive error"    
# 
# =============== End of Summary ===============
```

## Example 2 (Oral case)
``` r
library(nlmixr2autoinit)
inits.out<-getPPKinits(dat = Oral_1CPT)

# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                oral           
# Dose Type                                 combined_doses 
# Total Number of Subjects                  120            
# Total Number of Observations              6947           
# Subjects with First-Dose Interval Data    120            
# Observations in the First-Dose Interval   2273           
# Subjects with Multiple-Dose Data          120            
# Observations after Multiple Doses         4674           
# ----------------------------------------  ------
#   Estimating half-life....................
# Half-life estimation complete: Estimated t½ = 11.65 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# (hybrid mode: parameters combined across sources)....................
# Base PK parameter analysis finished. Estimated ka: 0.888, estimated CL: 4.2, estimated Vd: 66.7 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters.................... 
# > inits.out                                                                     
# ===============Initial Parameter Estimation Summary ===============
#   
#   Recommended initial estimates :
#   Parameters            Methods    Values
# 1                  Ka    Graphic methods     0.888
# 2                  CL    Graphic methods     4.200
# 3                  Vd    Graphic methods    66.700
# 4                Vmax Parameter sweeping 18079.559
# 5                  Km Parameter sweeping  4251.513
# 6           Vc(2CMPT) Parameter sweeping    57.983
# 7           Vp(2CMPT) Parameter sweeping    11.597
# 8            Q(2CMPT) Parameter sweeping     4.200
# 9           Vc(3CMPT) Parameter sweeping    57.983
# 10          Vp(3CMPT) Parameter sweeping     5.798
# 11         Vp2(3CMPT) Parameter sweeping     5.798
# 12           Q(3CMPT) Parameter sweeping     2.100
# 13          Q2(3CMPT) Parameter sweeping     2.100
# 14     Sigma additive        Model-based     9.529
# 15 Sigma proportional        Model-based     0.105
# 
# Time spent :
# [1] "39.635s"
# 
# ETA variances and derived covariances:
#   Parameters                  Methods Values
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
# 
# Parameter descriptions:
# [1] "Ka: absorption constant rate"                                                       
# [2] "CL: clearance"                                                                      
# [3] "Vd: volumn of distribution"                                                         
# [4] "Vmax: maximum metobolic rate"                                                       
# [5] "Km: Michaelis constant"                                                             
# [6] "Vc: volume of distribution of the central compartment"                              
# [7] "Vp: volume of distribution of the peripheral compartment"                           
# [8] "Vp: volume of distribution of the second peripheral compartment"                    
# [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# [11] "Sigma additive: stanadard deviation of additive residual additive error"            
# [12] "Sigma proportional: stanadard deviation of proportional residual additive error"    
# 
# =============== End of Summary ===============
```
