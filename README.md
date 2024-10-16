# nlmixr2autoinit
Automated Generation of Initial Estimates for Population Pharmacokinetic Modelling

Installation nlmixr2autoinit
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
```
IV Case Example
``` r
library(nlmixr2autoinit)
inits.out<-getppkinit(dat = Bolus_1CPT)
# Processing data....................
# Estimating half-life....................
# Estimated half-life : 12.19
# Run approximate PK calculation on individual data....................
# Try approximate calculation using estimated half-life: 12.19....................
# Try approximate calculation without no prior estimated half-life....................
# No half life detected, change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered.
# Run non-compartmental analysis on naive pooling data....................
# Run graphical analysis on naive pooling data only after the first dose....................
# Run hybrid calculation, combining approximate calculation methods and NCA ....................
# Base PK parameter analysis finished. Estimated ka :NA, estimated CL : 3.57, estimated Vd : 53.2
# Run sensitivity analysis on nonlinear eliminiation kinetics PK parameters....................
# Run sensitivity analysis on multi-compartmental PK parameters....................
# ===============Initial Parameter Estimation Summary ===============
#   Data information:
#   [1] "No. of subjects: 120, No. of observations: 6951, Dose route: bolus, Is single dose? N, Data within the first dosing interval available? Y"
# 
# Recommended initial estimates by using non-model fitting methods:
#   Parameters                             Methods    Values
# 1          Ka                       Wanger_nelson        NA
# 2          CL      Hybrid approximate calculation     3.570
# 3          Vd      Hybrid approximate calculation    53.200
# 4        Vmax Sensitivity analysis by simulation  14725.358
# 5          Km Sensitivity analysis by simulation   2999.818
# 6   Vc(2CMPT) Sensitivity analysis by simulation      4.840
# 7   Vp(2CMPT) Sensitivity analysis by simulation     48.400
# 8    Q(2CMPT) Sensitivity analysis by simulation    100.000
# 9   Vc(3CMPT) Sensitivity analysis by simulation      7.600
# 10  Vp(3CMPT) Sensitivity analysis by simulation     38.000
# 11 Vp2(3CMPT) Sensitivity analysis by simulation      7.600
# 12   Q(3CMPT) Sensitivity analysis by simulation    100.000
# 13  Q2(3CMPT) Sensitivity analysis by simulation    100.000
# 
# Recommended initial estimates by naive pooled data compartmental analysis:
#   [1] "No model fitting by naive pool data approach was conducted"
# 
# Parameter descriptions:
#   [1] "CL: clearance"                                                                      
# [2] "Vd: volumn of distribution"                                                         
# [3] "Vmax: maximum metobolic rate"                                                       
# [4] "Km: Michaelis constant"                                                             
# [5] "Vc: volume of distribution of the central compartment"                              
# [6] "Vp: volume of distribution of the peripheral compartment"                           
# [7] "Vp: volume of distribution of the second peripheral compartment"                    
# [8] "Q: inter-compartmental clearance"                                                   
# [9] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# 
# =============== End of Summary ===============
```
Oral case example

``` r
library(nlmixr2autoinit)
inits.out<-getppkinit(dat = Oral_1CPT)
# Processing data....................
# Administration site detected to differ from measurement site; extravascular (oral) administration assumed.
# Estimating half-life....................
# Estimated half-life : 12.18
# Run approximate PK calculation on individual data....................
# Try approximate calculation using estimated half-life: 12.18....................
# Try approximate calculation without no prior estimated half-life....................
# No half life detected, change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered.
# Run non-compartmental analysis on naive pooling data....................
# Run graphical analysis on naive pooling data only after the first dose....................
# Run hybrid calculation, combining approximate calculation methods and NCA ....................
# Base PK parameter analysis finished. Estimated ka :0.89, estimated CL : 4.18, estimated Vd : 73.5
# Run sensitivity analysis on nonlinear eliminiation kinetics PK parameters....................
# Run sensitivity analysis on multi-compartmental PK parameters....................
# ===============Initial Parameter Estimation Summary ===============
#   Data information:
#   [1] "No. of subjects: 120, No. of observations: 6947, Dose route: oral, Is single dose? N, Data within the first dosing interval available? Y"
# 
# Recommended initial estimates by using non-model fitting methods:
#   Parameters                             Methods   Values
# 1          Ka                       Wanger_nelson    0.890
# 2          CL      Hybrid approximate calculation    4.180
# 3          Vd      Hybrid approximate calculation   73.500
# 4        Vmax Sensitivity analysis by simulation  8855.661
# 5          Km Sensitivity analysis by simulation  1210.617
# 6   Vc(2CMPT) Sensitivity analysis by simulation    49.000
# 7   Vp(2CMPT) Sensitivity analysis by simulation    24.500
# 8    Q(2CMPT) Sensitivity analysis by simulation     0.100
# 9   Vc(3CMPT) Sensitivity analysis by simulation    29.400
# 10  Vp(3CMPT) Sensitivity analysis by simulation    14.700
# 11 Vp2(3CMPT) Sensitivity analysis by simulation    29.400
# 12   Q(3CMPT) Sensitivity analysis by simulation   100.000
# 13  Q2(3CMPT) Sensitivity analysis by simulation   100.000
# 
# Recommended initial estimates by naive pooled data compartmental analysis:
#   [1] "No model fitting by naive pool data approach was conducted"
# 
# Parameter descriptions:
#   [1] "CL: clearance"                                                                      
# [2] "Vd: volumn of distribution"                                                         
# [3] "Vmax: maximum metobolic rate"                                                       
# [4] "Km: Michaelis constant"                                                             
# [5] "Vc: volume of distribution of the central compartment"                              
# [6] "Vp: volume of distribution of the peripheral compartment"                           
# [7] "Vp: volume of distribution of the second peripheral compartment"                    
# [8] "Q: inter-compartmental clearance"                                                   
# [9] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# 
# =============== End of Summary ===============
```
