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
# Processing data....................
# 
# 
# Infometrics                               Values 
# ----------------------------------------  -------
# Dose Route                                bolus  
# Total Number of Subjects                  120    
# Total Number of Observations              6951   
# Subjects with First-Dose Interval Data    120    
# Observations in the First-Dose Interval   2276   
# Subjects with Multiple-Dose Data          120    
# Observations after Multiple Doses         4675   
# 
# Note: EVID=2 found in the dataset. Rows with EVID=2 are removed.
# ------------------------------------  ------
# Run adaptive single-point method to calculate PK parameters....................
# Estimating half-life....................
# Estimated half-life : 10.39
# Run non-compartmental analysis on naive pooling data....................
# Run graphical analysis on naive pooling data only after the first dose....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 4.05, estimated Vd: 65.2 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# |===============================================================================================================================================| 100%
# Run parameter sweeping on multi-compartmental PK parameters....................
# |===============================================================================================================================================| 100%
# |===============================================================================================================================================| 100%

inits.out
# ===============Initial Parameter Estimation Summary ===============
# Data information:
#   
#   
# Infometrics                               Values 
# ----------------------------------------  -------
# Dose Route                                bolus  
# Total Number of Subjects                  120    
# Total Number of Observations              6951   
# Subjects with First-Dose Interval Data    120    
# Observations in the First-Dose Interval   2276   
# Subjects with Multiple-Dose Data          120    
# Observations after Multiple Doses         4675   
# 
# Note: EVID=2 found in the dataset. Rows with EVID=2 are removed for the entire analysis.
# 
# Recommended initial estimates by using non-model fitting methods:
# Parameters                            Methods    Values
# 1          Ka               Wanger-nelson method        NA
# 2          CL NCA (data exclude first-dose part)     4.050
# 3          Vd NCA (data exclude first-dose part)    65.200
# 4        Vmax                 Parameter sweeping 25817.187
# 5          Km                 Parameter sweeping  5999.637
# 6   Vc(2CMPT)                 Parameter sweeping    65.200
# 7   Vp(2CMPT)                 Parameter sweeping     6.520
# 8    Q(2CMPT)                 Parameter sweeping     1.000
# 9   Vc(3CMPT)                 Parameter sweeping    65.200
# 10  Vp(3CMPT)                 Parameter sweeping     6.520
# 11 Vp2(3CMPT)                 Parameter sweeping    65.200
# 12   Q(3CMPT)                 Parameter sweeping     1.000
# 13  Q2(3CMPT)                 Parameter sweeping     1.000
# 
# Recommended initial estimates by naive pooled data compartmental analysis:
#   [1] "No model fitting by naive pool data approach was conducted"
# 
# Parameter descriptions:
#   [1] "Ka: absorption constant rate"                                                       
# [2] "CL: clearance"                                                                      
# [3] "Vd: volumn of distribution"                                                         
# [4] "Vmax: maximum metobolic rate"                                                       
# [5] "Km: Michaelis constant"                                                             
# [6] "Vc: volume of distribution of the central compartment"                              
# [7] "Vp: volume of distribution of the peripheral compartment"                           
# [8] "Vp: volume of distribution of the second peripheral compartment"                    
# [9] "Q: inter-compartmental clearance"                                                   
# [10] "Q2: inter-compartmental clearance between central and second peripheral compartment"
# 
# =============== End of Summary ===============
```

## Example 2 (Oral case)
``` r
library(nlmixr2autoinit)
inits.out<-getPPKinits(dat = Oral_1CPT)

# Processing data....................
# Administration site detected to differ from measurement site; extravascular (oral) administration assumed.
# 
# 
# Infometrics                               Values 
# ----------------------------------------  -------
# Dose Route                                oral   
# Total Number of Subjects                  120    
# Total Number of Observations              6947   
# Subjects with First-Dose Interval Data    120    
# Observations in the First-Dose Interval   2273   
# Subjects with Multiple-Dose Data          120    
# Observations after Multiple Doses         4674   
# 
# Note: EVID=2 found in the dataset. Rows with EVID=2 are removed for the entire analysis.
# ------------------------------------  ------
# Run adaptive single-point method to calculate PK parameters....................
# Estimating half-life....................
# Estimated half-life : 10.66
# Insufficient data points to support Vd calculation; derived from clearance and estimated half-life instead.....................
# Run non-compartmental analysis on naive pooling data....................
# Run graphical analysis on naive pooling data only after the first dose....................
# Base PK parameter analysis finished. Estimated ka: 0.573, estimated CL: 3.61, estimated Vd: 55.5 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# |==============================================================================================================================================================================| 100%
# Run parameter sweeping on multi-compartmental PK parameters....................
# |================================================================================================================| 100%
# |================================================================================================================| 100%

inits.out
# ===============Initial Parameter Estimation Summary ===============
# Data information:
#   
#   
# Infometrics                               Values 
# ----------------------------------------  -------
# Dose Route                                oral   
# Total Number of Subjects                  120    
# Total Number of Observations              6947   
# Subjects with First-Dose Interval Data    120    
# Observations in the First-Dose Interval   2273   
# Subjects with Multiple-Dose Data          120    
# Observations after Multiple Doses         4674   
# 
# Note: EVID=2 found in the dataset. Rows with EVID=2 are removed for the entire analysis.
# 
# Recommended initial estimates by using non-model fitting methods:
# Parameters              Methods    Values
# 1          Ka Methods of residuals     0.573
# 2          CL     Graphic analysis     3.610
# 3          Vd     Graphic analysis    55.500
# 4        Vmax   Parameter sweeping 19666.468
# 5          Km   Parameter sweeping  4842.467
# 6   Vc(2CMPT)   Parameter sweeping    55.500
# 7   Vp(2CMPT)   Parameter sweeping     5.550
# 8    Q(2CMPT)   Parameter sweeping     1.000
# 9   Vc(3CMPT)   Parameter sweeping    55.500
# 10  Vp(3CMPT)   Parameter sweeping     5.550
# 11 Vp2(3CMPT)   Parameter sweeping    55.500
# 12   Q(3CMPT)   Parameter sweeping     1.000
# 13  Q2(3CMPT)   Parameter sweeping     1.000
# 
# Recommended initial estimates by naive pooled data compartmental analysis:
# [1] "No model fitting by naive pool data approach was conducted"
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
```
