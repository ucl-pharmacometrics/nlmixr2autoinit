# nlmixr2autoinit
Automated Generation of Initial Estimates for Population Pharmacokinetic Modelling

Installation nlmixr2autoinit and nlmixr2auto:
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
install_github("ucl-pharmacometrics/nlmixr2auto")
```
Example
``` r
library(nlmixr2autoinit)
d1<-pheno_sd
getppkinit(dat = d1,runnpd = 0)
# Warning: Half-life not provided. Change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered.$Datainfo
# [1] "No. of subjects: 59, No. of observations: 155, Is infusion? N, Is single dose? N, First-dose data available? Y"
# 
# $Recommended_initial_estimates
# Parameters                             Methods     Values
# 1         CL              Simplified calculation   0.009410
# 2         Vd              Simplified calculation   1.270000
# 3       Vmax Sensitivity analysis by simulation    2.126512
# 4         Km Sensitivity analysis by simulation  220.472441
# 5  Vc(2CMPT) Sensitivity analysis by simulation    1.150000
# 6  Vp(2CMPT) Sensitivity analysis by simulation    0.115000
# 7  Vc(3CMPT) Sensitivity analysis by simulation    1.020000
# 8  Vp(3CMPT) Sensitivity analysis by simulation    0.205000
# 9 Vp2(3CMPT) Sensitivity analysis by simulation    0.041000
# 
# $Message
# NULL
# 
# $Run.history
# $Run.history$base.out
# Method Calculated CL Calculated Vd Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1        Simplified calculation      0.009410        1.2700                        1499.937                                    38 0.1313 secs
# 2           Graphic calculation      0.127000        0.6920                        3264.040                                    79 0.0041 secs
# 3         NCA (only first dose)      0.144000        0.7820                        3200.957                                  1422 0.0084 secs
# 4              NCA (all pooled)      0.000481        0.0994                       54222.785                                    77 0.0029 secs
# 5 Hybrid simplified calculation      0.009410        0.0511                       20239.870                                   617 0.1313 secs
# 
# $Run.history$sim.vmax.km
# Simulated Vmax Simulated Km Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1      2.12651181   220.472441                        1565.268                              61.42196 0.4299 secs
# 2      1.08918898   110.236220                        1615.050                              63.11101 0.5519 secs
# 3      0.57052756    55.118110                        1675.077                              65.38328 0.6660 secs
# 4      0.31119685    27.559055                        1773.154                              70.13415 0.7770 secs
# 5      0.18153150    13.779528                        1904.442                              77.37268 0.8834 secs
# 6      0.11669882     6.889764                        2040.223                              85.01391 1.0027 secs
# 7      0.10373228     5.511811                        2081.696                              87.25051 1.1103 secs
# 8      0.07779921     2.755906                        2174.681                              92.17428 1.2276 secs
# 9      2.20431102   220.472441                        1568.976                              61.55170 1.3423 secs
# 10     1.16698819   110.236220                        1622.008                              63.36649 1.4565 secs
# 11     0.64832677    55.118110                        1685.555                              65.62749 1.5679 secs
# 12     0.38899606    27.559055                        1749.288                              68.11595 1.6915 secs
# 13     0.25933071    13.779528                        1812.801                              70.78504 1.8127 secs
# 14     0.19449803     6.889764                        1866.408                              73.19215 1.9267 secs
# 15     0.18153150     5.511811                        1883.686                              74.10325 2.0383 secs
# 16     0.15559843     2.755906                        1920.674                              76.05906 2.1603 secs
# 17     2.33397638   220.472441                        1579.479                              61.95466 2.2724 secs
# 18     1.29665354   110.236220                        1641.569                              64.09946 2.3841 secs
# 19     0.77799213    55.118110                        1724.017                              67.02896 2.4936 secs
# 20     0.51866142    27.559055                        1819.272                              70.54786 2.6017 secs
# 21     0.38899606    13.779528                        1927.569                              74.63918 2.7194 secs
# 22     0.32416339     6.889764                        2015.615                              78.00270 2.8539 secs
# 23     0.31119685     5.511811                        2041.935                              78.99525 2.9584 secs
# 24     0.28526378     2.755906                        2115.960                              81.72119 3.0709 secs
# 25     2.46364173   220.472441                        1594.971                              62.58957 3.1773 secs
# 26     1.42631890   110.236220                        1678.217                              65.54361 3.2780 secs
# 27     0.90765748    55.118110                        1792.052                              69.59923 3.5134 secs
# 28     0.64832677    27.559055                        1958.472                              75.62486 3.6153 secs
# 29     0.51866142    13.779528                        2147.157                              82.41616 3.7315 secs
# 30     0.45382874     6.889764                        2333.111                              88.94069 3.8359 secs
# 31     0.44086220     5.511811                        2384.674                              90.72025 3.9484 secs
# 32     0.41492913     2.755906                        2529.903                              95.69928 4.0521 secs
# 
# $Run.history$sim.2cmpt
# Simulated Vc Simulated Vp Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1        1.150        0.115                        1550.253                                  39.9 0.4607 secs
# 2        1.060        0.212                        1616.998                                  42.1 0.5709 secs
# 3        0.847        0.423                        1908.218                                  51.6 0.6885 secs
# 4        0.635        0.635                        2454.753                                  69.3 0.7938 secs
# 5        0.423        0.847                        3544.244                                 103.6 0.9038 secs
# 6        0.212        1.060                        6190.860                                 185.8 1.0116 secs
# 7        0.115        1.150                        9394.424                                 283.2 1.1132 secs
# 
# $Run.history$sim.3cmpt
# Simulated Vc Simulated Vp Simulated Vp2 Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1        1.0200        0.205        0.0410                        1526.905                                  39.3 0.4553 secs
# 2        0.7940        0.397        0.0794                        1553.006                                  40.2 0.5851 secs
# 3        0.5770        0.577        0.1150                        1588.827                                  41.4 0.7025 secs
# 4        0.3740        0.747        0.1490                        1619.240                                  42.5 0.8230 secs
# 5        0.2190        0.876        0.1750                        1647.325                                  43.4 0.9400 secs
# 6        0.9770        0.195        0.0977                        1567.835                                  40.7 1.0587 secs
# 7        0.7260        0.363        0.1810                        1655.515                                  43.7 1.1785 secs
# 8        0.5080        0.508        0.2540                        1752.666                                  47.0 1.3002 secs
# 9        0.3180        0.635        0.3180                        1847.408                                  50.2 1.4205 secs
# 10       0.1810        0.726        0.3630                        1935.684                                  53.2 1.5390 secs
# 11       0.9070        0.181        0.1810                        1660.857                                  43.9 1.6666 secs
# 12       0.6350        0.318        0.3180                        1848.474                                  50.3 1.8225 secs
# 13       0.4230        0.423        0.4230                        2076.634                                  57.9 1.9384 secs
# 14       0.2540        0.508        0.5080                        2326.503                                  65.9 2.0614 secs
# 15       0.1410        0.564        0.5640                        2516.009                                  72.0 2.1741 secs
# 16       0.7940        0.159        0.3180                        1850.363                                  50.3 2.3025 secs
# 17       0.5080        0.254        0.5080                        2328.175                                  66.0 2.4164 secs
# 18       0.3180        0.318        0.6350                        2830.351                                  81.9 2.5357 secs
# 19       0.1810        0.363        0.7260                        3300.779                                  96.7 2.6706 secs
# 20       0.0977        0.391        0.7820                        3674.837                                 108.4 2.7911 secs
# 21       0.6350        0.127        0.5080                        2314.308                                  65.6 2.9215 secs
# 22       0.3630        0.181        0.7260                        3303.177                                  96.7 3.0355 secs
# 23       0.2120        0.212        0.8470                        4183.407                                 124.4 3.1593 secs
# 24       0.1150        0.231        0.9240                        5005.234                                 150.3 3.2777 secs
# 25       0.0605        0.242        0.9680                        5590.756                                 168.8 3.4131 secs
# 
# $Run.history$npd.out.vmax.km
# [1] "No model fitting by naive pool data approach was conducted"
# 
# $Run.history$npd.out.2cmpt
# [1] "No model fitting by naive pool data approach was conducted"
# 
# $Run.history$npd.out.3cmpt
# [1] "No model fitting by naive pool data approach was conducted"
# 
# 
# $Parameter.descriptions
# [1] "CL: clearance"                                                                      
# [2] "Vd: volumn of distribution"                                                         
# [3] "Vmax: maximum metobolic rate"                                                       
# [4] "Km: Michaelis constant"                                                             
# [5] "Vc: volume of distribution of the central compartment"                              
# [6] "Vp: volume of distribution of the peripheral compartment"                           
# [7] "Vp: volume of distribution of the second peripheral compartment"                    
# [8] "Q: inter-compartmental clearance"                                                   
# [9] "Q2: inter-compartmental clearance between central and second peripheral compartment"
``` 
