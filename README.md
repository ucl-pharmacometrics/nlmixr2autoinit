# nlmixr2autoinit
Automated Generation of Initial Estimates for Population Pharmacokinetic Modelling

Installation nlmixr2autoinit
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
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
Oral case example
```
library(nlmixr2autoinit)
d1<-Oral_1CPT
getppkinit(dat = d1,runnpd = 0)
# Settings of running nlmixr2autoinit
# 
# Provided reference half-life------------------------------------------------------------------------------- NA
# Number of plasma samples selected for linear regression on terminal phase slope---------------------------- 4
# Trapezoidal rule method----------------------------------------------------------------------------------- 1
# Number of bins during the naive pooled median data processing--------------------------------------------- 8
# Estimated method for naive pooled approach data---------------------------------------------------------- nls
# Non-intravenous case(s) detected, Oral administration was assumed---------------------------------------- 
#   Performed linear regression on the terminal phase of pooled dose-normalized data, estimated half-life: NA
# Run quick calculation with estimated half-life: 12.18
# Run quick calculation without estimated half-life: 
#   Change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered.Base parameter estimation finished. Estimated ka :0.89, estimated CL : 4.18, estimated Vd : 73.5
# $Datainfo
# [1] "No. of subjects: 120, No. of observations: 6947, Is infusion? N, Is single dose? N, First-dose data available? Y, is Oral case? Y"
# 
# $Recommended_initial_estimates
# Parameters                             Methods    Values
# 1          Ka             wanger_nelson_ (median)     0.890
# 2          CL       Hybrid simplified calculation     4.180
# 3          Vd       Hybrid simplified calculation    73.500
# 4        Vmax Sensitivity analysis by simulation  11942.857
# 5          Km Sensitivity analysis by simulation   1632.653
# 6   Vc(2CMPT) Sensitivity analysis by simulation     66.800
# 7   Vp(2CMPT) Sensitivity analysis by simulation      6.680
# 8   Vc(3CMPT) Sensitivity analysis by simulation     45.900
# 9   Vp(3CMPT) Sensitivity analysis by simulation      9.190
# 10 Vp2(3CMPT) Sensitivity analysis by simulation     18.400
# 
# $Message
# NULL
# 
# $Run.history
# $Run.history$base.out
# Method Calculated Ka Calculated CL Calculated Vd Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1             Simplified calculation            NA          4.18            NA                             Inf                                   Inf 3.3159 secs
# 2                Graphic calculation     0.4240000          4.05          71.2                       1086908.5                                    84 0.0047 secs
# 3              NCA (only first dose)     0.8868166          3.99          70.0                        922741.4                                    78 0.0040 secs
# 4 NCA (data exclude first-dose part)     0.9939569          2.95          53.4                       1425113.3                                   122 0.0046 secs
# 5                   NCA (all pooled)     0.3573872          3.31          49.9                       1101035.0                                    80 0.0052 secs
# 6      Hybrid simplified calculation     0.8900000          4.18          73.5                        925834.2                                    74 3.3159 secs
# 
# $Run.history$sim.vmax.km
# Simulated Vmax Simulated Km Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1       27980.408   6530.61224                        951595.3                              5565.710 0.1618 secs
# 2       14331.429   3265.30612                       1098983.8                              6522.694 0.3177 secs
# 3        7506.939   1632.65306                       2132256.4                             13750.077 0.4996 secs
# 4        4094.694    816.32653                       4033459.0                              7477.732 0.6581 secs
# 5        2388.571    408.16327                       4424244.4                              7810.829 0.8167 secs
# 6        1535.510    204.08163                       5030780.9                              8626.398 0.9851 secs
# 7        1364.898    163.26531                       5158936.8                              8800.526 1.1403 secs
# 8        1023.673     81.63265                       5711261.4                              9418.194 1.3018 secs
# 9       29004.082   6530.61224                        935220.4                              4970.290 1.4701 secs
# 10      15355.102   3265.30612                       1006146.2                              5204.291 1.7504 secs
# 11       8530.612   1632.65306                       1416037.7                              6826.666 1.9259 secs
# 12       5118.367    816.32653                       3394329.4                              6737.945 2.0810 secs
# 13       3412.245    408.16327                       4485674.2                              7752.525 2.2447 secs
# 14       2559.184    204.08163                       4302926.1                              7648.851 2.4034 secs
# 15       2388.571    163.26531                       4407828.0                              7781.501 2.6718 secs
# 16       2047.347     81.63265                       4657926.9                              8116.441 2.8417 secs
# 17      30710.204   6530.61224                        931918.6                              4362.223 3.0065 secs
# 18      17061.224   3265.30612                        952640.8                              3950.194 3.1823 secs
# 19      10236.735   1632.65306                       1035960.2                              3616.722 3.3544 secs
# 20       6824.490    816.32653                       1366089.6                              4284.614 3.5221 secs
# 21       5118.367    408.16327                       3352264.3                              6705.271 3.6836 secs
# 22       4265.306    204.08163                       3856210.3                              7195.465 3.9441 secs
# 23       4094.694    163.26531                       3975330.2                              7310.829 4.1034 secs
# 24       3753.469     81.63265                       4231412.5                              7559.313 4.2493 secs
# 25      32416.327   6530.61224                        947329.1                              3805.215 4.4114 secs
# 26      18767.347   3265.30612                        978652.8                              3261.544 4.5885 secs
# 27      11942.857   1632.65306                       1056757.5                              2955.124 4.7451 secs
# 28       8530.612    816.32653                       1191361.5                              3354.992 5.0094 secs
# 29       6824.490    408.16327                       1342934.8                              3904.319 5.1666 secs
# 30       5971.429    204.08163                       1457355.4                              4241.281 5.3179 secs
# 31       5800.816    163.26531                       1507334.7                              4330.352 5.4789 secs
# 32       5459.592     81.63265                       1627005.2                              4552.210 5.7563 secs
# 
# $Run.history$sim.2cmpt
# Simulated Vc Simulated Vp Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1        66.80         6.68                        919033.1                                  74.4 0.1721 secs
# 2        61.20        12.20                        921137.6                                  75.8 0.3306 secs
# 3        49.00        24.50                       1001457.9                                  84.8 0.4771 secs
# 4        36.80        36.80                       1265321.1                                 100.8 0.6231 secs
# 5        24.50        49.00                       1933480.1                                 125.8 0.7686 secs
# 6        12.30        61.20                       3406427.5                                 173.6 0.9294 secs
# 7         6.68        66.80                       4626678.7                                 215.8 1.0767 secs
# 
# $Run.history$sim.3cmpt
# Simulated Vc Simulated Vp Simulated Vp2 Absolute prediction error (APE) Mean absolute prediction error (MAPE)  Time spent
# 1         59.30        11.90          2.37                         1413220                                  89.4 0.2086 secs
# 2         45.90        23.00          4.59                         1699182                                  98.8 0.4114 secs
# 3         33.40        33.40          6.68                         2177591                                 115.3 0.6174 secs
# 4         21.60        43.20          8.65                         2943860                                 144.2 0.8120 secs
# 5         12.70        50.70         10.10                         3832827                                 181.8 1.0272 secs
# 6         56.50        11.30          5.65                         1463454                                  85.7 1.2290 secs
# 7         42.00        21.00         10.50                         1819565                                  93.8 1.4370 secs
# 8         29.40        29.40         14.70                         2393269                                 110.9 1.6365 secs
# 9         18.40        36.80         18.40                         3222846                                 140.5 1.8473 secs
# 10        10.50        42.00         21.00                         4087036                                 177.0 2.0522 secs
# 11        52.50        10.50         10.50                         1545589                                  82.1 2.2571 secs
# 12        36.80        18.40         18.40                         2028070                                  91.3 2.4570 secs
# 13        24.50        24.50         24.50                         2736667                                 112.1 2.6627 secs
# 14        14.70        29.40         29.40                         3640429                                 144.9 2.8714 secs
# 15         8.17        32.70         32.70                         4400369                                 180.5 3.0825 secs
# 16        45.90         9.19         18.40                         1738391                                  79.6 3.3898 secs
# 17        29.40        14.70         29.40                         2496958                                  97.8 3.5817 secs
# 18        18.40        18.40         36.80                         3420177                                 129.3 3.7886 secs
# 19        10.50        21.00         42.00                         4386271                                 170.0 4.1017 secs
# 20         5.65        22.60         45.20                         4885716                                 201.3 4.3192 secs
# 21        36.80         7.35         29.40                         2182821                                  85.7 4.5097 secs
# 22        21.00        10.50         42.00                         3460693                                 127.7 4.7169 secs
# 23        12.20        12.20         49.00                         4692050                                 175.7 4.9974 secs
# 24         6.68        13.40         53.50                         5592591                                 219.2 5.2128 secs
# 25         3.50        14.00         56.00                         5701954                                 235.7 5.4225 secs
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
# [1] "CL: clearance"                                                                       "Vd: volumn of distribution"                                                         
# [3] "Vmax: maximum metobolic rate"                                                        "Km: Michaelis constant"                                                             
# [5] "Vc: volume of distribution of the central compartment"                               "Vp: volume of distribution of the peripheral compartment"                           
# [7] "Vp: volume of distribution of the second peripheral compartment"                     "Q: inter-compartmental clearance"                                                   
# [9] "Q2: inter-compartmental clearance between central and second peripheral compartment"
```
