# Run graphical analysis of pharmacokinetic parameters

Performs graphical estimation of pharmacokinetic parameters based on
pooled concentration–time data and the specified route of
administration.

## Usage

``` r
run_graphcal(
  dat,
  route,
  dose_type = c("first_dose", "repeated_doses", "combined_doses"),
  pooled = NULL,
  pooled_ctrl = pooled_control(),
  ...
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- route:

  Route of administration. Must be one of bolus, oral, or infusion.

- dose_type:

  Specifies the dosing context of the pharmacokinetic observations.
  Classified as first_dose, repeated_doses, or combined_doses based on
  whether observed concentrations occur following the first
  administration, during repeated dosing, or across both contexts.

- pooled:

  Optional pooled dataset. If NULL, pooling is performed internally.

- pooled_ctrl:

  Control settings created by
  [`pooled_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md)
  for time binning and pooling.

- ...:

  Additional arguments passed to graphical calculation functions.

## Value

A list containing graphical estimates of key pharmacokinetic parameters.

## Details

The function pools individual profiles using
[`get_pooled_data()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)
when needed, and then applies route-specific graphical methods
(`graphcal_iv` or `graphcal_oral`) to estimate parameters such as
clearance, volume of distribution, terminal slope, and absorption rate
constant (for oral data).

## See also

[`graphcal_iv`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_iv.md),
[`graphcal_oral`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_oral.md),
[`get_pooled_data`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)

## Author

Zhonghui Huang

## Examples

``` r
# Example 1 (iv case)
dat <- Bolus_1CPT
dat <- processData(dat)$dat
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
run_graphcal(dat, route="bolus")
#> $cl
#> [1] 4.19866
#> 
#> $vd
#> [1] 66.6356
#> 
#> $slope
#> [1] -0.06300926
#> 
#> $C0exp
#> [1] 0.01500699
#> 
#> $method
#> [1] "find_best_lambdaz"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>     -4.19924      -0.06301  
#> 
#> 
#> $time.spent
#> [1] 0.003
#> 

# Example 2 (oral case)
dat <- Oral_1CPT
dat <- processData(dat)$dat
#> 
#> 
#> Infometrics                               Value          
#> ----------------------------------------  ---------------
#> Dose Route                                oral           
#> Dose Type                                 combined_doses 
#> Number of Subjects                        120            
#> Number of Observations                    6947           
#> Subjects with First-Dose Interval Data    120            
#> Observations in the First-Dose Interval   2273           
#> Subjects with Multiple-Dose Data          120            
#> Observations after Multiple Doses         4674           
#> ----------------------------------------  ------
run_graphcal(dat, route="oral")
#> $ka
#> [1] 0.8880668
#> 
#> $kel
#> [1] 0.06293042
#> 
#> $slope
#> [1] -0.06293042
#> 
#> $C0exp
#> [1] 0.01614338
#> 
#> $cl
#> [1] 4.195522
#> 
#> $vd
#> [1] 66.66922
#> 
#> $method
#> [1] "find_best_lambdaz"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>     -4.12625      -0.06293  
#> 
#> 
#> $time.spent
#> [1] 0.003
#> 

# Example 3 (infusion case).
# Approximate calculation. only use when the infusion duration is very short

dat <- Infusion_1CPT
dat <- processData(dat)$dat
#> 
#> 
#> Infometrics                               Value          
#> ----------------------------------------  ---------------
#> Dose Route                                infusion       
#> Dose Type                                 combined_doses 
#> Number of Subjects                        120            
#> Number of Observations                    6953           
#> Subjects with First-Dose Interval Data    120            
#> Observations in the First-Dose Interval   2276           
#> Subjects with Multiple-Dose Data          120            
#> Observations after Multiple Doses         4677           
#> ----------------------------------------  ------
run_graphcal(dat, route="infusion")
#> $cl
#> [1] 4.219652
#> 
#> $vd
#> [1] 68.13998
#> 
#> $slope
#> [1] -0.06192623
#> 
#> $C0exp
#> [1] 0.01467567
#> 
#> $method
#> [1] "find_best_lambdaz"
#> 
#> $slopefit
#> 
#> Call:
#> lm(formula = log(conc[subset]) ~ time[subset])
#> 
#> Coefficients:
#>  (Intercept)  time[subset]  
#>     -4.22156      -0.06193  
#> 
#> 
#> $time.spent
#> [1] 0.002
#> 
```
