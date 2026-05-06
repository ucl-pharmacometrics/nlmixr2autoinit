# Estimate the absorption rate constant using pointwise methods

It implements pointwise estimation of absorption rate constants for
single-dose and multiple-dose pharmacokinetic models.

## Usage

``` r
run_ka_solution(df, cl, ke, Fbio = 1)
```

## Arguments

- df:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- cl:

  A numeric value for drug clearance. It is assumed constant across
  subjects unless pre-specified per subject.

- ke:

  A numeric value for the elimination rate constant. This is assumed
  known or estimated from the terminal phase.

- Fbio:

  A numeric value for bioavailability (F). Default is 1.

## Value

A list with the following elements:

- ka_calc_median:

  Median ka value across all valid observations.

- ka_calc_dat_sd:

  Data frame containing absorption-phase single-dose data, with
  estimated ka and diagnostic messages.

- ka_calc_dat_md:

  Data frame containing absorption-phase steady-state multiple-dose
  data, with ka estimates.

- ka_calc_dat:

  Combined data frame (single-dose and multiple-dose) containing all ka
  estimates.

## Details

For each subject, the time of maximum observed concentration (Tmax) is
identified as the time corresponding to the highest DV. Only records
with TIME less than or equal to Tmax are retained, representing the
absorption phase.

Two scenario-specific calculations are implemented: single-dose and
multiple-dose at steady state.

- For single-dose data (dose_number == 1 and SteadyState == FALSE), the
  function uses
  [`ka_calculation_sd()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_sd.md),
  which applies a one-compartment oral absorption model under
  first-order absorption and elimination.

- For steady-state multiple-dose data (dose_number \> 1 and SteadyState
  == TRUE), the function uses
  [`ka_calculation_md()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_md.md),
  which accounts for accumulation using the dosing interval
  (dose_interval).

This function does not perform model fitting. The median is recommended
for use in pharmacokinetic modeling.

## Author

Zhonghui Huang

## Examples

``` r
# Single-dose
df <- Oral_1CPT[Oral_1CPT$SD == 1, ]
df <- processData(df)$dat
#> 
#> 
#> Infometrics                               Value      
#> ----------------------------------------  -----------
#> Dose Route                                oral       
#> Dose Type                                 first_dose 
#> Number of Subjects                        120        
#> Number of Observations                    2273       
#> Subjects with First-Dose Interval Data    120        
#> Observations in the First-Dose Interval   2273       
#> Subjects with Multiple-Dose Data          0          
#> Observations after Multiple Doses         0          
#> ----------------------------------------  ------
df <- is_ss(df)
run_ka_solution(df = df, cl = 4, ke = 4/70, Fbio = 1)$ka_calc_median
#> [1] 0.8409131

# Mixed doses
dat <- Oral_1CPT
df_ss <- processData(dat)$dat
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
df_ss <- is_ss(df_ss)
run_ka_solution(df = df_ss, cl = 4, ke = 4/70, Fbio = 1)$ka_calc_median
#> [1] 0.6948711
```
