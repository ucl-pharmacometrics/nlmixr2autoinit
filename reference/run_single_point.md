# Run full adaptive single-point PK analysis

Runs full adaptive single-point pharmacokinetic analysis to estimate CL,
Vd, and ka across bolus, oral, and infusion dosing scenarios.

## Usage

``` r
run_single_point(
  dat,
  route = c("bolus", "oral", "infusion"),
  half_life = NULL,
  dose_type = NULL,
  pooled_ctrl = pooled_control(),
  ssctrl = ss_control()
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- route:

  Character string specifying the route of administration. Must be one
  of "bolus", "oral", or "infusion".

- half_life:

  Optional numeric value representing the elimination half-life of the
  drug. If not provided, half-life is estimated within
  run_single_point_base() using get_hf() applied to the pharmacokinetic
  observations.

- dose_type:

  Classified as "first_dose", "repeated_doses", or "combined_doses"
  based on whether observed concentrations occur following the first
  administration, during repeated dosing, or across both contexts. This
  parameter is passed directly to run_single_point_base().

- pooled_ctrl:

  A list of pooled-analysis control options created using
  pooled_control(). These control time binning and time-after-dose
  rounding when pooled summaries are required.

- ssctrl:

  A list of steady-state control options created using ss_control().
  These govern assumptions and thresholds used when interpreting
  steady-state behavior in single-point logic.

## Value

An object of class "single.point.lst" with results from the base and
extended steps.

## See also

[`run_single_point_base`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md),
[`run_single_point_extra`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_extra.md)

## Author

Zhonghui Huang

## Examples

``` r
# Example: Adaptive single-point PK analysis for bolus data
# Step 1: Preprocess the data
dat <- Oral_1CPT
out <- processData(dat)
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
# Step 2: Extract route and dose type info
froute <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
fdose_type <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Type"]
# Step 3: Estimate half life
half_life <- get_hf(dat = out$dat)$half_life_median
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 11.08 h
# Step 4: Run single-point analysis (CL, Vd, Ka if oral)
result <- run_single_point(
  dat = out$dat,
  route = froute,
  dose_type = fdose_type,
  half_life = half_life
)
# Step 5: View results
print(result$singlepoint.results)
#>      ka   cl   vd           starttime time.spent
#> 1 0.889 3.28 52.4 2026-05-10 21:26:51      1.481
#>                                                         single_point.message
#> 1 CL and Vd were calculated directly from steady-state and single-dose data.
```
