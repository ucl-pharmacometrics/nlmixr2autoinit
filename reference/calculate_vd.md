# Calculates volume of distribution from concentration data

Calculates the volume of distribution (Vd) using an adaptive
single-point approach

## Usage

``` r
calculate_vd(
  dat,
  half_life = NULL,
  dose_type = NULL,
  pooled_ctrl = pooled_control(),
  route = c("bolus", "oral", "infusion")
)
```

## Arguments

- dat:

  A data frame containing raw time–concentration data in the standard
  nlmixr2 format.

- half_life:

  Optional numeric value for the drug's half-life. If not provided, it
  will be estimated using
  [`get_hf()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_hf.md)
  from pooled observations.

- dose_type:

  Specifies the dosing context of the pharmacokinetic observations.
  Required when half_life is not provided. Classified as first_dose,
  repeated_doses, or combined_doses based on whether observed
  concentrations occur following the first administration, during
  repeated dosing, or across both contexts.

- pooled_ctrl:

  Optional list of control parameters used by
  [`get_pooled_data()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)
  for pooling observations. Defaults to output from
  [`pooled_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md).

- route:

  Character string specifying the route of administration. Must be one
  of bolus, oral, or infusion. Currently, oral is not implemented.

## Value

A list with two elements:

- vd_df: individual volume of distribution estimates

- trimmed_mean_vd: population volume of distribution estimated as a
  trimmed geometric mean using a 5 percent trimming level

## Details

The function uses a concentration observed within the first 20% of the
elimination half-life after dosing as the early point for estimating the
volume of distribution.

\$\$Vd = \frac{\text{Dose}}{C_0}\$\$ For infusion: \$\$Vd =
\frac{\text{Rate} \times \min(\text{TIME}, \text{durationobs})}{C_0}\$\$

Here, \\C_0\\ represents the early concentration observed within the
first 20% of the elimination half-life after dosing, which is used as an
approximation of the initial concentration for estimating volume of
distribution (Vd). `TIME` refers to time after dose; `durationobs` is
the actual infusion duration.

When half_life is not provided, it is estimated from pooled data using
the functions
[`get_pooled_data()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)
and
[`get_hf()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_hf.md).

## See also

[get_pooled_data](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md),
[get_hf](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_hf.md),
[trimmed_geom_mean](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trimmed_geom_mean.md)

## Author

Zhonghui Huang

## Examples

``` r

dat <- Bolus_1CPT
out <- processData(dat)
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
fdat<- out$dat
froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
half_life <- get_hf(dat = fdat)$half_life_median
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 11 h
calculate_vd(dat = fdat, half_life = half_life,route=froute)$trimmed_mean_vd
#> [1] 67.57943
```
