# Approximate volume of distribution from observed Cmax

Estimates the volume of distribution (Vd) from observed peak
concentrations (Cmax) in single-dose, multiple-dose, or mixed datasets.

## Usage

``` r
approx.vc(
  dat = NULL,
  half_life = NULL,
  single_point_base.lst = NULL,
  route = c("bolus", "oral", "infusion"),
  dose_type = NULL,
  pooled_ctrl = pooled_control(),
  ssctrl = ss_control()
)
```

## Arguments

- dat:

  A data frame containing pharmacokinetic data, including observed
  concentrations (DV), time after dose (tad), dose, and route
  information.

- half_life:

  The elimination half-life (t1/2) of the compound, used to identify
  early-phase Cmax values.

- single_point_base.lst:

  Optional list object returned by
  [run_single_point_base](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md)().
  If not supplied, the function will generate it internally.

- route:

  Route of administration. One of "bolus", "oral", or "infusion"
  (default = "bolus").

- dose_type:

  Optional string specifying the dosing type, passed to
  [run_single_point_base](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md)().

- pooled_ctrl:

  Control object created by
  [pooled_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md)(),
  defining data pooling options.

- ssctrl:

  Control object created by
  [ss_control](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md)(),
  defining steady-state control options.

## Value

A list containing individual and population Vd estimates and related
dose-level data.

## Details

Estimates individual apparent volumes of distribution from observed peak
concentrations. Individual estimates are then summarized to obtain a
population-level value.

For single-dose data, Vd is calculated according to the route of
administration:

- Bolus: \\V_d = \mathrm{Dose} / C\_{\mathrm{max}}\\

- Infusion: \\V_d = (\mathrm{Rate} \times t\_{\mathrm{inf}}) /
  C\_{\mathrm{max}}\\

- Oral: \\V_d = (\mathrm{Dose} \times F) / C\_{\mathrm{max}}\\, where
  \\F = 1 - e^{-k_a t}\\

For multiple-dose data, observed Cmax values are adjusted to single-dose
equivalents using the accumulation ratio: \$\$R\_{\mathrm{ac}} =
\frac{1}{1 - e^{-k_e \tau}}, \quad k_e = \ln(2)/t\_{1/2}\$\$ Adjusted
values are used to estimate Vd using the same route-specific equations.

## See also

[`run_single_point_base`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md),
[`trimmed_geom_mean`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trimmed_geom_mean.md),
[`pooled_control`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md),
[`ss_control`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md)

## Author

Zhonghui Huang

## Examples

``` r
# Process dataset
out <- processData(Bolus_1CPT)
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
# Get half-life and dose route
hf <- get_hf(dat = out$dat)$half_life_median
#> Estimating half-life....................
#> Half-life estimation complete: Estimated t1/2 = 11 h
rt <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
# Estimate Vd
approx.vc(dat = out$dat, half_life = hf, route = rt)$approx.vc.value
#> [1] 48.33771
```
