# Estimate individual-level residual error from the elimination phase

Performs log-linear regression on the elimination phase of a single
individual's or one group's pharmacokinetic concentration–time data to
estimate additive and proportional residual standard deviations.

## Usage

``` r
getsigmas(group_df, nlastpoints = 3)
```

## Arguments

- group_df:

  A data frame for a single group (e.g., one subject or dose),
  containing columns: EVID (event ID), DV (observed concentration), TIME
  (time after dose), and routeobs (administration route).

- nlastpoints:

  Integer specifying the number of terminal data points used for
  regression.

## Value

A tibble with the following columns:

- intercept: Intercept of the log-linear regression line

- slope: Estimate of the terminal elimination rate constant

- residual_sd_additive: Standard deviation of additive residuals

- residual_sd_proportional: Standard deviation of proportional residuals

## Details

Residuals are computed from individual-predicted concentrations (IPRED)
and observed concentrations (DV) using the following definitions: \$\$
\sigma\_{add} = \sqrt{Var(C\_{obs} - C\_{pred})} \$\$

\$\$ \sigma\_{prop} = \sqrt{Var\left(\frac{C\_{obs}}{C\_{pred}} -
1\right)} \$\$

where \\C\_{obs}\\ is the observed concentration and \\C\_{pred}\\ is
the model-predicted concentration obtained by back-transformation of the
log-linear regression. The additive residual standard deviation
(\\\sigma\_{add}\\) and proportional residual standard deviation
(\\\sigma\_{prop}\\) are calculated per individual.

## Author

Zhonghui Huang

## Examples

``` r
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
getsigmas(dat[dat$ID == 1 & dat$dose_number == 1 & dat$resetflag == 1 &
              dat$EVID == 0, ])
#> # A tibble: 1 × 4
#>   intercept   slope residual_sd_additive residual_sd_proportional
#>       <dbl>   <dbl>                <dbl>                    <dbl>
#> 1      6.97 -0.0679                 1.70                   0.0750
```
