# Generate ETA variance and covariance table

This function constructs a combined table containing:

- **ETA variances** (e.g., `eta.cl`, `eta.vc`), which represent
  inter-individual variability (IIV) in pharmacokinetic parameters;

- **Derived covariances** (e.g., `cov.eta_cl_vc`) computed from ETA
  variances and assumed pairwise correlations.

## Usage

``` r
getOmegas()
```

## Value

A `data.frame` with columns: `Parameters`, `Methods`, and `Values`.

## Details

ETA variances are initialized to 0.1 by default. Correlations within
defined omega blocks (block 1: `eta.vmax`, `eta.km`; block 2: `eta.cl`,
`eta.vc`, `eta.vp`, `eta.vp2`, `eta.q`, `eta.q2`) are assumed to be 0.1
and used to compute covariances as: \$\$Cov(i, j) = sqrt(Var_i) \*
sqrt(Var_j) \* Corr(i, j)\$\$

The resulting output format aligns with `Recommended_initial_estimates`.

## Examples

``` r
getOmegas()
#>         Parameters                  Methods Values
#> 1           eta.ka             fixed_values    0.1
#> 2           eta.cl             fixed_values    0.1
#> 3           eta.vc             fixed_values    0.1
#> 4           eta.vp             fixed_values    0.1
#> 5            eta.q             fixed_values    0.1
#> 6          eta.vp2             fixed_values    0.1
#> 7           eta.q2             fixed_values    0.1
#> 8         eta.vmax             fixed_values    0.1
#> 9           eta.km             fixed_values    0.1
#> 10 cor.eta_vmax_km eta_corr_derived (r=0.1)   0.01
#> 11   cor.eta_cl_vc eta_corr_derived (r=0.1)   0.01
#> 12   cor.eta_cl_vp eta_corr_derived (r=0.1)   0.01
#> 13  cor.eta_cl_vp2 eta_corr_derived (r=0.1)   0.01
#> 14    cor.eta_cl_q eta_corr_derived (r=0.1)   0.01
#> 15   cor.eta_cl_q2 eta_corr_derived (r=0.1)   0.01
#> 16   cor.eta_vc_vp eta_corr_derived (r=0.1)   0.01
#> 17  cor.eta_vc_vp2 eta_corr_derived (r=0.1)   0.01
#> 18    cor.eta_vc_q eta_corr_derived (r=0.1)   0.01
#> 19   cor.eta_vc_q2 eta_corr_derived (r=0.1)   0.01
#> 20  cor.eta_vp_vp2 eta_corr_derived (r=0.1)   0.01
#> 21    cor.eta_vp_q eta_corr_derived (r=0.1)   0.01
#> 22   cor.eta_vp_q2 eta_corr_derived (r=0.1)   0.01
#> 23   cor.eta_vp2_q eta_corr_derived (r=0.1)   0.01
#> 24  cor.eta_vp2_q2 eta_corr_derived (r=0.1)   0.01
#> 25    cor.eta_q_q2 eta_corr_derived (r=0.1)   0.01
```
