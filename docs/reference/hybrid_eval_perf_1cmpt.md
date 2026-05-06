# Generate Unique Mixture Parameter Grid (with Deduplication and NA Removal)

Constructs a grid of all combinations of ka, cl, and vd parameters from
different sources (e.g., simpcal, graph, NCA methods), and removes
combinations that are redundant based on relative tolerance. Only oral
routes consider ka value for deduplication. Any parameter value set that
includes NA is removed up front.

## Usage

``` r
hybrid_eval_perf_1cmpt(
  route = "bolus",
  dat,
  sp_out_ka,
  sp_out_cl,
  sp_out_vd,
  graph_out_ka,
  graph_out_cl,
  graph_out_vd,
  nca_fd_ka,
  nca_fd_cl,
  nca_fd_vd,
  nca_efd_ka,
  nca_efd_cl,
  nca_efd_vd,
  nca_all_ka,
  nca_all_cl,
  nca_all_vd,
  verbose = TRUE
)
```

## Arguments

- route:

  Route of administration. Must be one of bolus, oral, or infusion.

- dat:

  A data.frame containing PK data with columns such as EVID and DV.

- sp_out_ka:

  Numeric; ka estimated from adaptive single-point methods.

- sp_out_cl:

  Numeric; clearance estimated from adaptive single-point methods.

- sp_out_vd:

  Numeric; volume of distribution estimated from adaptive single-point
  methods.

- graph_out_ka:

  Numeric; ka estimated from naive pooled graphic methods.

- graph_out_cl:

  Numeric; clearance estimated from naive pooled graphic methods.

- graph_out_vd:

  Numeric; volume of distribution estimated from naive pooled graphic
  methods.

- nca_fd_ka:

  Numeric; ka estimated from naive pooled NCA using first-dose data.

- nca_fd_cl:

  Numeric; clearance estimated from naive pooled NCA using first-dose
  data.

- nca_fd_vd:

  Numeric; volume of distribution estimated from naive pooled NCA using
  first-dose data.

- nca_efd_ka:

  Numeric; ka estimated from naive pooled NCA using repeated-dose data.

- nca_efd_cl:

  Numeric; clearance estimated from naive pooled NCA using repeated-dose
  data.

- nca_efd_vd:

  Numeric; volume of distribution estimated from naive pooled NCA using
  repeated-dose data.

- nca_all_ka:

  Numeric; ka estimated from naive pooled NCA using combined first- and
  repeated-dose data.

- nca_all_cl:

  Numeric; clearance estimated from naive pooled NCA using combined
  first- and repeated-dose data.

- nca_all_vd:

  Numeric; volume of distribution estimated from naive pooled NCA using
  combined first- and repeated-dose data.

- verbose:

  Logical; if TRUE (default), displays a textual progress bar during
  model evaluation using the 'progressr' package. Set to FALSE to run
  silently without showing progress information, which is recommended
  for automated analyses or CRAN checks.

## Value

A `data.frame` of unique parameter combinations with source labels and
values.

## Examples

``` r
dat <- Bolus_1CPT
# Example parameter estimates from different methods
sp_out_ka <- 1.2; sp_out_cl <- 3.5; sp_out_vd <- 50
graph_out_ka <- 1.1; graph_out_cl <- 3.6; graph_out_vd <- 52
nca_fd_ka <- 1.3; nca_fd_cl <- 3.4; nca_fd_vd <- 49
nca_efd_ka <- NA;  nca_efd_cl <- NA;  nca_efd_vd <- NA
nca_all_ka <- 1.25; nca_all_cl <- 3.55; nca_all_vd <- 51
# Run hybrid evaluation (silent)
 hybrid_eval_perf_1cmpt(
  route = "oral",
  dat = dat,
  sp_out_ka = sp_out_ka, sp_out_cl = sp_out_cl, sp_out_vd = sp_out_vd,
  graph_out_ka = graph_out_ka, graph_out_cl = graph_out_cl, graph_out_vd = graph_out_vd,
  nca_fd_ka = nca_fd_ka, nca_fd_cl = nca_fd_cl, nca_fd_vd = nca_fd_vd,
  nca_efd_ka = nca_efd_ka, nca_efd_cl = nca_efd_cl, nca_efd_vd = nca_efd_vd,
  nca_all_ka = nca_all_ka, nca_all_cl = nca_all_cl, nca_all_vd = nca_all_vd,
  verbose = FALSE
)
#>   ka_source cl_source vd_source ka_value cl_value vd_value     APE     MAE
#> 1   simpcal   simpcal   simpcal     1.20     3.50       50 1434780 206.413
#> 2     graph     graph     graph     1.10     3.60       52 1424299 204.906
#> 3    nca_fd    nca_fd    nca_fd     1.30     3.40       49 1452438 208.954
#> 4   nca_all   nca_all   nca_all     1.25     3.55       51 1407249 202.453
#>     MAPE    RMSE rRMSE1 rRMSE2
#> 1 63.235 363.167 58.743 54.581
#> 2 62.148 367.619 59.463 54.864
#> 3 66.146 361.330 58.446 54.361
#> 4 61.619 358.910 58.055 54.249
```
