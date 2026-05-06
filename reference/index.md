# Package index

## Main Function

Primary function and control options for automated initial parameter
estimation

- [`getPPKinits()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getPPKinits.md)
  : Automated pipeline for generating initial estimates in population PK
  models

- [`print(`*`<getPPKinits>`*`)`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/print.getPPKinits.md)
  :

  Print method for `getPPKinits` objects

- [`initsControl()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/initsControl.md)
  : Create full control list for initial parameter estimation

- [`fallback_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/fallback_control.md)
  : Control settings for fallback rules in parameter estimation

- [`nca_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/nca_control.md)
  : Control options for non-compartmental analysis

- [`pooled_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md)
  : Control settings for pooled data analysis

- [`ss_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md)
  : Internal control builder for steady-state evaluation

## Graphic and NCA Methods

Naive pooled graphic methods and noncompartmental analysis

- [`graphcal_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_iv.md)
  : Graphical calculation of clearance and volume of distribution (IV
  route)
- [`graphcal_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_oral.md)
  : Graphical calculation of pharmacokinetic parameters for oral
  administration
- [`run_graphcal()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_graphcal.md)
  : Run graphical analysis of pharmacokinetic parameters
- [`getnca()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getnca.md)
  : Perform non-compartmental pharmacokinetic analysis
- [`run_pooled_nca()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_pooled_nca.md)
  : Performs non-compartmental analysis on pooled data
- [`get_pooled_data()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)
  : Generate pooled data for pharmacokinetic analysis

## Single-Point Methods

Adaptive single-point estimation functions

- [`run_single_point()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point.md)
  : Run full adaptive single-point PK analysis
- [`run_single_point_base()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md)
  : Run adaptive single-point pharmacokinetic analysis
- [`run_single_point_extra()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_extra.md)
  : Perform extended single-point pharmacokinetic calculations
- [`run_ka_solution()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_ka_solution.md)
  : Estimate the absorption rate constant using pointwise methods
- [`ka_calculation_md()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_md.md)
  : Calculate absorption rate constant (ka) in a multiple-dose
  one-compartment model
- [`ka_calculation_sd()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_sd.md)
  : Estimate absorption rate constant in a one-compartment oral model
- [`ka_wanger_nelson()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_wanger_nelson.md)
  : Calculate the absorption rate constant using the Wagner-Nelson
  method

## Parameter Sweeping

Model fitting and parameter sweeping across PK structures

- [`Fit_1cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a one-compartment linear
  elimination model
- [`Fit_1cmpt_mm_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_iv.md)
  : Fit intravenous pharmacokinetic data to a one-compartment model with
  Michaelis-Menten elimination
- [`Fit_1cmpt_mm_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_oral.md)
  : Fit oral pharmacokinetic data to a one-compartment model with
  Michaelis-Menten elimination
- [`Fit_1cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_oral.md)
  : Fit oral pharmacokinetic data to a one-compartment linear
  elimination model
- [`Fit_2cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a two-compartment linear
  elimination model
- [`Fit_2cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_oral.md)
  : Fit oral pharmacokinetic data to a two-compartment model
- [`Fit_3cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a three-compartment linear
  elimination model
- [`Fit_3cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_oral.md)
  : Fit oral pharmacokinetic data to a three-compartment linear
  elimination model
- [`run_npd_1cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_iv.md)
  : Run and evaluate a one-compartment IV model
- [`run_npd_1cmpt_mm_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_mm_iv.md)
  : Run and evaluate a one-compartment IV Michaelis-Menten model
- [`run_npd_1cmpt_mm_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_mm_oral.md)
  : Run and evaluate a one-compartment oral model with Michaelis-Menten
  kinetics
- [`run_npd_1cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_oral.md)
  : Run and evaluate a one-compartment oral model
- [`run_npd_2cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_2cmpt_iv.md)
  : Run and evaluate a two-compartment IV model
- [`run_npd_2cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_2cmpt_oral.md)
  : Run and evaluate a two-compartment oral model
- [`run_npd_3cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_3cmpt_iv.md)
  : Run and evaluate a three-compartment IV model
- [`run_npd_3cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_3cmpt_oral.md)
  : Run and evaluate a three-compartment oral model
- [`sim_sens_1cmpt_mm()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_1cmpt_mm.md)
  : Parameter sweeping for a one-compartment Michaelis-Menten model
- [`sim_sens_2cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_2cmpt.md)
  : Parameter sweeping for a two-compartment pharmacokinetic model
- [`sim_sens_3cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_3cmpt.md)
  : Parameter sweeping for a three-compartment pharmacokinetic model

## Variability Estimation

Residual and inter-individual variability estimation

- [`getsigma()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getsigma.md)
  : Compute overall residual variability from elimination phase
- [`getsigmas()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getsigmas.md)
  : Estimate individual-level residual error from the elimination phase
- [`getOmegas()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getOmegas.md)
  : Generate ETA variance and covariance table

## All Other Functions

Data processing and supporting utility functions

- [`Fit_1cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a one-compartment linear
  elimination model

- [`Fit_1cmpt_mm_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_iv.md)
  : Fit intravenous pharmacokinetic data to a one-compartment model with
  Michaelis-Menten elimination

- [`Fit_1cmpt_mm_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_mm_oral.md)
  : Fit oral pharmacokinetic data to a one-compartment model with
  Michaelis-Menten elimination

- [`Fit_1cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_1cmpt_oral.md)
  : Fit oral pharmacokinetic data to a one-compartment linear
  elimination model

- [`Fit_2cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a two-compartment linear
  elimination model

- [`Fit_2cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_2cmpt_oral.md)
  : Fit oral pharmacokinetic data to a two-compartment model

- [`Fit_3cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_iv.md)
  : Fit intravenous pharmacokinetic data to a three-compartment linear
  elimination model

- [`Fit_3cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/Fit_3cmpt_oral.md)
  : Fit oral pharmacokinetic data to a three-compartment linear
  elimination model

- [`approx.vc()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/approx.vc.md)
  : Approximate volume of distribution from observed Cmax

- [`bin.time()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/bin.time.md)
  : Bin time-concentration data using quantile or algorithmic binning

- [`calculate_cl()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/calculate_cl.md)
  : Calculate clearance using an adaptive single-point method

- [`calculate_tad()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/calculate_tad.md)
  : Calculate time after dose for pharmacokinetic data

- [`calculate_vd()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/calculate_vd.md)
  : Calculates volume of distribution from concentration data

- [`eval_perf_1cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/eval_perf_1cmpt.md)
  : Evaluates predictive performance of a one-compartment model

- [`fallback_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/fallback_control.md)
  : Control settings for fallback rules in parameter estimation

- [`find_best_lambdaz()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/find_best_lambdaz.md)
  : Find the best terminal elimination rate constant (lambdaz)

- [`force_find_lambdaz()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/force_find_lambdaz.md)
  : Forceful estimation of terminal slope

- [`getOmegas()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getOmegas.md)
  : Generate ETA variance and covariance table

- [`getPPKinits()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getPPKinits.md)
  : Automated pipeline for generating initial estimates in population PK
  models

- [`get_hf()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_hf.md)
  : Estimate half-life from pooled pharmacokinetic data

- [`get_pooled_data()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/get_pooled_data.md)
  : Generate pooled data for pharmacokinetic analysis

- [`getnca()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getnca.md)
  : Perform non-compartmental pharmacokinetic analysis

- [`getsigma()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getsigma.md)
  : Compute overall residual variability from elimination phase

- [`getsigmas()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/getsigmas.md)
  : Estimate individual-level residual error from the elimination phase

- [`graphcal_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_iv.md)
  : Graphical calculation of clearance and volume of distribution (IV
  route)

- [`graphcal_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/graphcal_oral.md)
  : Graphical calculation of pharmacokinetic parameters for oral
  administration

- [`hybrid_eval_perf_1cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/hybrid_eval_perf_1cmpt.md)
  : Generate Unique Mixture Parameter Grid (with Deduplication and NA
  Removal)

- [`initsControl()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/initsControl.md)
  : Create full control list for initial parameter estimation

- [`is_ss()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/is_ss.md)
  : Determine steady state for pharmacokinetic observations

- [`ka_calculation_md()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_md.md)
  : Calculate absorption rate constant (ka) in a multiple-dose
  one-compartment model

- [`ka_calculation_sd()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_calculation_sd.md)
  : Estimate absorption rate constant in a one-compartment oral model

- [`ka_wanger_nelson()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ka_wanger_nelson.md)
  : Calculate the absorption rate constant using the Wagner-Nelson
  method

- [`mark_dose_number()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/mark_dose_number.md)
  : Mark dose number

- [`metrics.()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/metrics..md)
  : Calculate metrics for model predictive performance evaluation

- [`nca_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/nca_control.md)
  : Control options for non-compartmental analysis

- [`nmpkconvert()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/nmpkconvert.md)
  : Expand additional dosing (ADDL) records for pharmacokinetic analysis

- [`pooled_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/pooled_control.md)
  : Control settings for pooled data analysis

- [`print(`*`<getPPKinits>`*`)`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/print.getPPKinits.md)
  :

  Print method for `getPPKinits` objects

- [`processData()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/processData.md)
  : Process time–concentration dataset for pharmacokinetic analysis

- [`run_graphcal()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_graphcal.md)
  : Run graphical analysis of pharmacokinetic parameters

- [`run_ka_solution()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_ka_solution.md)
  : Estimate the absorption rate constant using pointwise methods

- [`run_npd_1cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_iv.md)
  : Run and evaluate a one-compartment IV model

- [`run_npd_1cmpt_mm_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_mm_iv.md)
  : Run and evaluate a one-compartment IV Michaelis-Menten model

- [`run_npd_1cmpt_mm_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_mm_oral.md)
  : Run and evaluate a one-compartment oral model with Michaelis-Menten
  kinetics

- [`run_npd_1cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_1cmpt_oral.md)
  : Run and evaluate a one-compartment oral model

- [`run_npd_2cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_2cmpt_iv.md)
  : Run and evaluate a two-compartment IV model

- [`run_npd_2cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_2cmpt_oral.md)
  : Run and evaluate a two-compartment oral model

- [`run_npd_3cmpt_iv()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_3cmpt_iv.md)
  : Run and evaluate a three-compartment IV model

- [`run_npd_3cmpt_oral()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_npd_3cmpt_oral.md)
  : Run and evaluate a three-compartment oral model

- [`run_pooled_nca()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_pooled_nca.md)
  : Performs non-compartmental analysis on pooled data

- [`run_single_point()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point.md)
  : Run full adaptive single-point PK analysis

- [`run_single_point_base()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_base.md)
  : Run adaptive single-point pharmacokinetic analysis

- [`run_single_point_extra()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/run_single_point_extra.md)
  : Perform extended single-point pharmacokinetic calculations

- [`sim_sens_1cmpt_mm()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_1cmpt_mm.md)
  : Parameter sweeping for a one-compartment Michaelis-Menten model

- [`sim_sens_2cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_2cmpt.md)
  : Parameter sweeping for a two-compartment pharmacokinetic model

- [`sim_sens_3cmpt()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/sim_sens_3cmpt.md)
  : Parameter sweeping for a three-compartment pharmacokinetic model

- [`ss_control()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/ss_control.md)
  : Internal control builder for steady-state evaluation

- [`trapezoidal_linear()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trapezoidal_linear.md)
  : Linear trapezoidal rule

- [`trapezoidal_linear_up_log_down()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trapezoidal_linear_up_log_down.md)
  : Linear-up and log-down trapezoidal rule

- [`trimmed_geom_mean()`](https://ucl-pharmacometrics.github.io/nlmixr2autoinit/reference/trimmed_geom_mean.md)
  : Computes the trimmed geometric mean
