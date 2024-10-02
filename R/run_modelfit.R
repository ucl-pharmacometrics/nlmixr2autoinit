#' Run and evaluate one-compartment IV Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a one-compartment intravenous model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vd A numeric value for the initial estimate of volume of distribution (Vd). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the 1-compartment IV model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL and Vd along with the APE (absolute prediction error) and MAPE
#' (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.1cmpt_results}{A data frame with the estimated CL and Vd, as well
#' as the time spent during estimation.}
#' \item{npd.1cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.1cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{nnpd.1cmpt.list}{The full output list from the `Fit_1cmpt_iv` function.}
#'
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_1cmpt_iv(dat = Bolus_1CPT,input.cl=4,input.vd=70)
#' }
#'
#' @export

run_npd_1cmpt_iv <- function(dat,
                             est.method="nls",
                             input.cl=1,
                             input.vd=1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list <-  Fit_1cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vd = input.vd,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  npd_results <-
    data.frame(
      cl = signif(npd.list$parFixedDf$`Back-transformed`[1], 3),
      vd = signif(npd.list$parFixedDf$`Back-transformed`[2], 3),
      timespent = time.spent
    )

  npd.APE <- sum(abs(npd.list$IRES), na.rm = T)
  npd.MAPE <- sum(abs(npd.list$IRES) / npd.list$DV) / nrow(npd.list) * 100

  # add penalty of rse
  if (max(npd.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.APE <- Inf
    npd.MAPE <- Inf
  }

  # Filter out failed runs
  if (npd.list$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.APE <- Inf
    npd.MAPE <- Inf
  }

  return(
    list(
      npd.1cmpt_results = npd_results,
      npd.1cmpt.APE = npd.APE,
      npd.1cmpt.MAPE = npd.MAPE,
      nnpd.1cmpt.list = npd.list
    )
  )
}


#' Run and evaluate one-compartment IV model with Michaelis-Menten Kinetics

#' Perform parameter estimation on naive pooled data for a one-compartment Michaelis-Menten IV model.
#'
#' Estimates the pharmacokinetic parameters (Vmax, Km, Vd) for a
#' 1-compartment Michaelis-Menten IV (Intravenous) model by fitting.
#' Optionally, an initial estimate of Km is based on the
#' maximum concentration in the dataset (if `km_threshold = TRUE`).
#'
#' @param dat A data frame containing pharmacokinetic data. Must include columns
#' such as `ID`, `EVID`, `DV`, `dose`, and `AMT`.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param npdmm_inputvmax A numeric value for the initial estimate of Vmax.
#' @param npdmm_inputkm A numeric value for the initial estimate of Km.
#' @param npdmm_inputcl A numeric value for the clearance.
#' @param npdmm_inputvd A numeric value for the initial estimate of volume of distribution (Vd).
#' @param km_threshold A logical value (`TRUE` or `FALSE`). If `TRUE`,
#' initial estimates for \eqn{V_{max}} and \eqn{K_m} will be set based on the
#' observed maximum concentration in the dataset and referenced clearance.
#'
#' @details
#' If `km_threshold = TRUE`, this function first calculates initial estimates for \eqn{V_{max}} and \eqn{K_m}
#' based on the observed maximum concentration in the dataset and referenced clearance. This approach ensures
#' that the value of \eqn{K_m} is adjusted to lie between the linear and nonlinear regimes, providing a more robust
#' starting estimate when there is uncertainty about whether the model follows linear or nonlinear kinetics.
#' This adjustment ensures that the fitting process does not result in \eqn{K_m} values too far from reality,
#' regardless of whether the dynamics are linear or nonlinear.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{npd.1cmpt.mm_results}{A data frame with the estimated values of \eqn{V_{max}},
#'   \eqn{K_m}, and \eqn{V_d}, as well as the time taken for the estimation.}
#'   \item{npd.1cmpt.mm.APE}{The absolute prediction error (APE).}
#'   \item{npd.1cmpt.mm.MAPE}{The mean absolute percentage error (MAPE).}
#'   \item{npd.1cmpt.mm.list}{The full output list from the `Fit_1cmpt_mm_iv` function.}
#' }
#'
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_1cmpt_mm_iv(dat = Bolus_1CPT,
#'                               est.method="foce",
#'                               npdmm_inputcl = 4,
#'                               npdmm_inputvd = 70,
#'                               km_threshold = TRUE)
#' }
#'
#' @export
run_npd_1cmpt_mm_iv <- function(dat,
                                est.method="nls",
                                npdmm_inputvmax=1,
                                npdmm_inputkm=1,
                                npdmm_inputcl=1,
                                npdmm_inputvd=1,
                                km_threshold=T) {
  start.time <- Sys.time()
  estvmax<-npdmm_inputvmax
  estkm<-npdmm_inputkm

  # Initial estimates of Vmax and Km will be set based on threshold if km_threshold=T
  # Determine the maximum concentration
  if (km_threshold){
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
    # max.dose <-
    #   max(dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0,]$dose)[1]
    # calc.cmax <- max.dose / npdmm_inputvd
    # fcmax <- max(c(mean.pop.cmax, calc.cmax))
    estmaxkm <- mean.pop.cmax * 4 # if km>>4cmax, it nearly fall into the linear range
    estkm<-mean.pop.cmax # initial km starts from cmax
    estvmax <-  estmaxkm * npdmm_inputcl
  }

  npdmm.list <- Fit_1cmpt_mm_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.vmax =  estvmax,
    input.km = estkm,
    input.add = 1,
    input.vd =   npdmm_inputvd
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 2)
  npdmm_results <-
    data.frame(
      vmax = signif(npdmm.list$parFixedDf$`Back-transformed`[1], 3),
      km = signif(npdmm.list$parFixedDf$`Back-transformed`[2], 3),
      vd = signif(npdmm.list$parFixedDf$`Back-transformed`[3], 3),
      timespent = time.spent
    )

  npdmm.APE <- sum(abs(npdmm.list$IRES), na.rm = T)
  npdmm.MAPE <-
    sum(abs(npdmm.list$IRES) / npdmm.list$DV) / nrow(npdmm.list) * 100

  countna <- is.na(npdmm.list$IRES)
  if (sum(countna) > (nrow(npdmm.list) / 2)) {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }

  if (max(npdmm.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }
  if (npdmm.list$message != "relative convergence (4)" &
      est.method == "focei") {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }

  return(
    list(
      npd.1cmpt.mm_results = npdmm_results,
      npd.1cmpt.mm.APE = npdmm.APE,
      npd.1cmpt.mm.MAPE = npdmm.MAPE,
      npd.1cmpt.mm.list = npdmm.list
    )
  )
}

#' Run and evaluate two-compartment IV Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a two-compartment intravenous model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vc2cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is 1.
#' @param input.vp2cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q2cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the two-compartment IV model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL, Vc, Vp, and Q along with the APE (absolute prediction error)
#' and MAPE (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.2cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.2cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.2cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.list.2cmpt}{The full output list from the `Fit_2cmpt_iv` function.}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_2cmpt_iv(dat = Bolus_2CPT,
#'                            input.cl = 4,
#'                            input.vc2cmpt = 35,
#'                            input.vp2cmpt = 35,
#'                            input.q2cmpt = 4)
#' }
#'
#' @export
#

run_npd_2cmpt_iv <- function(dat,
                             est.method="nls",
                             input.cl=1,
                             input.vc2cmpt=1,
                             input.vp2cmpt=1,
                             input.q2cmpt=1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.2cmpt <-  Fit_2cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vc2cmpt = input.vc2cmpt,
    input.vp2cmpt = input.vp2cmpt,
    input.q2cmpt = input.q2cmpt,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)

  npd_results_2cmpt <-
    data.frame(
      cl = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[1], 3),
      vc = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[2], 3),
      vp = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[3], 3),
      q = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[4], 3),
      timespent = time.spent
    )

  npd.2cmpt.APE <- sum(abs(npd.list.2cmpt$IRES), na.rm = T)
  npd.2cmpt.MAPE  <-
    sum(abs(npd.list.2cmpt$IRES) / npd.list.2cmpt$DV) / nrow(npd.list.2cmpt) *
    100

  if (max(npd.list.2cmpt$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.2cmpt.APE <- Inf
    npd.2cmpt.MAPE <- Inf
  }
  if (npd.list.2cmpt$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.2cmpt.APE <- Inf
    npd.2cmpt.MAPE <- Inf
  }

  return(
    list(
      npd.2cmpt_results = npd_results_2cmpt,
      npd.2cmpt.APE = npd.2cmpt.APE,
      npd.2cmpt.MAPE = npd.2cmpt.MAPE,
      npd.list.2cmpt = npd.list.2cmpt
    )
  )
}


#' Run and evaluate three-compartment IV Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a three-compartment intravenous model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vc3cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is 1.
#' @param input.vp3cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q3cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.vp23cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q23cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the three-compartment IV model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL, Vc, Vp, Vp2, Q and Q2 along with the APE (absolute prediction error)
#' and MAPE (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.3cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.3cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.3cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.list.3cmpt}{The full output list from the `Fit_3cmpt_iv` function.}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_3cmpt_iv(dat = Bolus_2CPT,
#'                            input.cl = 4,
#'                            input.vc3cmpt = 10,
#'                            input.vp3cmpt = 10,
#'                            input.vp23cmpt = 10,
#'                            input.q3cmpt = 4,
#'                            input.q23cmpt = 4
#'                            )
#' }
#'
#' @export
#

run_npd_3cmpt_iv <- function(dat,
                             est.method="nls",
                             input.cl=1,
                             input.vc3cmpt = 1,
                             input.vp3cmpt = 1,
                             input.vp23cmpt =1,
                             input.q3cmpt =  1,
                             input.q23cmpt = 1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.3cmpt <-  Fit_3cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vc3cmpt = input.vc3cmpt,
    input.vp3cmpt = input.vp3cmpt,
    input.vp23cmpt =   input.vp23cmpt,
    input.q3cmpt =  input.q3cmpt,
    input.q23cmpt =  input.q23cmpt,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)

  npd_results_3cmpt <-
    data.frame(
      cl = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[1], 3),
      vc = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[2], 3),
      vp = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[3], 3),
      vp2 = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[4], 3),
      q = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[5], 3),
      q2 = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[6], 3),
      timespent = time.spent
    )

  npd.3cmpt.APE <- sum(abs(npd.list.3cmpt$IRES), na.rm = T)
  npd.3cmpt.MAPE  <-
    sum(abs(npd.list.3cmpt$IRES) / npd.list.3cmpt$DV) / nrow(npd.list.3cmpt) *
    100

  if (max(npd.list.3cmpt$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.3cmpt.APE <- Inf
    npd.3cmpt.MAPE  <- Inf
  }
  if (npd.list.3cmpt$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.3cmpt.APE <- Inf
    npd.3cmpt.MAPE <- Inf
  }

  return(
    list(
      npd.3cmpt_results = npd_results_3cmpt,
      npd.3cmpt.APE = npd.3cmpt.APE,
      npd.3cmpt.MAPE = npd.3cmpt.MAPE,
      npd.list.3cmpt =  npd.list.3cmpt
    )
  )
}


#' Run and evaluate one-compartment Oral Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a one-compartment oral model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is 1.
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vd A numeric value for the initial estimate of volume of distribution (Vd). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the 1-compartment IV model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL and Vd along with the APE (absolute prediction error) and MAPE
#' (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.1cmpt_results}{A data frame with the estimated CL and Vd, as well
#' as the time spent during estimation.}
#' \item{npd.1cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.1cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{nnpd.1cmpt.list}{The full output list from the `Fit_1cmpt_oral` function.}
#'
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_1cmpt_oral(dat = Oral_1CPT,input.ka=1,input.cl=4,input.vd=70)
#' }
#'
#' @export

run_npd_1cmpt_oral <- function(dat,
                             est.method="nls",
                             input.ka=1,
                             input.cl=1,
                             input.vd=1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list <-  Fit_1cmpt_oral(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.ka = input.ka,
    input.cl = input.cl,
    input.vd = input.vd,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  npd_results <-
    data.frame(
      ka= signif(npd.list$parFixedDf$`Back-transformed`[1], 3),
      cl = signif(npd.list$parFixedDf$`Back-transformed`[2], 3),
      vd = signif(npd.list$parFixedDf$`Back-transformed`[3], 3),
      timespent = time.spent
    )

  npd.APE <- sum(abs(npd.list$IRES), na.rm = T)
  npd.MAPE <- sum(abs(npd.list$IRES) / npd.list$DV) / nrow(npd.list) * 100

  # add penalty of rse
  if (max(npd.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.APE <- Inf
    npd.MAPE <- Inf
  }

  # Filter out failed runs
  if (npd.list$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.APE <- Inf
    npd.MAPE <- Inf
  }

  return(
    list(
      npd.1cmpt_results = npd_results,
      npd.1cmpt.APE = npd.APE,
      npd.1cmpt.MAPE = npd.MAPE,
      nnpd.1cmpt.list = npd.list
    )
  )
}


#' Run and evaluate one-compartment oral model with Michaelis-Menten Kinetics

#' Perform parameter estimation on naive pooled data for a one-compartment Michaelis-Menten oral model.
#'
#' Estimates the pharmacokinetic parameters (Vmax, Km, Vd) for a
#' 1-compartment Michaelis-Menten Oral (Intravenous) model by fitting.
#' Optionally, an initial estimate of Km is based on the
#' maximum concentration in the dataset (if `km_threshold = TRUE`).
#'
#' @param dat A data frame containing pharmacokinetic data. Must include columns
#' such as `ID`, `EVID`, `DV`, `dose`, and `AMT`.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is 1.
#' @param npdmm_inputvmax A numeric value for the initial estimate of Vmax.
#' @param npdmm_inputkm A numeric value for the initial estimate of Km.
#' @param npdmm_inputcl A numeric value for the clearance.
#' @param npdmm_inputvd A numeric value for the initial estimate of volume of distribution (Vd).
#' @param km_threshold A logical value (`TRUE` or `FALSE`). If `TRUE`,
#' initial estimates of Vmax and Km are adjusted based on a threshold related
#' to the maximum observed concentration in the data.
#'
#' @details
#' If `km_threshold = TRUE`, this function first calculates initial estimates for \eqn{V_{max}} and \eqn{K_m}
#' based on the observed maximum concentration in the dataset and referenced clearance. This approach ensures
#' that the value of \eqn{K_m} is adjusted to lie between the linear and nonlinear regimes, providing a more robust
#' starting estimate when there is uncertainty about whether the model follows linear or nonlinear kinetics.
#' This adjustment ensures that the fitting process does not result in \eqn{K_m} values too far from reality,
#' regardless of whether the dynamics are linear or nonlinear.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{npd.1cmpt.mm_results}{A data frame with the estimated values of  \eqn{K_a},\eqn{V_{max}},
#'   \eqn{K_m}, and \eqn{V_d}, as well as the time taken for the estimation.}
#'   \item{npd.1cmpt.mm.APE}{The absolute prediction error (APE).}
#'   \item{npd.1cmpt.mm.MAPE}{The mean absolute percentage error (MAPE).}
#'   \item{npd.1cmpt.mm.list}{The full output list from the `Fit_1cmpt_mm_oral` function.}
#' }
#'
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' \dontrun{
#' # Example 1 (Linear kinetics):
#' result <- run_npd_1cmpt_mm_oral(dat = Oral_1CPT,
#'                               est.method="foce",
#'                               npdmm_inputcl = 4,
#'                               npdmm_inputvd = 70,
#'                               km_threshold = TRUE)
#' result
#'
#' # Example 2 ( nonlinear kinetics):
#' result2 <- run_npd_1cmpt_mm_oral(dat = Oral_1CPTMM,
#'                               est.method="foce",
#'                               npdmm_inputcl = 4,
#'                               npdmm_inputvd = 70,
#'                               km_threshold = TRUE)
#' result2
#' }
#'
#' @export
run_npd_1cmpt_mm_oral <- function(dat,
                                est.method="nls",
                                npdmm_inputka=1,
                                npdmm_inputvmax=1,
                                npdmm_inputkm=1,
                                npdmm_inputcl=1,
                                npdmm_inputvd=1,
                                km_threshold=T) {
  start.time <- Sys.time()
  estvmax<-npdmm_inputvmax
  estkm<-npdmm_inputkm

  # Initial estimates of Vmax and Km will be set based on threshold if km_threshold=T
  # Determine the maximum concentration
  if (km_threshold){
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
    # max.dose <-
    #   max(dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0,]$dose)[1]
    # calc.cmax <- max.dose / npdmm_inputvd
    # fcmax <- max(c(mean.pop.cmax, calc.cmax))
    estmaxkm <- mean.pop.cmax * 4 # if km>>4cmax, it nearly fall into the linear range
    estkm<-fcmax # initial km starts from cmax
    estvmax <-  estmaxkm * npdmm_inputcl
  }

  npdmm.list <- Fit_1cmpt_mm_oral(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.ka = input.ka,
    input.vmax =  estvmax,
    input.km = estkm,
    input.add = 1,
    input.vd =   npdmm_inputvd
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 2)
  npdmm_results <-
    data.frame(
      ka = signif(npdmm.list$parFixedDf$`Back-transformed`[1], 3),
      vmax = signif(npdmm.list$parFixedDf$`Back-transformed`[2], 3),
      km = signif(npdmm.list$parFixedDf$`Back-transformed`[3], 3),
      vd = signif(npdmm.list$parFixedDf$`Back-transformed`[4], 3),
      timespent = time.spent
    )

  npdmm.APE <- sum(abs(npdmm.list$IRES), na.rm = T)
  npdmm.MAPE <-
    sum(abs(npdmm.list$IRES) / npdmm.list$DV) / nrow(npdmm.list) * 100

  countna <- is.na(npdmm.list$IRES)
  if (sum(countna) > (nrow(npdmm.list) / 2)) {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }

  if (max(npdmm.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }
  if (npdmm.list$message != "relative convergence (4)" &
      est.method == "focei") {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }

  return(
    list(
      npd.1cmpt.mm_results = npdmm_results,
      npd.1cmpt.mm.APE = npdmm.APE,
      npd.1cmpt.mm.MAPE = npdmm.MAPE,
      npd.1cmpt.mm.list = npdmm.list
    )
  )
}

#' Run and evaluate two-compartment oral Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a two-compartment oral model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is 1.
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vc2cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is 1.
#' @param input.vp2cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q2cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the two-compartment IV model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL, Vc, Vp, and Q along with the APE (absolute prediction error)
#' and MAPE (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.2cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.2cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.2cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.list.2cmpt}{The full output list from the `Fit_2cmpt_oral` function.}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_2cmpt_oral(dat = Oral_2CPT,
#'                            input.ka = 1,
#'                            input.cl = 4,
#'                            input.vc2cmpt = 35,
#'                            input.vp2cmpt = 35,
#'                            input.q2cmpt = 4)
#' }
#'
#' @export
#

run_npd_2cmpt_oral <- function(dat,
                             est.method="nls",
                             input.ka=1,
                             input.cl=1,
                             input.vc2cmpt=1,
                             input.vp2cmpt=1,
                             input.q2cmpt=1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.2cmpt <-  Fit_2cmpt_oral(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.ka = input.ka,
    input.cl = input.cl,
    input.vc2cmpt = input.vc2cmpt,
    input.vp2cmpt = input.vp2cmpt,
    input.q2cmpt = input.q2cmpt,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)

  npd_results_2cmpt <-
    data.frame(
      ka = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[1], 3),
      cl = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[2], 3),
      vc = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[3], 3),
      vp = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[4], 3),
      q = signif(npd.list.2cmpt$parFixedDf$`Back-transformed`[5], 3),
      timespent = time.spent
    )

  npd.2cmpt.APE <- sum(abs(npd.list.2cmpt$IRES), na.rm = T)
  npd.2cmpt.MAPE  <-
    sum(abs(npd.list.2cmpt$IRES) / npd.list.2cmpt$DV) / nrow(npd.list.2cmpt) *
    100

  if (max(npd.list.2cmpt$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.2cmpt.APE <- Inf
    npd.2cmpt.MAPE <- Inf
  }
  if (npd.list.2cmpt$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.2cmpt.APE <- Inf
    npd.2cmpt.MAPE <- Inf
  }

  return(
    list(
      npd.2cmpt_results = npd_results_2cmpt,
      npd.2cmpt.APE = npd.2cmpt.APE,
      npd.2cmpt.MAPE = npd.2cmpt.MAPE,
      npd.list.2cmpt = npd.list.2cmpt
    )
  )
}


#' Run and evaluate three-compartment oral model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a three-compartment oral model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is 1.
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is 1.
#' @param input.vc3cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is 1.
#' @param input.vp3cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q3cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.vp23cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is 1.
#' @param input.q23cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is 1.
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the three-compartment IV model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`, as these
#' represent non-observation events. After fitting the model, the function returns
#' the estimated CL, Vc, Vp, Vp2, Q and Q2 along with the APE (absolute prediction error)
#' and MAPE (mean absolute percentage error) to assess model fit.
#'
#' The function also applies penalties to the results if certain conditions are
#' met, such as high relative standard error (`%RSE > 50`) or if the fitting method
#' does not converge successfully for the `focei` method.
#'
#' @return A list containing the following elements:
#' \item{npd.3cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.3cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.3cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.list.3cmpt}{The full output list from the `Fit_3cmpt_oral` function.}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- run_npd_3cmpt_oral(dat = Oral_2CPT,
#'                            input.ka = 1,
#'                            input.cl = 4,
#'                            input.vc3cmpt = 10,
#'                            input.vp3cmpt = 10,
#'                            input.vp23cmpt = 10,
#'                            input.q3cmpt = 4,
#'                            input.q23cmpt = 4
#'                            )
#' }
#'
#' @export
#

run_npd_3cmpt_oral <- function(dat,
                             est.method="nls",
                             input.ka=1,
                             input.cl=1,
                             input.vc3cmpt = 1,
                             input.vp3cmpt = 1,
                             input.vp23cmpt =1,
                             input.q3cmpt =  1,
                             input.q23cmpt = 1,
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.3cmpt <-  Fit_3cmpt_oral(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.ka = input.ka,
    input.cl = input.cl,
    input.vc3cmpt = input.vc3cmpt,
    input.vp3cmpt = input.vp3cmpt,
    input.vp23cmpt =   input.vp23cmpt,
    input.q3cmpt =  input.q3cmpt,
    input.q23cmpt =  input.q23cmpt,
    input.add = input.add
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)

  npd_results_3cmpt <-
    data.frame(
      ka = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[1], 3),
      cl = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[2], 3),
      vc = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[3], 3),
      vp = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[4], 3),
      vp2 = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[5], 3),
      q = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[6], 3),
      q2 = signif(npd.list.3cmpt$parFixedDf$`Back-transformed`[7], 3),
      timespent = time.spent
    )

  npd.3cmpt.APE <- sum(abs(npd.list.3cmpt$IRES), na.rm = T)
  npd.3cmpt.MAPE  <-
    sum(abs(npd.list.3cmpt$IRES) / npd.list.3cmpt$DV) / nrow(npd.list.3cmpt) *
    100

  if (max(npd.list.3cmpt$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.3cmpt.APE <- Inf
    npd.3cmpt.MAPE  <- Inf
  }
  if (npd.list.3cmpt$message != "relative convergence (4)" &
      est.method == "focei") {
    npd.3cmpt.APE <- Inf
    npd.3cmpt.MAPE <- Inf
  }

  return(
    list(
      npd.3cmpt_results = npd_results_3cmpt,
      npd.3cmpt.APE = npd.3cmpt.APE,
      npd.3cmpt.MAPE = npd.3cmpt.MAPE,
      npd.list.3cmpt =  npd.list.3cmpt
    )
  )
}
