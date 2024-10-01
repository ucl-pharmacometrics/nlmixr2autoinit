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
                             est.method,
                             input.cl,
                             input.vd) {
  start.time <- Sys.time()
  npd.list <-  Fit_1cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vd = input.vd,
    input.add = 1
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

  if (max(npd.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npd.APE <- Inf
    npd.MAPE <- Inf
  }
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


#' Conduct and evaluate one-compartment with Michaelis-Menten Kinetics modelling
#'
#' Perform parameter estimation with naive pooled data approach analysis for a one-compartment intravenous model, and use provided or calculated Vmax, Km and volume of distribution estimates to fit the model and calculate the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param npdmm_inputcl Pre-estimated clearance for the intravenous data, used to determine the boundary values for Vmax and Km.
#' @param npdmm_inputvd The initial estimate of volume of distribution.
#' @return A list containing the results of the model fitting, including the fitted parameters (`npd.1cmpt.mm_results` and `npd.1cmpt.mm_results.2`), APE (`npd.1cmpt.mm.APE` and `npd.1cmpt.mm.APE.2`), MAPE (`npd.1cmpt.mm.MAPE` and `npd.1cmpt.mm.MAPE.2`), and the model fitting lists (`npd.1cmpt.mm.list` and `npd.1cmpt.mm.list2`).
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPTMM
#' run_npd_1cmpt_mm_iv(dat=dat,est.method="nls",input.cl=4,input.vd=70)
#' @export
#'
#'
run_npd_1cmpt_mm_iv <- function(dat,
                                est.method,
                                npdmm_inputcl,
                                npdmm_inputvd) {
  # Scenario 1 start with default 1
  start.time <- Sys.time()
  npdmm.list <- Fit_1cmpt_mm_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.vmax = 1,
    input.km = 1,
    input.add = 1,
    input.vd =  npdmm_inputvd
  )

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
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
  # add penalty of rse
  if (max(npdmm.list$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }
  if (npdmm.list$message != "relative convergence (4)" &
      est.method == "focei") {
    npdmm.APE <- Inf
    npdmm.MAPE <- Inf
  }

  # Scenario 2 start with threshold value
  # determine the maximum concentration
  dat.obs <- dat[dat$EVID == 0, ]
  pop.cmax <- aggregate(dat.obs$DV,
                        list(dat.obs$ID),
                        FUN = max,
                        na.rm = T)
  mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
  max.dose <-
    max(dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0,]$dose)[1]
  #calculated cmax (c0) based on volume of distribution
  calc.cmax <- max.dose / npdmm_inputvd
  fcmax <- max(c(mean.pop.cmax, calc.cmax))

  linear_minkm <-
    fcmax * 4 # if km>>4cmax, it nearly fall into the linear range
  linear_maxvmax <-  linear_minkm * npdmm_inputcl

  start.time <- Sys.time()
  npdmm.list2 <- Fit_1cmpt_mm_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.vmax =  linear_maxvmax,
    input.km = linear_minkm / 4,
    # start from 1 fold of cmax
    input.add = 1,
    input.vd =   npdmm_inputvd
  )


  end.time <- Sys.time()
  time.spent2 <- round(difftime(end.time, start.time), 2)
  npdmm_results2 <-
    data.frame(
      vmax = signif(npdmm.list2$parFixedDf$`Back-transformed`[1], 3),
      km = signif(npdmm.list2$parFixedDf$`Back-transformed`[2], 3),
      vd = signif(npdmm.list2$parFixedDf$`Back-transformed`[3], 3),
      timespent = time.spent2
    )

  npdmm2.APE <- sum(abs(npdmm.list2$IRES), na.rm = T)
  npdmm2.MAPE <-
    sum(abs(npdmm.list2$IRES) / npdmm.list2$DV) / nrow(npdmm.list2) * 100

  countna <- is.na(npdmm.list2$IRES)
  if (sum(countna) > (nrow(npdmm.list2) / 2)) {
    npdmm2.APE <- Inf
    npdmm2.MAPE <- Inf
  }

  if (max(npdmm.list2$parFixedDf$`%RSE`, na.rm = T) > 50) {
    npdmm2.APE <- Inf
    npdmm2.MAPE <- Inf
  }
  if (npdmm.list2$message != "relative convergence (4)" &
      est.method == "focei") {
    npdmm2.APE <- Inf
    npdmm2.MAPE <- Inf
  }


  return(
    list(
      npd.1cmpt.mm_results = npdmm_results,
      npd.1cmpt.mm.APE = npdmm.APE,
      npd.1cmpt.mm.MAPE = npdmm.MAPE,
      npd.1cmpt.mm_results.2 = npdmm_results2,
      npd.1cmpt.mm.APE.2 = npdmm2.APE,
      npd.1cmpt.mm.MAPE.2 = npdmm2.MAPE,
      npd.1cmpt.mm.list = npdmm.list,
      npd.1cmpt.mm.list2 = npdmm.list2
    )
  )
}



#' Conduct and evaluate two-compartment modelling
#'
#' Perform parameter estimation with naive pooled data approach analysis for a two-compartment intravenous model, and use provided clearance and volume of distribution estimates to fit the model and calculate the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl The initial estimate for clearance.
#' @param input.vd The initial estimate for volume of distribution.
#' @return A list containing the results of the model fitting, including the fitted parameters (`npd.2cmpt_results`), APE (`npd.2cmpt.APE`), MAPE (`npd.2cmpt.MAPE`), and the model fitting list (`npd.list.2cmpt`).
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' run_npd_2cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vd=70)
#' @export

run_npd_2cmpt_iv <- function(dat,
                             est.method,
                             input.cl,
                             input.vd) {
  start.time <- Sys.time()
  npd.list.2cmpt <-  Fit_2cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vc2cmpt = input.vd / 2,
    input.vp2cmpt = input.vd / 2,
    input.q2cmpt = input.cl,
    input.add = 1
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


#' Conduct and evaluate three-compartment modelling
#'
#' Perform parameter estimation with naive pooled data approach analysis for a three-compartment intravenous model, and use provided clearance and volume of distribution estimates to fit the model and calculate the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl The initial estimate for clearance.
#' @param input.vd The initial estimate for volume of distribution.
#' @return A list containing the results of the model fitting, including the fitted parameters (`npd.3cmpt_results`), APE (`npd.3cmpt.APE`), MAPE (`npd.3cmpt.MAPE`), and the model fitting list (`npd.list.3cmpt`).
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' dat <- Bolus_2CPT
#' run_npd_3cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vd=70)
#' @export

run_npd_3cmpt_iv <- function(dat,
                             est.method,
                             input.cl,
                             input.vd) {
  start.time <- Sys.time()
  npd.list.3cmpt <-  Fit_3cmpt_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.cl = input.cl,
    input.vc3cmpt = input.vd / 3,
    input.vp3cmpt = input.vd / 3,
    input.vp23cmpt =  input.vd / 3,
    input.q3cmpt = input.cl,
    input.q23cmpt = input.cl,
    input.add = 1
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
