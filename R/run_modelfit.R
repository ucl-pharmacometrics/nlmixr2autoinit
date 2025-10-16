#' Run and evaluate one-compartment IV Model
#'
#' Perform parameter estimation with naive pooled data approach analysis for a one-compartment intravenous model, and evaluate the model fitting performance by the absolute prediction error (APE) and mean absolute prediction error (MAPE).
#'
#' @param dat A data frame containing intravenous pharmacokinetic data. It should include
#' columns such as `ID`, `EVID`, `DV`, `dose`, and `AMT` and other necessary variables.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is exp(1).
#' @param input.vd A numeric value for the initial estimate of volume of distribution (Vd). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is Default is 1.
#'
#' @details
#' This function fits the one-compartment IV model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`.
#' After fitting the model, the function returns the estimated CL and Vd along
#' with metrics to assess model fit.
#'
#' @return A list containing the following elements:
#' \item{npd.1cmpt_results}{A data frame with the estimated CL and Vd, as well
#' as the time spent during estimation.}
#' \item{npd.1cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.1cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.1cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.1cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.1cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.cl=exp(1),
                             input.vd=exp(1),
                             input.add=1) {
  start.time <- Sys.time()
  npd.list <-NA
  npd_results<-data.frame(cl = NA,
                          vd = NA,
                          timespent = NA)

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV  )[4],1)
  npd.rRMSE  <- round( metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV  )[5],1)

  return(
    list(
      npd.1cmpt_results = npd_results,
      npd.1cmpt.APE = npd.APE,
      npd.1cmpt.MAE = npd.MAE,
      npd.1cmpt.MAPE = npd.MAPE,
      npd.1cmpt.RMSE = npd.RMSE,
      npd.1cmpt.rRMSE = npd.rRMSE,
      nnpd.1cmpt.list = npd.list
    )
  )
}


#' Run and evaluate one-compartment IV model with Michaelis-Menten Kinetics

#' Perform parameter estimation on naive pooled data for a one-compartment Michaelis-Menten IV model.
#'
#' Estimates the pharmacokinetic parameters (Vmax, Km, Vd) for a
#' 1-compartment Michaelis-Menten IV (Intravenous) model by fitting.
#' Optionally, an initial estimate of Km is set based on the
#' maximum concentration in the dataset (if `km_threshold = TRUE`).
#'
#' @param dat A data frame containing pharmacokinetic data. Must include columns
#' such as `ID`, `EVID`, `DV`, `dose`, and `AMT`.
#' @param est.method The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".
#' @param npdmm_inputvmax A numeric value for the initial estimate of Vmax. Default is exp(1).
#' @param npdmm_inputkm A numeric value for the initial estimate of Km. Default is exp(1).
#' @param npdmm_inputcl A numeric value for the clearance. Default is exp(1).
#' @param npdmm_inputvd A numeric value for the initial estimate of volume of distribution (Vd). Default is exp(1).
#' @param input.add A numeric value representing the additive residual error component
#' included in the model fitting. Default is 1.
#' @param km_threshold A logical value (`TRUE` or `FALSE`). If `TRUE`,
#' initial estimates for \eqn{V_{max}} and \eqn{K_m} will be set based on the
#' observed maximum concentration in the dataset and referenced clearance.
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
#'   \item{npd.1cmpt.mm..APE}{The absolute prediction error (APE).}
#'   \item{npd.1cmpt.mm..MAE}{The mean absolute error (MAE).}
#'   \item{npd.1cmpt.mm..MAPE}{The mean absolute percentage error (MAPE).}
#'   \item{npd.1cmpt.mm..RMSE}{The root mean square error (RMSE).}
#'   \item{npd.1cmpt.mm.rRMSE}{The relative root mean square error (rRMSE).}
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
                                npdmm_inputvmax=exp(1),
                                npdmm_inputkm=exp(1),
                                npdmm_inputcl=exp(1),
                                npdmm_inputvd=exp(1),
                                input.add=1,
                                km_threshold=F) {
  start.time <- Sys.time()
  estvmax<-npdmm_inputvmax
  estkm<-npdmm_inputkm

  npdmm.list <-NA
  npdmm_results<-data.frame(  vmax = NA,
                              km =  NA,
                              vd =  NA,
                              timespent = NA)

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

  # Initial estimates of Vmax and Km will be set based on threshold if km_threshold=T
  # Determine the maximum concentration
  if (km_threshold){
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
    estmaxkm <- mean.pop.cmax * 4 # if km>>4cmax, it nearly fall into the linear range
    estkm<-mean.pop.cmax # initial km starts from cmax
    estvmax <-  estmaxkm * npdmm_inputcl
  }

  npdmm.list <- Fit_1cmpt_mm_iv(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.vmax =  estvmax,
    input.km = estkm,
    input.add = input.add,
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

  npd.APE <-  round(metrics.(pred.x = npdmm.list$cp, obs.y = npdmm.list$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x = npdmm.list$cp, obs.y = npdmm.list$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x = npdmm.list$cp, obs.y = npdmm.list$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =npdmm.list$cp, obs.y = npdmm.list$DV  )[5],1)


  return(
    list(
      npd.1cmpt.mm_results = npdmm_results,
      npd.1cmpt.mm.APE = npd.APE,
      npd.1cmpt.mm.MAE = npd.MAE,
      npd.1cmpt.mm.MAPE = npd.MAPE,
      npd.1cmpt.mm.RMSE = npd.RMSE,
      npd.1cmpt.mm.rRMSE = npd.rRMSE,
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
#' @param input.cl A numeric value for the initial estimate of the clearance (CL). Default is exp(1).
#' @param input.vc2cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is exp(1).
#' @param input.vp2cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is exp(1).
#' @param input.q2cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the two-compartment IV model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`.
#' After fitting the model, the function returns the estimated CL and Vd along
#' with metrics to assess model fit.
#'
#' @return A list containing the following elements:
#' \item{npd.2cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.2cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.2cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.2cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.2cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.2cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.cl=exp(1),
                             input.vc2cmpt=exp(1),
                             input.vp2cmpt=exp(1),
                             input.q2cmpt=exp(1),
                             input.add=1) {
  start.time <- Sys.time()

  npd.list.2cmpt <-NA
  npd_results_2cmpt <-
    data.frame(
      cl = NA,
      vc = NA,
      vp = NA,
      q = NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[5],1)

  return(
    list(
      npd.2cmpt_results = npd_results_2cmpt,
      npd.2cmpt.APE= npd.APE,
      npd.2cmpt.MAE  = npd.MAE,
      npd.2cmpt.MAPE = npd.MAPE,
      npd.2cmpt.RMSE = npd.RMSE,
      npd.2cmpt.rRMSE  = npd.rRMSE,
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
#' @param input.cl A numeric value for the initial estimate of log-transformed clearance (CL). Default is exp(1).
#' @param input.vc3cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is exp(1).
#' @param input.vp3cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is exp(1).
#' @param input.q3cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is exp(1).
#' @param input.vp23cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp2). Default is exp(1).
#' @param input.q23cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q2). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the three-compartment IV model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`.
#' After fitting the model, the function returns the estimated CL and Vd along
#' with metrics to assess model fit.
#'
#' @return A list containing the following elements:
#' \item{npd.3cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.3cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.3cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.3cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.3cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.3cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.cl= exp(1),
                             input.vc3cmpt = exp(1),
                             input.vp3cmpt = exp(1),
                             input.vp23cmpt =exp(1),
                             input.q3cmpt =  exp(1),
                             input.q23cmpt = exp(1),
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.3cmpt<-NA
  npd_results_3cmpt <-
    data.frame(
      cl = NA,
      vc = NA,
      vp = NA,
      vp2 = NA,
      q = NA,
      q2 = NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[5],1)


  return(
    list(
      npd.3cmpt_results = npd_results_3cmpt,
      npd.3cmpt.APE= npd.APE,
      npd.3cmpt.MAE  = npd.MAE,
      npd.3cmpt.MAPE = npd.MAPE,
      npd.3cmpt.RMSE = npd.RMSE,
      npd.3cmpt.rRMSE  = npd.rRMSE,
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
#' @param input.ka A numeric value for the initial estimate of the absorption rate constant. Default is exp(1).
#' @param input.cl A numeric value for the initial estimate of the clearance (CL). Default is exp(1).
#' @param input.vd A numeric value for the initial estimate of the volume of distribution (Vd). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the one-compartment oral model to the given dataset using
#' the specified estimation methods. It excludes rows where `EVID == 2`.
#' After fitting the model, the function returns the estimated CL and Vd
#' along with the metrics to assess model fit.
#'
#' @return A list containing the following elements:
#' \item{npd.1cmpt_results}{A data frame with the estimated CL and Vd, as well
#' as the time spent during estimation.}
#' \item{npd.1cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.1cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.1cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.1cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.1cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.ka= exp(1),
                             input.cl= exp(1),
                             input.vd= exp(1),
                             input.add=1) {
  start.time <- Sys.time()
  npd.list <-NA
  npd_results <-
    data.frame(
      ka= NA,
      cl = NA,
      vd = NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV  )[5],1)

  return(
    list(
      npd.1cmpt_results = npd_results,

      npd.1cmpt.APE = npd.APE,
      npd.1cmpt.MAE = npd.MAE,
      npd.1cmpt.MAPE = npd.MAPE,
      npd.1cmpt.RMSE = npd.RMSE,
      npd.1cmpt.rRMSE = npd.rRMSE,

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
#' @param est.method The estimation method in nlmixr2 to use. The default value is "nls".
#' @param input.ka A numeric value for the initial estimate of the log-transformed absorption rate constant.
#' Default is 1.
#' @param input.vmax A numeric value for the initial estimate of Vmax (maximum rate of metabolism).
#' Default is \code{exp(1)}. This serves as the initial value for \eqn{V_{max}} in the
#' nonlinear Michaelis-Menten oral absorption model.
#' @param input.km A numeric value for the initial estimate of Km (Michaelis constant).
#' Default is \code{exp(1)}. This defines the concentration at which the rate of metabolism
#' is half of \eqn{V_{max}}.
#' @param input.cl A numeric value for the clearance (CL) parameter. Default is \code{exp(1)}.
#' This value is also used for Vmax and Km calculation when \code{km_threshold = TRUE}.
#' @param input.vd A numeric value for the initial estimate of the apparent volume of distribution (Vd).
#' Default is \code{exp(1)}.
#' @param input.add A numeric value for the additive error model. Default is 1.
#' Adjust this if observed data show a different level of additive residual variability.
#' @param km_threshold A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},
#' initial estimates for \eqn{V_{max}} and \eqn{K_m} are adjusted based on the
#' observed maximum concentration and clearance.
#'
#' @details
#' If `km_threshold = TRUE`, this function first calculates initial estimates for \eqn{V_{max}} and \eqn{K_m}
#' based on the observed maximum concentration in the dataset and referenced clearance. This approach ensures
#' that the value of \eqn{K_m} is adjusted to lie between the linear and nonlinear regimes, providing a more robust
#' starting estimate when there is uncertainty about whether the model follows linear or nonlinear kinetics.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{npd.1cmpt.mm_results}{A data frame with the estimated values of  \eqn{K_a},\eqn{V_{max}},
#'   \eqn{K_m}, and \eqn{V_d}, as well as the time taken for the estimation.}
#'   \item{npd.1cmpt.mm..APE}{The absolute prediction error (APE).}
#'   \item{npd.1cmpt.mm..MAE}{The mean absolute error (MAE).}
#'   \item{npd.1cmpt.mm..MAPE}{The mean absolute percentage error (MAPE).}
#'   \item{npd.1cmpt.mm..RMSE}{The root mean square error (RMSE).}
#'   \item{npd.1cmpt.mm.rRMSE}{The relative root mean square error (rRMSE).}
#'   \item{npd.1cmpt.mm.list}{The full output list from the `Fit_1cmpt_mm_oral` function.}
#' }
#'
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' \dontrun{
#' # Example 1 (Linear kinetics case):
#' result <- run_npd_1cmpt_mm_oral(dat = Oral_1CPT,
#'                               est.method="foce",
#'                               input.cl = 4,
#'                               input.vd = 70,
#'                               km_threshold = TRUE)
#' result
#'
#' # Example 2 ( nonlinear-kinetics case):
#' result2 <- run_npd_1cmpt_mm_oral(dat = Oral_1CPTMM,
#'                               input.cl = 4,
#'                               input.vd = 70,
#'                               km_threshold = TRUE)
#' result2
#' }
#'
#' @export
run_npd_1cmpt_mm_oral <- function(dat,
                                est.method="nls",
                                input.ka=exp(1),
                                input.vmax= exp(1),
                                input.km= exp(1),
                                input.cl= exp(1),
                                input.vd= exp(1),
                                input.add= 1,
                                km_threshold=F) {
  start.time <- Sys.time()
  estvmax<-input.vmax
  estkm<-input.km

  npdmm.list <- NA
  npdmm_results <-
    data.frame(
      ka = NA,
      vmax =  NA,
      km =  NA,
      vd =  NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

  # Initial estimates of Vmax and Km will be set based on threshold if km_threshold=T
  # Determine the maximum concentration
  if (km_threshold){
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
    estmaxkm <- mean.pop.cmax * 4 # if km>>4cmax, it nearly fall into the linear range
    estkm<-mean.pop.cmax # initial km starts from cmax
    estvmax <-  estmaxkm * input.cl
  }

  npdmm.list <- Fit_1cmpt_mm_oral(
    data = dat[dat$EVID != 2, ],
    est.method = est.method,
    input.ka = input.ka,
    input.vmax =  estvmax,
    input.km = estkm,
    input.add =  input.add,
    input.vd =   input.vd
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

  npd.APE <-  round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x =  npdmm.list$cp, obs.y =npdmm.list$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV  )[5],1)

  return(
    list(
      npd.1cmpt.mm_results = npdmm_results,
      npd.1cmpt.mm.APE = npd.APE,
      npd.1cmpt.mm.MAE = npd.MAE,
      npd.1cmpt.mm.MAPE = npd.MAPE,
      npd.1cmpt.mm.RMSE = npd.RMSE,
      npd.1cmpt.mm.rRMSE = npd.rRMSE,
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
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is exp(1).
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is exp(1).
#' @param input.vc2cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is exp(1).
#' @param input.vp2cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is exp(1).
#' @param input.q2cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the two-compartment oral model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`
#' After fitting the model, the function returns the estimated CL, Vc, Vp, and Q
#' along with metrics to assess model fit.
#'
#'
#' @return A list containing the following elements:
#' \item{npd.2cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.2cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.2cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.2cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.2cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.2cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.ka= exp(1),
                             input.cl= exp(1),
                             input.vc2cmpt= exp(1),
                             input.vp2cmpt= exp(1),
                             input.q2cmpt= exp(1),
                             input.add= 1) {
  start.time <- Sys.time()
  npd.list.2cmpt<-NA
  npd_results_2cmpt <-
    data.frame(
      ka = NA,
      cl = NA,
      vc = NA,
      vp = NA,
      q = NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x =   npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV  )[5],1)

  return(
    list(
      npd.2cmpt_results = npd_results_2cmpt,
      npd.2cmpt.APE = npd.APE,
      npd.2cmpt.MAE = npd.MAE,
      npd.2cmpt.MAPE = npd.MAPE,
      npd.2cmpt.RMSE = npd.RMSE,
      npd.2cmpt.rRMSE = npd.rRMSE,
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
#' @param input.ka A numeric value for the initial estimate of absorption rate constant. Default is exp(1).
#' @param input.cl A numeric value for the initial estimate of clearance (CL). Default is exp(1).
#' @param input.vc3cmpt A numeric value for the initial estimate of the volume of distribution in the central compartment (Vc). Default is exp(1).
#' @param input.vp3cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is exp(1).
#' @param input.q3cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is exp(1).
#' @param input.vp23cmpt A numeric value for the initial estimate of the volume of distribution in the peripheral compartment (Vp). Default is exp(1).
#' @param input.q23cmpt A numeric value for the initial estimate of the intercompartmental clearance (Q). Default is exp(1).
#' @param input.add A numeric value for the additive error model. Default is 1.
#'
#' @details
#' This function fits the three-compartment oral model to the given dataset using
#' the specified estimation method. It excludes rows where `EVID == 2`
#' After fitting the model, the function returns the estimated CL, Vc, Vp, Vp2, Q
#' and Q2 along with metrics to assess model fit.
#'
#'
#' @return A list containing the following elements:
#' \item{npd.3cmpt_results}{A data frame with the estimated CL, Vc, Vp, Q, and
#' the time spent during estimation.}
#' \item{npd.3cmpt.APE}{The absolute prediction error (APE).}
#' \item{npd.3cmpt.MAE}{The mean absolute error (MAE).}
#' \item{npd.3cmpt.MAPE}{The mean absolute percentage error (MAPE).}
#' \item{npd.3cmpt.RMSE}{The root mean square error (RMSE).}
#' \item{npd.3cmpt.rRMSE}{The relative root mean square error (rRMSE).}
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
                             input.ka= exp(1),
                             input.cl= exp(1),
                             input.vc3cmpt = exp(1),
                             input.vp3cmpt = exp(1),
                             input.vp23cmpt = exp(1),
                             input.q3cmpt =  exp(1),
                             input.q23cmpt =  exp(1),
                             input.add=1) {
  start.time <- Sys.time()
  npd.list.3cmpt <-NA
  npd_results_3cmpt <-
    data.frame(
      ka = NA,
      cl = NA,
      vc = NA,
      vp =  NA,
      vp2 = NA,
      q = NA,
      q2 = NA,
      timespent = NA
    )

  npd.APE<-NA
  npd.MAE <-NA
  npd.MAPE<-NA
  npd.RMSE <-NA
  npd.rRMSE<-NA

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

  npd.APE <-  round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV  )[1],1)
  npd.MAE <-  round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV  )[2],1)
  npd.MAPE <- round( metrics.(pred.x =   npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV  )[3],1)
  npd.RMSE <-  round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV  )[4],1)
  npd.rRMSE <- round( metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV  )[5],1)


  return(
    list(
      npd.3cmpt_results = npd_results_3cmpt,
      npd.3cmpt.APE = npd.APE,
      npd.3cmpt.MAE = npd.MAE,
      npd.3cmpt.MAPE = npd.MAPE,
      npd.3cmpt.RMSE = npd.RMSE,
      npd.3cmpt.rRMSE = npd.rRMSE,
      npd.list.3cmpt =  npd.list.3cmpt
    )
  )
}
