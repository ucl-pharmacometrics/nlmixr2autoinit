#' Run and evaluate a one-compartment IV model
#'
#' Fits a one-compartment intravenous pharmacokinetic model using a naive pooled
#' data approach and evaluates model performance based on prediction error
#' metrics.
#'
#' @param dat A data frame containing raw time–concentration data in standard
#'   nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param input.vd Initial estimate for volume of distribution (Vd). Defaults to
#'   exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded prior to model fitting. The model is fitted
#' using `Fit_1cmpt_iv`, and prediction-based metrics are calculated to assess
#' performance.
#'
#' @return A list containing the fitted parameter estimates and prediction error
#' metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_1cmpt_iv(dat = Bolus_1CPT, input.cl = 4, input.vd = 70)
#' }
#'
#' @seealso \code{\link{Fit_1cmpt_iv}}
#' @export

run_npd_1cmpt_iv <- function(dat,
                             est.method = "nls",
                             input.cl = exp(1),
                             input.vd = exp(1),
                             input.add = 1) {
  start.time <- Sys.time()
  npd.list <- NA
  npd_results <- data.frame(cl = NA,
                            vd = NA,
                            timespent = NA)

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list <-  Fit_1cmpt_iv(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV)[1], 1)
  npd.MAE <-
    round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV)[2], 1)
  npd.MAPE <-
    round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV)[3], 1)
  npd.RMSE <-
    round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV)[4], 1)
  npd.rRMSE  <-
    round(metrics.(pred.x = npd.list$cp, obs.y = npd.list$DV)[5], 1)

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


#' Run and evaluate a one-compartment IV Michaelis-Menten model
#'
#' Fits a one-compartment intravenous pharmacokinetic model with
#' Michaelis-Menten elimination using a naive pooled data approach and evaluates
#' model performance based on prediction error metrics.
#'
#' @param dat A data frame containing pharmacokinetic data in standard nlmixr2
#'   format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param npdmm_inputvmax Initial estimate for Vmax. Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param npdmm_inputkm Initial estimate for Km. Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param npdmm_inputcl Initial estimate for clearance (CL). Defaults to exp(1)
#'  , corresponding to a log-scale value of 1.
#' @param npdmm_inputvd Initial estimate for volume of distribution (Vd).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#' @param km_threshold Logical value. If TRUE, initial estimates for Vmax and Km
#'   are calculated based on the maximum observed concentration.
#'
#' @details
#' Rows where `EVID == 2` are excluded before model fitting. The model is
#' fitted using `Fit_1cmpt_mm_iv`. When `km_threshold = TRUE`, initial estimates
#' for Vmax and Km are derived from the dataset to provide a representative
#' starting point for nonlinear elimination.
#'
#' @return A list containing parameter estimates and prediction error metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_1cmpt_mm_iv(
#'   dat = Bolus_1CPT,
#'   est.method = "focei",
#'   npdmm_inputcl = 4,
#'   npdmm_inputvd = 70,
#'   km_threshold = TRUE
#' )
#' }
#'
#' @seealso \code{\link{Fit_1cmpt_mm_iv}}
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

#' Run and evaluate a two-compartment IV model
#'
#' Fits a two-compartment intravenous pharmacokinetic model using a naive pooled
#' data approach and evaluates model performance based on prediction error
#' metrics.
#'
#' @param dat A data frame containing raw time–concentration data in standard
#'   nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param input.vc2cmpt Initial estimate for central compartment volume (Vc).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp2cmpt Initial estimate for peripheral compartment volume (Vp).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.q2cmpt Initial estimate for intercompartmental clearance (Q).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded prior to model fitting. The model is fitted
#' using `Fit_2cmpt_iv`, and prediction-based metrics are calculated to assess
#' performance.
#'
#' @return A list containing parameter estimates and prediction error metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_2cmpt_iv(dat = Bolus_2CPT,
#'                            input.cl = 4,
#'                            input.vc2cmpt = 35,
#'                            input.vp2cmpt = 35,
#'                            input.q2cmpt = 4)
#' }
#'
#' @seealso \code{\link{Fit_2cmpt_iv}}
#' @export

run_npd_2cmpt_iv <- function(dat,
                             est.method = "nls",
                             input.cl = exp(1),
                             input.vc2cmpt = exp(1),
                             input.vp2cmpt = exp(1),
                             input.q2cmpt = exp(1),
                             input.add = 1) {
  start.time <- Sys.time()

  npd.list.2cmpt <- NA
  npd_results_2cmpt <-
    data.frame(
      cl = NA,
      vc = NA,
      vp = NA,
      q = NA,
      timespent = NA
    )

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list.2cmpt <-  Fit_2cmpt_iv(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[1],
          1)
  npd.MAE <-
    round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[2],
          1)
  npd.MAPE <-
    round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[3],
          1)
  npd.RMSE <-
    round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[4],
          1)
  npd.rRMSE <-
    round(metrics.(pred.x = npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[5],
          1)

  return(
    list(
      npd.2cmpt_results = npd_results_2cmpt,
      npd.2cmpt.APE = npd.APE,
      npd.2cmpt.MAE  = npd.MAE,
      npd.2cmpt.MAPE = npd.MAPE,
      npd.2cmpt.RMSE = npd.RMSE,
      npd.2cmpt.rRMSE  = npd.rRMSE,
      npd.list.2cmpt = npd.list.2cmpt
    )
  )
}


#' Run and evaluate a three-compartment IV model
#'
#' Fits a three-compartment intravenous pharmacokinetic model using a naive
#' pooled data approach and evaluates model performance based on prediction
#' error metrics.
#'
#' @param dat A data frame containing raw intravenous concentration–time data in
#'   standard nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param input.vc3cmpt Initial estimate for central volume of distribution (Vc).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp3cmpt Initial estimate for first peripheral volume (Vp1).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp23cmpt Initial estimate for second peripheral volume (Vp2).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.q3cmpt Initial estimate for intercompartmental clearance between
#'   central and first peripheral compartments (Q1). Defaults to exp(1),
#'   corresponding to a log-scale value of 1.
#' @param input.q23cmpt Initial estimate for intercompartmental clearance between
#'   central and second peripheral compartments (Q2). Defaults to exp(1),
#'   corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded prior to model fitting. The model is
#' fitted using `Fit_3cmpt_iv`, and prediction-based metrics are calculated to
#' evaluate performance.
#'
#' @return A list containing fitted parameter estimates and model prediction
#' error metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_3cmpt_iv(
#'   dat = Bolus_2CPT,
#'   input.cl = 4,
#'   input.vc3cmpt = 10,
#'   input.vp3cmpt = 10,
#'   input.vp23cmpt = 10,
#'   input.q3cmpt = 4,
#'   input.q23cmpt = 4
#' )
#' }
#'
#' @seealso \code{\link{Fit_3cmpt_iv}}
#' @export

run_npd_3cmpt_iv <- function(dat,
                             est.method = "nls",
                             input.cl = exp(1),
                             input.vc3cmpt = exp(1),
                             input.vp3cmpt = exp(1),
                             input.vp23cmpt = exp(1),
                             input.q3cmpt =  exp(1),
                             input.q23cmpt = exp(1),
                             input.add = 1) {
  start.time <- Sys.time()
  npd.list.3cmpt <- NA
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

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list.3cmpt <-  Fit_3cmpt_iv(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[1],
          1)
  npd.MAE <-
    round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[2],
          1)
  npd.MAPE <-
    round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[3],
          1)
  npd.RMSE <-
    round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[4],
          1)
  npd.rRMSE <-
    round(metrics.(pred.x =  npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[5],
          1)


  return(
    list(
      npd.3cmpt_results = npd_results_3cmpt,
      npd.3cmpt.APE = npd.APE,
      npd.3cmpt.MAE  = npd.MAE,
      npd.3cmpt.MAPE = npd.MAPE,
      npd.3cmpt.RMSE = npd.RMSE,
      npd.3cmpt.rRMSE  = npd.rRMSE,
      npd.list.3cmpt =  npd.list.3cmpt
    )
  )
}


#' Run and evaluate a one-compartment oral model
#'
#' Fits a one-compartment oral pharmacokinetic model using a naive pooled data
#' approach and evaluates model performance based on prediction error metrics.
#'
#' @param dat A data frame containing raw time–concentration data in standard
#'   nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.ka Initial estimate for the absorption rate constant (Ka).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param input.vd Initial estimate for volume of distribution (Vd). Defaults to
#'   exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded before model fitting. The model is fitted
#' using `Fit_1cmpt_oral`, and prediction-based metrics are calculated to assess
#' performance.
#'
#' @return A list containing the fitted parameter estimates and prediction error
#' metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_1cmpt_oral(dat = Oral_1CPT, input.ka = 1, input.cl = 4, input.vd = 70)
#' }
#'
#' @seealso \code{\link{Fit_1cmpt_oral}}
#' @export


run_npd_1cmpt_oral <- function(dat,
                               est.method = "nls",
                               input.ka = exp(1),
                               input.cl = exp(1),
                               input.vd = exp(1),
                               input.add = 1) {
  start.time <- Sys.time()
  npd.list <- NA
  npd_results <-
    data.frame(
      ka = NA,
      cl = NA,
      vd = NA,
      timespent = NA
    )

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list <-  Fit_1cmpt_oral(
    data = dat[dat$EVID != 2,],
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
      ka = signif(npd.list$parFixedDf$`Back-transformed`[1], 3),
      cl = signif(npd.list$parFixedDf$`Back-transformed`[2], 3),
      vd = signif(npd.list$parFixedDf$`Back-transformed`[3], 3),
      timespent = time.spent
    )

  npd.APE <-
    round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV)[1], 1)
  npd.MAE <-
    round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV)[2], 1)
  npd.MAPE <-
    round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV)[3], 1)
  npd.RMSE <-
    round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV)[4], 1)
  npd.rRMSE <-
    round(metrics.(pred.x =  npd.list$cp, obs.y = npd.list$DV)[5], 1)

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


#' Run and evaluate a one-compartment oral model with Michaelis-Menten kinetics
#'
#' Fits a one-compartment oral pharmacokinetic model with Michaelis-Menten
#' elimination using a naive pooled data approach, and evaluates model
#' performance using prediction error metrics.
#'
#' @param dat A data frame containing raw time–concentration data in standard
#'   nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.ka Initial estimate for the absorption rate constant (ka).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vmax Initial estimate for the maximum metabolic rate (Vmax).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.km Initial estimate for the Michaelis constant (Km). Defaults to
#'   exp(1), corresponding to a log-scale value of 1.
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1. This value
#'   is used to derive initial Vmax and Km when `km_threshold = TRUE`.
#' @param input.vd Initial estimate for volume of distribution (Vd). Defaults to
#'   exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#' @param km_threshold Logical indicating whether initial Vmax and Km should be
#'   automatically adjusted based on observed maximum concentration and
#'   clearance. Defaults to FALSE.
#'
#' @details
#' The function excludes dosing records (`EVID == 2`) prior to model fitting.
#' When `km_threshold = TRUE`, initial estimates for Vmax and Km are derived
#' using the observed maximum concentration and clearance. The model is then
#' fitted using `Fit_1cmpt_mm_oral`, and prediction-based metrics are calculated
#' to assess performance.
#'
#' @return A list containing the fitted parameter estimates and prediction error
#' metrics.
#'
#' @examples
#' \dontrun{
#' run_npd_1cmpt_mm_oral(dat = Oral_1CPTMM, input.cl = 4, input.vd = 70,est.method="focei")
#' }
#'
#' @seealso \code{\link{Fit_1cmpt_mm_oral}}
#' @export

run_npd_1cmpt_mm_oral <- function(dat,
                                  est.method = "nls",
                                  input.ka = exp(1),
                                  input.vmax = exp(1),
                                  input.km = exp(1),
                                  input.cl = exp(1),
                                  input.vd = exp(1),
                                  input.add = 1,
                                  km_threshold = F) {
  start.time <- Sys.time()
  estvmax <- input.vmax
  estkm <- input.km

  npdmm.list <- NA
  npdmm_results <-
    data.frame(
      ka = NA,
      vmax =  NA,
      km =  NA,
      vd =  NA,
      timespent = NA
    )

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  # Initial estimates of Vmax and Km will be set based on threshold if km_threshold=T
  # Determine the maximum concentration
  if (km_threshold) {
    dat.obs <- dat[dat$EVID == 0,]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
    estmaxkm <-
      mean.pop.cmax * 4 # if km>>4cmax, it nearly fall into the linear range
    estkm <- mean.pop.cmax # initial km starts from cmax
    estvmax <-  estmaxkm * input.cl
  }

  npdmm.list <- Fit_1cmpt_mm_oral(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV)[1], 1)
  npd.MAE <-
    round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV)[2], 1)
  npd.MAPE <-
    round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV)[3], 1)
  npd.RMSE <-
    round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV)[4], 1)
  npd.rRMSE <-
    round(metrics.(pred.x =  npdmm.list$cp, obs.y = npdmm.list$DV)[5], 1)

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

#' Run and evaluate a two-compartment oral model
#'
#' Fits a two-compartment oral pharmacokinetic model using a naive pooled
#' data approach and evaluates model performance using prediction error metrics.
#'
#' @param dat A data frame containing time–concentration data in standard
#'   nlmixr2 format.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.ka Initial estimate for the absorption rate constant (Ka).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1)
#'  , corresponding to a log-scale value of 1.
#' @param input.vc2cmpt Initial estimate for the central volume of
#'   distribution (Vc). Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp2cmpt Initial estimate for the peripheral volume of
#'   distribution (Vp). Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.q2cmpt Initial estimate for intercompartmental clearance (Q).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded before fitting the model. The model is
#' fitted using `Fit_2cmpt_oral`, and prediction-based metrics are computed to
#' evaluate performance.
#'
#' @return A list containing the fitted parameter estimates and prediction error
#' metrics.
#'
#' @examples
#' \dontrun{
#' result <- run_npd_2cmpt_oral(
#'   dat = Oral_2CPT,
#'   input.ka = 1,
#'   input.cl = 4,
#'   input.vc2cmpt = 35,
#'   input.vp2cmpt = 35,
#'   input.q2cmpt = 4
#' )
#' }
#'
#' @seealso \code{\link{Fit_2cmpt_oral}}
#' @export

run_npd_2cmpt_oral <- function(dat,
                               est.method = "nls",
                               input.ka = exp(1),
                               input.cl = exp(1),
                               input.vc2cmpt = exp(1),
                               input.vp2cmpt = exp(1),
                               input.q2cmpt = exp(1),
                               input.add = 1) {
  start.time <- Sys.time()
  npd.list.2cmpt <- NA
  npd_results_2cmpt <-
    data.frame(
      ka = NA,
      cl = NA,
      vc = NA,
      vp = NA,
      q = NA,
      timespent = NA
    )

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list.2cmpt <-  Fit_2cmpt_oral(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV)[1],
          1)
  npd.MAE <-
    round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV)[2],
          1)
  npd.MAPE <-
    round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y = npd.list.2cmpt$DV)[3],
          1)
  npd.RMSE <-
    round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV)[4],
          1)
  npd.rRMSE <-
    round(metrics.(pred.x =   npd.list.2cmpt$cp, obs.y =  npd.list.2cmpt$DV)[5],
          1)

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


#' Run and evaluate a three-compartment oral model
#'
#' Fits a three-compartment oral pharmacokinetic model using a naive pooled data
#' approach and evaluates model performance using prediction error metrics.
#'
#' @param dat A data frame containing pharmacokinetic data in the standard
#'   nlmixr2 format, including required columns such as `ID`, `EVID`, `DV`, and
#'   `dose`.
#' @param est.method Estimation method used in nlmixr2. Defaults to "nls".
#' @param input.ka Initial estimate for the absorption rate constant (Ka).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.cl Initial estimate for clearance (CL). Defaults to exp(1),
#' corresponding to a log-scale value of 1.
#' @param input.vc3cmpt Initial estimate for central volume of distribution (Vc).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp3cmpt Initial estimate for first peripheral volume (Vp1).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.vp23cmpt Initial estimate for second peripheral volume (Vp2).
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.q3cmpt Initial estimate for intercompartmental clearance Q1.
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.q23cmpt Initial estimate for intercompartmental clearance Q2.
#'   Defaults to exp(1), corresponding to a log-scale value of 1.
#' @param input.add Additive error term. Defaults to 1.
#'
#' @details
#' Rows with `EVID == 2` are excluded prior to model fitting. The model is
#' fitted using `Fit_3cmpt_oral`, and prediction-based metrics are calculated to
#' assess model performance.
#'
#' @return A list containing parameter estimates and prediction error metrics.
#'
#' @examples
#' \dontrun{
#' result <- run_npd_3cmpt_oral(dat = Oral_3CPT)
#' }
#'
#' @seealso \code{\link{Fit_3cmpt_oral}}
#' @export

run_npd_3cmpt_oral <- function(dat,
                               est.method = "nls",
                               input.ka = exp(1),
                               input.cl = exp(1),
                               input.vc3cmpt = exp(1),
                               input.vp3cmpt = exp(1),
                               input.vp23cmpt = exp(1),
                               input.q3cmpt =  exp(1),
                               input.q23cmpt =  exp(1),
                               input.add = 1) {
  start.time <- Sys.time()
  npd.list.3cmpt <- NA
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

  npd.APE <- NA
  npd.MAE <- NA
  npd.MAPE <- NA
  npd.RMSE <- NA
  npd.rRMSE <- NA

  npd.list.3cmpt <-  Fit_3cmpt_oral(
    data = dat[dat$EVID != 2,],
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

  npd.APE <-
    round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV)[1],
          1)
  npd.MAE <-
    round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV)[2],
          1)
  npd.MAPE <-
    round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y = npd.list.3cmpt$DV)[3],
          1)
  npd.RMSE <-
    round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV)[4],
          1)
  npd.rRMSE <-
    round(metrics.(pred.x =   npd.list.3cmpt$cp, obs.y =  npd.list.3cmpt$DV)[5],
          1)


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
