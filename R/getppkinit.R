#' Get initial estimates for a population pharmacokinetic modelling
#'
#' Computes initial pharmacokinetic parameters using multiple methods including approximate calculation, non-compartmental analysis, graphical methods, and parameter estimation with naive pooled data approach if specified.
#' @param dat A data frame containing the pharmacokinetic data.
#' @param run.npdcmpt An integer indicating whether to run the naive pooled data approach (default is 0).
#' @param getinit.settings A list or data frame containing calculation settings (optional). The following settings can be provided:
#' \itemize{
#'   \item \code{half_life}: Numeric value for the half-life of the drug (default is \code{NA}). If not provided, it will be estimated based on the data.
#'   \item \code{nlastpoints}: Numeric value specifying the number of last data points used for linear regression to obtain the slope in the terminal phase. (default is \code{4}).
#'   \item \code{trap.rule.method}: Numeric value indicating the trapezoidal rule method to use (default is \code{1}).
#'   \code{1} refers to the linear trapezoidal method, while \code{2} refers to the linear-up/log-down trapezoidal method, which uses a linear method for the increasing phase and a logarithmic method for the decreasing phase.
#'   \item \code{nbins}: Numeric value specifying the number of bins to use when performing naive pooling of data (default is \code{8}).
#'   \item \code{est.method}: Character string indicating the estimation method for naive pooled data compartmental analysis to use (default is \code{"nls"}) or other methods ((e.g., "nls", "nlm", ""foce", "focei") depending on the analysis.
#'   \item \code{selection.criteria}: Character string indicating the selection criteria for method comparison (default is \code{"MAPE"}). The selection method can be set to either Mean Absolute Percentage Error (\code{"MAPE"}) or Absolute Percentage Error (\code{"APE"}), depending on the user requirement.
#'   \item \code{npdcmpt.inits.strategy}: Numeric value indicating the strategy for setting initial estimates in the model (default is \code{0}).  \code{0} means all initial estimates are set to \code{1} in this step, while \code{1} means that the initial estimates are based on parameters derived from non-model-fitting calculation methods.
#'
#' }
#' If any of these settings are not provided, default values will be used.
#' @return A list containing data information, initial parameter estimates, messages, and run history.
#' @importFrom dplyr %>% mutate filter select
#' @import nlmixr2
#' @importFrom tidyr fill
#' @import crayon
#' @examples
#' inits.out<-getppkinit(dat = Bolus_1CPT)
#' inits.out
#' @export

getppkinit <- function(dat,
                       run.npdcmpt=0,
                       getinit.settings=NULL) {

 # Default settings for getppkinit
  getinit.settings0 <- list(
    half_life = NA,
    nlastpoints = 4,
    trap.rule.method = 1,
    nbins = 8,
    est.method = "nls",
    selection.criteria="MAPE",
    npdcmpt.inits.strategy= 1
  )

  if (!is.null(getinit.settings)) {
    getinit.settings <- as.list(getinit.settings)

    for (name in names(getinit.settings)) {
      getinit.settings0[[name]] <- getinit.settings[[name]]
    }
  }


  # Extract the setting values
  half_life <- as.numeric(getinit.settings0$half_life)
  nlastpoints <- as.numeric(getinit.settings0$nlastpoints)
  trap.rule.method <- as.numeric(getinit.settings0$trap.rule.method)
  nbins <- as.numeric(getinit.settings0$nbins)
  est.method <- getinit.settings0$est.method
  selection.criteria<-getinit.settings0$selection.criteria
  npdcmpt.inits.strategy<-getinit.settings0$npdcmpt.inits.strategy
  ################## Data preprocessing and information summary#################

  message(black(
    paste0("Processing data", strrep(".", 20))
  ))

  # Calculate tad and dose number
  column_names <- toupper(colnames(dat))
  colnames(dat) <- toupper(colnames(dat))

  # Convert pharmacokinetic dataset to depreciated format with additional dosing
  if ("ADDL" %in% column_names) {
    dat <- nmpkconvert(dat)
  }

  # Check whether infusion case, currently
  infusion_flag <- 0
  infusion_flag_c <- "N"
  route<-"bolus"

  if ("DUR" %in% column_names) {
    if (!"RATE" %in% column_names) {
      dat$RATE <- 0
      dat[dat$DUR > 0,]$RATE <-
        dat[dat$DUR > 0,]$AMT / dat[dat$DUR > 0,]$DUR
    }
  }

  if ("RATE" %in% column_names) {
    infusion_flag <- 1
    infusion_flag_c <- "Y"
    route="infusion"
    message(black(
      paste0(
        "Infusion administration detected.",strrep(".", 20)
      )
    ))
  }

  if ("DOSE" %in% column_names) {
    dat$DOSE_PRE <- dat$DOSE
    dat$DOSE<-NULL
  }

  dat <- calculate_tad(dat, infusion_flag)

  # check whether non intravenous case
  noniv_flag <- 0
  noniv_flag_c <- "N"

  if ("CMT" %in% column_names) {
    if (length(unique(dat$CMT)) > 1) {
      noniv_flag <- 1
      noniv_flag_c <- "Y"
      route="oral"
      message(black(
        paste0(
          "Administration site detected to differ from measurement site; extravascular (oral) administration assumed."
        )
      ))
    }
  }

  sdflag <- 0
  sdflag_c <- "N"
  # check whether it is a single dose case
  if (nrow(dat[dat$EVID %in% c(1, 4, 101),]) == length(unique(dat[dat$EVID %in% c(1, 4, 101),]$ID))) {
    sdflag <- 1
    sdflag_c <- "Y"
  }

  # Check points with first dose interval
  fdobsflag <- 0
  fdobsflag_c <-  "N"
  if (nrow(dat[dat$dose_number == 1 & dat$EVID == 0,]) > nbins) {
    fdobsflag <- 1
    fdobsflag_c <- "Y"
  }

  # check characteristics of dataset
  nids <- nrow(dat[!duplicated(dat$ID),])
  nobs <- nrow(dat[dat$EVID == 0,])

  Datainfo <-
    paste0(
      "No. of subjects: ",
      nids,
      ", No. of observations: ",
      nobs,
      ", Dose route: ",
      route,
      ", Is single dose? ",
      sdflag_c,
      ", Data within the first dosing interval available? ",
      fdobsflag_c
    )

  ########################Half-life estimated ##############################

  message(black(
    paste0("Estimating half-life",strrep(".", 20))))

  # message(black(
  #   paste0("Performed linear regression on the terminal phase of pooled dose-normalized data, estimated half-life: ",
  #          half_life)))

   half_life<-half_life_estimated(dat = dat,
                                  sdflag = sdflag,
                                  fdobsflag = fdobsflag,
                                  nlastpoints = nlastpoints,
                                  nbins=nbins)

   message(black(
     paste0("Estimated half-life : ", half_life)))


  ##############Approximate calculation############################################

  message(black(
     paste0("Run approximate PK calculation on individual data",strrep(".", 20))))

  simpcal.APE <- Inf
  simpcal.MAPE <- Inf

  # iv case
  if (noniv_flag ==0){
    # half_life is estimated
    message(black(
      paste0("Try approximate calculation using estimated half-life: ",half_life,strrep(".", 20))))

    simpcal.results <- run_simpcal_iv(
    dat = dat,
    route = route,
    sdflag = sdflag,
    fdobsflag = fdobsflag,
    half_life = half_life)

    message(black(
      paste0("Try approximate calculation without no prior estimated half-life",strrep(".", 20))))

    # use most commonly used dose-interval to replace half life
    simpcal.results.2 <- run_simpcal_iv(
      dat = dat,
      route = route,
      sdflag = sdflag,
      fdobsflag = fdobsflag,
      half_life = NA
  )

  simpcal.out <- simpcal.results$simpcal.results
  simpcal.out.2 <- simpcal.results.2$simpcal.results

 # Evaluate the provided initial values of cl and vd by their goodness of fit
  if (is.na(simpcal.out$cl)==F & is.na(simpcal.out$vd)==F) {
    simpcal_sim <-
      Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = simpcal.out$cl,
        input.vd = simpcal.out$vd,
        input.add = 0
      )
    simpcal_sim$DV <- dat[dat$EVID == 0,]$DV
    # Absolute prediction error
    simpcal.APE <- sum(abs(simpcal_sim$cp - simpcal_sim$DV))
    simpcal.MAPE <-
      round(sum(abs(simpcal_sim$cp - simpcal_sim$DV) / simpcal_sim$DV) / nrow(simpcal_sim) *
              100, 1)
   }

  # Evaluate the provided initial values of cl and vd by their goodness of fit
  if (is.na(simpcal.out.2$cl)==F & is.na(simpcal.out.2$vd)==F) {
    simpcal_sim_2 <-
      Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = simpcal.out.2$cl,
        input.vd = simpcal.out.2$vd,
        input.add = 0
      )
    simpcal_sim_2$DV <- dat[dat$EVID == 0,]$DV
    # Absolute prediction error
    simpcal2.APE <- sum(abs(simpcal_sim_2$cp - simpcal_sim_2$DV))
    simpcal2.MAPE <-
      round(sum(abs(simpcal_sim_2$cp - simpcal_sim_2$DV) / simpcal_sim_2$DV) / nrow(simpcal_sim_2) *
              100, 1)
  }
}

  # Extravascular/oral  case

  if (noniv_flag ==1){
    # half_life is estimated
    message(black(
      paste0("Try approximate calculation using estimated half-life: ",half_life,strrep(".", 20))))

    simpcal.results <- run_simpcal_iv(
      dat = dat,
      route = route,
      sdflag = sdflag,
      fdobsflag = 0, # do not run vd part
      half_life = half_life)

    message(black(
      paste0("Try approximate calculation without no prior estimated half-life",strrep(".", 20))))

    # use most commonly used dose-interval to replace half life
    simpcal.results.2 <- run_simpcal_iv(
      dat = dat,
      route = route,
      sdflag = sdflag,
      fdobsflag = 0,
      half_life = NA
    )
    # obtained clearance
    simpcal.out <- simpcal.results$simpcal.results
    simpcal.out.2 <- simpcal.results.2$simpcal.results

  }

#################Naive pooled Non-compartmental analysis ##################################

  message(black(
    paste0("Run non-compartmental analysis on naive pooling data", strrep(".", 20))
  ))

 nca.results <- run_nca.normalised(
    dat = dat,
    nlastpoints = nlastpoints,
    trap.rule.method=trap.rule.method,
    nbins = nbins,
    fdobsflag = fdobsflag,
    sdflag=sdflag
  )


# Run extra ka part if oral case
  if ( noniv_flag==1 ){
    if (!is.null(nrow(nca.results$datpooled_fd$test.pool.normalised))){
     ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_fd$test.pool.normalised,
                      nlastpoints = nlastpoints,
                      nca.out = unlist(nca.results$nca.fd.results, use.names = FALSE))

     ka_method_1_fd <-ka_wanger_nelson_result$ka
     ka_method_1_out_fd<-ka_wanger_nelson_result$dat_out_wanger_nelson

  }

  if (!is.null(nrow(nca.results$datpooled_efd$test.pool.normalised))){
      ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_efd$test.pool.normalised,
                                                nlastpoints = nlastpoints,
                                                nca.out = unlist(nca.results$nca.efd.results, use.names = FALSE))
      ka_method_1_efd <-ka_wanger_nelson_result$ka
      ka_method_1_out_efd<-ka_wanger_nelson_result$dat_out_wanger_nelson

    }

   if (!is.null(nrow(nca.results$datpooled_all$test.pool.normalised))){
      ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_all$test.pool.normalised,
                                                nlastpoints = nlastpoints,
                                                nca.out = unlist(nca.results$nca.all.results, use.names = FALSE))
      ka_method_1_all <-ka_wanger_nelson_result$ka
      ka_method_1_out_all<-ka_wanger_nelson_result$dat_out_wanger_nelson

   }

    # can be used for later hybrid method
    ka_values<-c(ka_method_1_fd, ka_method_1_efd, ka_method_1_all)
    # Remove negative numbers
    positive_ka_values <-    ka_values[   ka_values > 0 &  is.na(ka_values)==F]
    ka_median <- round(median(positive_ka_values),2)

  }

  # Simulation and evaluation
  if (noniv_flag==0){
  # All pooled

  nca.results_all <- nca.results$nca.all.results
  nca.results_fd <- nca.results$nca.fd.results
  nca.results_efd <- nca.results$nca.efd.results

  nca.APE <- Inf
  nca.MAPE <- Inf

  if (is.na(nca.results_all$cl)==F & is.na(nca.results_all$vd )==F) {
  if (nca.results_all$cl > 0 & nca.results_all$vd > 0) {
    nca_sim <- Fit_1cmpt_iv(
      data = dat[dat$EVID != 2,],
      est.method = "rxSolve",
      input.cl = nca.results_all$cl,
      input.vd = nca.results_all$vd,
      input.add = 0
    )

    nca_sim$DV <- dat[dat$EVID == 0,]$DV
    nca.APE <- sum(abs(nca_sim$cp - nca_sim$DV))
    nca.MAPE <-
      round(sum(abs(nca_sim$cp - nca_sim$DV) / nca_sim$DV) / nrow(nca_sim) *
              100, 1)
  }
  }

  nca_fd.APE <- Inf
  nca_fd.MAPE <- Inf

  if (is.na(nca.results_fd$cl)==F & is.na(nca.results_fd$vd)==F) {
    if (nca.results_fd$cl > 0 & nca.results_fd$vd > 0) {
      nca_fd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = nca.results_fd$cl,
        input.vd = nca.results_fd$vd,
        input.add = 0
      )

      nca_fd_sim$DV <- dat[dat$EVID == 0,]$DV
      nca_fd.APE <- sum(abs(nca_fd_sim$cp - nca_fd_sim$DV))
      nca_fd.MAPE <-
        round(sum(abs(nca_fd_sim$cp - nca_fd_sim$DV) / nca_fd_sim$DV) / nrow(nca_fd_sim) *
                100, 1)
    }
  }


  nca_efd.APE <- Inf
  nca_efd.MAPE <- Inf

  if (is.na(nca.results_efd$cl)==F & is.na(nca.results_efd$vd)==F) {
    if (nca.results_efd$cl > 0 & nca.results_efd$vd > 0) {
      nca_efd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = nca.results_efd$cl,
        input.vd = nca.results_efd$vd,
        input.add = 0
      )

      nca_efd_sim$DV <- dat[dat$EVID == 0,]$DV
      nca_efd.APE <- sum(abs(nca_efd_sim$cp - nca_efd_sim$DV))
      nca_efd.MAPE <-
        round(sum(abs(nca_efd_sim$cp - nca_efd_sim$DV) / nca_efd_sim$DV) / nrow(nca_efd_sim) *
                100, 1)
    }
  }

}

  if (noniv_flag==1){
    # All pooled

    nca.results_all <- nca.results$nca.all.results
    nca.results_fd <- nca.results$nca.fd.results
    nca.results_efd <- nca.results$nca.efd.results

    nca.APE <- Inf
    nca.MAPE <- Inf

    if (is.na(nca.results_all$cl)==F & is.na(nca.results_all$vd)==F & is.na(ka_method_1_all)==F) {
    if (nca.results_all$cl > 0 & nca.results_all$vd > 0 &  ka_method_1_all>0 ) {
      nca_sim <- Fit_1cmpt_oral(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.ka = ka_method_1_all,
        input.cl = nca.results_all$cl,
        input.vd = nca.results_all$vd,
        input.add = 0
      )

      nca_sim$DV <- dat[dat$EVID == 0,]$DV
      nca.APE <- sum(abs(nca_sim$cp - nca_sim$DV))
      nca.MAPE <-
        round(sum(abs(nca_sim$cp - nca_sim$DV) / nca_sim$DV) / nrow(nca_sim) *
                100, 1)
    }
    }

    nca_fd.APE <- Inf
    nca_fd.MAPE <- Inf

    if (is.na(nca.results_fd$cl)==F & is.na(nca.results_fd$vd)==F & is.na(ka_method_1_fd)==F) {
      if (nca.results_fd$cl > 0 & nca.results_fd$vd > 0  &  ka_method_1_fd>0 ) {
        nca_fd_sim <- Fit_1cmpt_oral(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.ka = ka_method_1_fd,
          input.cl = nca.results_fd$cl,
          input.vd = nca.results_fd$vd,
          input.add = 0
        )

        nca_fd_sim$DV <- dat[dat$EVID == 0,]$DV
        nca_fd.APE <- sum(abs(nca_fd_sim$cp - nca_fd_sim$DV))
        nca_fd.MAPE <-
          round(sum(abs(nca_fd_sim$cp - nca_fd_sim$DV) / nca_fd_sim$DV) / nrow(nca_fd_sim) *
                  100, 1)
      }

    }

    nca_efd.APE <- Inf
    nca_efd.MAPE <- Inf

    if (is.na(nca.results_efd$cl)==F & is.na(nca.results_efd$vd)==F &  is.na(ka_method_1_efd)==F) {
      if (nca.results_efd$cl > 0 & nca.results_efd$vd > 0 &  ka_method_1_efd>0) {
        nca_efd_sim <- Fit_1cmpt_oral(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.ka = ka_method_1_efd,
          input.cl = nca.results_efd$cl,
          input.vd = nca.results_efd$vd,
          input.add = 0
        )

        nca_efd_sim$DV <- dat[dat$EVID == 0,]$DV
        nca_efd.APE <- sum(abs(nca_efd_sim$cp - nca_efd_sim$DV))
        nca_efd.MAPE <-
          round(sum(abs(nca_efd_sim$cp - nca_efd_sim$DV) / nca_efd_sim$DV) / nrow(nca_efd_sim) *
                  100, 1)
      }
    }

  }

##############Graphic analysis##################################################

# iv case
if ( noniv_flag==0 & fdobsflag==1 ){
 message(black(
    paste0("Run graphical analysis on naive pooling data only after the first dose", strrep(".", 20))
  ))

# infusion case, the infusion rate should be very fast
 # if (route=="infusion"){
 #    duration_value_counts<- table(dat[dat$EVID %in% c(1,4) &!duplicated(dat$ID), ]$AMT/dat[dat$EVID %in% c(1,4)&!duplicated(dat$ID), ]$RATE)
 #    most_commonly_used_infusion_duration<- as.numeric(names(which.max(duration_value_counts)))
 #  }
  graph_fd.APE <- Inf
  graph_fd.MAPE <- Inf

  graph.results_fd <- run_graphcal(
    dat = dat,
    route =   route ,
    nbins = nbins,
    nlastpoints = nlastpoints
  )

  if (is.na(graph.results_fd$cl)==F & is.na(graph.results_fd$vd)==F) {
    if (graph.results_fd$cl > 0 & graph.results_fd$vd > 0) {
      graph_fd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = graph.results_fd$cl,
        input.vd = graph.results_fd$vd,
        input.add = 0
      )

      graph_fd_sim$DV <- dat[dat$EVID == 0,]$DV
      graph_fd.APE <- sum(abs(graph_fd_sim$cp - graph_fd_sim$DV))
      graph_fd.MAPE <-
        round(sum(abs(graph_fd_sim$cp - graph_fd_sim$DV) / graph_fd_sim$DV) / nrow(graph_fd_sim) *
                100, 1)
    }

  }
}



# oral case
if ( noniv_flag==1  & fdobsflag==1  ){

  message(black(
    paste0("Run graphical analysis on naive pooling data only after the first dose", strrep(".", 20))
  ))

    graph.results_fd <- run_graphcal(
      dat = dat,
      route =   route ,
      nbins = nbins,
      nlastpoints = nlastpoints
    )

    graph_fd.APE <- Inf
    graph_fd.MAPE <- Inf
    if (is.na(graph.results_fd$cl)==F & is.na(graph.results_fd$vd)==F) {
      if (graph.results_fd$cl > 0 & graph.results_fd$vd > 0) {
        graph_fd_sim <- Fit_1cmpt_oral(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.ka = graph.results_fd$ka,
          input.cl = graph.results_fd$cl,
          input.vd = graph.results_fd$vd,
          input.add = 0
        )

        graph_fd_sim$DV <- dat[dat$EVID == 0,]$DV
        graph_fd.APE <- sum(abs(graph_fd_sim$cp - graph_fd_sim$DV))
        graph_fd.MAPE <-
          round(sum(abs(graph_fd_sim$cp - graph_fd_sim$DV) / graph_fd_sim$DV) / nrow(graph_fd_sim) *
                  100, 1)
      }

    }
  }

  ############Hybrid calculation##########################
  # volume of distribution was calculated by slope
  # Find an appropriate slope
  # if (nca_fd.MAPE < nca.MAPE) {
  #   usedslope <- nca.results_fd$slope
  # }
  #
  # if (nca_fd.MAPE >= nca.MAPE) {
  #   usedslope <- nca.results_all$slope
  # }

 message(black(
   paste0("Run hybrid calculation, combining approximate calculation methods and NCA ", strrep(".", 20))
 ))

  if (min(c(nca_fd.MAPE,nca_efd.MAPE,nca.MAPE))==nca_fd.MAPE) {
    usedslope <- nca.results_fd$slope
  }

  if (min(c(nca_fd.MAPE,nca_efd.MAPE,nca.MAPE))==nca_efd.MAPE) {
    usedslope <- nca.results_efd$slope
  }

  if (min(c(nca_fd.MAPE,nca_efd.MAPE,nca.MAPE))==nca.MAPE) {
    usedslope <- nca.results_all$slope
  }


  hybrid_cl <- NA
  hybrid_vd <- NA

  # Rapid calculation clearance + slope estimated
  if (is.na(simpcal.out$cl)==F) {
    hybrid_cl <- simpcal.out$cl
    hybrid_vd <-  signif(-simpcal.out$cl / usedslope, 3)
  }

#iv case
 if ( noniv_flag==0 ){
  hybrid.APE <- Inf
  hybrid.MAPE <- Inf
  if (is.na(hybrid_vd)==F & is.na(simpcal.out$cl)==F) {
    hybrid_sim <- Fit_1cmpt_iv(
      data = dat[dat$EVID != 2,],
      est.method = "rxSolve",
      input.cl = simpcal.out$cl,
      input.vd =   hybrid_vd,
      input.add = 0
    )

    hybrid_sim$DV <- dat[dat$EVID == 0,]$DV
    hybrid.APE <- sum(abs(hybrid_sim$cp - hybrid_sim$DV))
    hybrid.MAPE <-
      round(sum(abs(hybrid_sim$cp - hybrid_sim$DV) / hybrid_sim$DV) / nrow(hybrid_sim) *
              100, 1)

  }
 }

 if ( noniv_flag==1 ){
    hybrid.APE <- Inf
    hybrid.MAPE <- Inf
    if (is.na(hybrid_vd)==F & is.na(simpcal.out$cl)==F) {
      hybrid_sim <- Fit_1cmpt_oral(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.ka=ka_median,
        input.cl = simpcal.out$cl,
        input.vd =   hybrid_vd,
        input.add = 0
      )

      hybrid_sim$DV <- dat[dat$EVID == 0,]$DV
      hybrid.APE <- sum(abs(hybrid_sim$cp - hybrid_sim$DV))
      hybrid.MAPE <-
        round(sum(abs(hybrid_sim$cp - hybrid_sim$DV) / hybrid_sim$DV) / nrow(hybrid_sim) *
                100, 1)

    }
  }

  # Output all of test
  if (noniv_flag==0){
     ka=c(NA,NA,NA,NA,NA,NA)
  }

  if (noniv_flag==1){
    ka=c(NA,graph.results_fd$ka,ka_method_1_fd,ka_method_1_efd,ka_method_1_all,ka_median)
  }

  all.out <- data.frame(
    method = c(
      "Approximate calculation",
      "Graphic calculation",
      "NCA (only first dose)",
      "NCA (data exclude first-dose part)",
      "NCA (all pooled)",
      "Hybrid approximate calculation"
    ),
    ka=ka,
    cl = c(
      simpcal.out$cl,
      graph.results_fd$cl,
      nca.results_fd$cl,
      nca.results_efd$cl,
      nca.results_all$cl,
      hybrid_cl
    ),
    vd = c(
      simpcal.out$vd,
      graph.results_fd$vd,
      nca.results_fd$vd,
      nca.results_efd$vd,
      nca.results_all$vd,
      hybrid_vd
    ),
    simAPE = c(simpcal.APE, graph_fd.APE, nca_fd.APE, nca_efd.APE,nca.APE, hybrid.APE),
    simMAPE = c(
      simpcal.MAPE,
      graph_fd.MAPE,
      nca_fd.MAPE,
      nca_efd.MAPE,
      nca.MAPE,
      hybrid.MAPE
    ),
    time.spent = c(
      simpcal.out$time.spent,
      graph.results_fd$time.spent,
      nca.results_fd$time.spent,
      nca.results_efd$time.spent,
      nca.results_all$time.spent,
      simpcal.out$time.spent
    )
  )

  # For single dose case, only report the pooled results of the single dose,
  # as it is the same as the results when all data are pooled together.

  if (sdflag == 1) {
    # Output all of test
    if (noniv_flag==0){
      ka=c(NA,NA,NA,NA)
    }

    if (noniv_flag==1){
      ka=c(NA,graph.results_fd$ka,ka_method_1_fd,ka_median)
    }

    all.out <- data.frame(
      method = c(
        "Approximate calculation",
        "Graphic calculation",
        "NCA (only first dose)",
        "Hybrid approximate calculation"
      ),
      ka=ka,
      cl = c(
        simpcal.out$cl,
        graph.results_fd$cl,
        nca.results_fd$cl,
        hybrid_cl
      ),
      vd = c(
        simpcal.out$vd,
        graph.results_fd$vd,
        nca.results_fd$vd,
        hybrid_vd
      ),
      simAPE = c(simpcal.APE, graph_fd.APE, nca.APE, hybrid.APE),
      simMAPE = c(simpcal.MAPE, graph_fd.MAPE, nca.MAPE, hybrid.MAPE),
      time.spent = c(
        simpcal.out$time.spent,
        graph.results_fd$time.spent,
        nca.results_all$time.spent,
        simpcal.out$time.spent
      )
    )
  }

  colnames(all.out) <-
    c(
      "Method",
      "Calculated Ka",
      "Calculated CL",
      "Calculated Vd",
      "Absolute Prediction Error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )


  all.out$`Mean absolute prediction error (MAPE)` <-
    round(all.out$`Mean absolute prediction error (MAPE)`, 1)

  all.out$`Absolute Prediction Error (APE)`<-
    round(all.out$`Absolute Prediction Error (APE)`, 1)

  if (selection.criteria=="APE"){
  base.ka.best <-
    all.out[all.out$`Absolute Prediction Error (APE)` == min(all.out$`Absolute Prediction Error (APE)`),]$`Calculated Ka`

  base.cl.best <-
    all.out[all.out$`Absolute Prediction Error (APE)` == min(all.out$`Absolute Prediction Error (APE)`),]$`Calculated CL`

  base.vd.best <-
    all.out[all.out$`Absolute Prediction Error (APE)`== min(all.out$`Absolute Prediction Error (APE)`),]$`Calculated Vd`
  }

  if (selection.criteria=="MAPE"){
  base.ka.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated Ka`

  base.cl.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated CL`

  base.vd.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated Vd`
}

message(black(
    paste0("Base PK parameter analysis finished. Estimated ka :",base.ka.best, ", estimated CL : ", base.cl.best, ", estimated Vd : ", base.vd.best )))


##############################Vmax and Km estimation##########################

  message(black(
    paste0("Run sensitivity analysis on nonlinear eliminiation kinetics PK parameters",strrep(".", 20))))

  # For vmax, km, if no model fitting for parameter estimation is used,
  # sensitivity analysis via simulation to select appropriate initial estimates
  max.dose <-
    max(dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0, ]$dose)[1]

  # sort out the cmax during the whole pk course
  # observed maximum cmax (median value of population)
  dat.obs <- dat[dat$EVID == 0,]
  pop.cmax <- aggregate(dat.obs$DV,
                        list(dat.obs$ID),
                        FUN = max,
                        na.rm = T)
  mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
  fcmax<-  mean.pop.cmax
  # calc.cmax <- max.dose / min(base.vd.best)[1]
  # fcmax <- max(c(mean.pop.cmax, calc.cmax))[1]
  # fcmax <- max(c(mean.pop.cmax, calc.cmax))[1]
  # if (fcmax == mean.pop.cmax) {
  #   cmax.message = "Cmax (mean maximum concontration) was obtained from mean observed cmax of analysis dataset"
  # }
  #
  # if (fcmax == calc.cmax) {
  #   cmax.message = "Cmax was calculated from maximum dose divided by volume of distribution"
  # }
  linear_minkm <- fcmax * 4 # if km>>4cmax, it nearly fall into the linear range

  sim.vmax.km.results.all <- NULL

  if (noniv_flag==0){
  for (besti in 1:length(base.cl.best)) {
    sim.vmax.km.results.all.i <- sim_sens_vmax_km(
      dat = dat,
      estcmax =  fcmax,
      estcl = base.cl.best[besti],
      estvd = base.vd.best[besti]
    )
    sim.vmax.km.results.all <-
      rbind(sim.vmax.km.results.all, sim.vmax.km.results.all.i)
    rownames(sim.vmax.km.results.all) <-
      seq(1, nrow(sim.vmax.km.results.all), 1)
   }
  }
  if (noniv_flag==1){
    for (besti in 1:length(base.cl.best)) {
      sim.vmax.km.results.all.i <- sim_sens_vmax_km(
        dat = dat,
        estcmax =  fcmax,
        estcl = base.cl.best[besti],
        estvd = base.vd.best[besti],
        estka = base.ka.best[besti],
        noniv_flag = 1
      )
      sim.vmax.km.results.all <-
        rbind(sim.vmax.km.results.all, sim.vmax.km.results.all.i)
      rownames(sim.vmax.km.results.all) <-
        seq(1, nrow(sim.vmax.km.results.all), 1)
    }

  }

  recommended_vmax_init <-
    sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE == min(sim.vmax.km.results.all$sim.mm.MAPE),]$vmax[1]
  recommended_km_init <-
    sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE == min(sim.vmax.km.results.all$sim.mm.MAPE),]$km[1]

  # message(black(
  #   paste0("Nonlinear elimination parameter analysis finished. Estimated Vmax : ",   recommended_vmax_init, ", estimated km : ", recommended_km_init  )))

  ###########Multi-Compartmental Model Parameter Analysis#######################

  message(black(
    paste0("Run sensitivity analysis on multi-compartmental PK parameters",strrep(".", 20))))

  # Default. simulation test
  sim.2cmpt.results.all <- NULL

  if (noniv_flag==0){
  for (besti in 1:length(base.cl.best)) {
    sim.2cmpt.results.all.i <- sim_sens_2cmpt(dat = dat,
                                              estcl = base.cl.best[besti],
                                              estvd = base.vd.best[besti])
    sim.2cmpt.results.all <-
      rbind(sim.2cmpt.results.all, sim.2cmpt.results.all.i)
    rownames(sim.2cmpt.results.all) <-
      seq(1, nrow(sim.2cmpt.results.all), 1)
  }
  }


  if (noniv_flag==1){
    for (besti in 1:length(base.cl.best)) {
      sim.2cmpt.results.all.i <- sim_sens_2cmpt(dat = dat,
                                                estcl = base.cl.best[besti],
                                                estvd = base.vd.best[besti],
                                                estka = base.ka.best[besti],
                                                noniv_flag = 1)
      sim.2cmpt.results.all <-
        rbind(sim.2cmpt.results.all, sim.2cmpt.results.all.i)
      rownames(sim.2cmpt.results.all) <-
        seq(1, nrow(sim.2cmpt.results.all), 1)
    }
  }

  if (selection.criteria=="MAPE"){
  recommended_vc2cmpt_init <-
    sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE),]$vc[1]
  recommended_vp2cmpt_init <-
    sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE),]$vp[1]
  recommended_q2cmpt_init <-
    sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE),]$q[1]
  }

  if (selection.criteria=="APE"){
    recommended_vc2cmpt_init <-
      sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE== min(sim.2cmpt.results.all$sim.2cmpt.APE),]$vc[1]
    recommended_vp2cmpt_init <-
      sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE== min(sim.2cmpt.results.all$sim.2cmpt.APE),]$vp[1]
    recommended_q2cmpt_init <-
      sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE == min(sim.2cmpt.results.all$sim.2cmpt.APE),]$q[1]
  }

  sim.3cmpt.results.all <- NULL

  if (noniv_flag==0){
  for (besti in 1:length(base.cl.best)) {
    sim.3cmpt.results.all.i <- sim_sens_3cmpt(dat = dat,
                                              estcl = base.cl.best[besti],
                                              estvd = base.vd.best[besti])
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results.all.i)
    rownames(sim.3cmpt.results.all) <-
      seq(1, nrow(sim.3cmpt.results.all), 1)
  }
  }

  if (noniv_flag==1){
  for (besti in 1:length(base.cl.best)) {
    sim.3cmpt.results.all.i <- sim_sens_3cmpt(dat = dat,
                                              estcl = base.cl.best[besti],
                                              estvd = base.vd.best[besti],
                                              estka = base.ka.best[besti],
                                              noniv_flag = 1)
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results.all.i)
    rownames(sim.3cmpt.results.all) <-
      seq(1, nrow(sim.3cmpt.results.all), 1)
  }
  }

  if (selection.criteria=="MAPE"){
  recommended_vc3cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vc[1]
  recommended_vp3cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vp[1]
  recommended_vp23cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vp2[1]
  recommended_q3cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$q[1]
  recommended_q23cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$q2[1]
  }

  if (selection.criteria=="APE"){
    recommended_vc3cmpt_init <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE == min(sim.3cmpt.results.all$sim.3cmpt.APE),]$vc[1]
    recommended_vp3cmpt_init <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE == min(sim.3cmpt.results.all$sim.3cmpt.APE),]$vp[1]
    recommended_vp23cmpt_init <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE == min(sim.3cmpt.results.all$sim.3cmpt.APE),]$vp2[1]
    recommended_q3cmpt_init <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE == min(sim.3cmpt.results.all$sim.3cmpt.APE),]$q[1]
    recommended_q23cmpt_init <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE == min(sim.3cmpt.results.all$sim.3cmpt.APE),]$q2[1]
  }

  # Remove these temporary global variables run before.
  # List of variables to remove
  vars_to_remove <-
    c(
      "input.add",
      "input.ka",
      "input.cl",
      "input.vc2cmpt",
      "input.vc3cmpt",
      "input.vp2cmpt",
      "input.vp3cmpt",
      "input.vp23cmpt",
      "input.q2cmpt",
      "input.q3cmpt",
      "input.q23cmpt",
      "input.vmax",
      "input.km",
      "input.vd"

    )

  # Check if variables exist and remove them
  vars_to_remove <-
    vars_to_remove[vars_to_remove %in% ls(envir = .GlobalEnv)]
  rm(list = vars_to_remove, envir = .GlobalEnv)

#####################Naive pooled data approach compartmental analysis #######
  if (run.npdcmpt==1){
    message(black(
      paste0("Run naive pooled data compartmental analysis ",strrep(".", 20))))

    if (npdcmpt.inits.strategy==0){
      message(red(
        paste0("Warning: there is no reference for intial estimates in current naive pooled data compartmental analysis, running maybe crashed down due to failure convergence")))
      input.cl = 1
      input.vd = 1
      input.vc2cmpt=  1
      input.vp2cmpt= 1
      input.vc3cmpt =  1
      input.vp3cmpt = 1
      input.vp23cmpt = 1
      input.q2cmpt=    1
      input.q3cmpt =   1
      input.q23cmpt = 1
      input.vmax=1
      input.km=1
      vmax_km_threshold=F
    }

    if (npdcmpt.inits.strategy==1){
    input.cl = mean(base.cl.best)
    input.vd = mean(base.vd.best)
    input.vc2cmpt=  recommended_vc2cmpt_init
    input.vp2cmpt= recommended_vp2cmpt_init
    input.vc3cmpt =  recommended_vc3cmpt_init
    input.vp3cmpt =  recommended_vp3cmpt_init
    input.vp23cmpt = recommended_vp23cmpt_init
    input.q2cmpt=    recommended_q2cmpt_init
    input.q3cmpt =   recommended_q3cmpt_init
    input.q23cmpt = recommended_q23cmpt_init
    input.vmax=  recommended_vmax_init
    input.km=  recommended_km_init
    vmax_km_threshold=T
    }


  if (noniv_flag==0){

    message(black(
      paste0("Run one-compartment model with first-order elimination", strrep(".", 20))
    ))


    npd_1cmpt_out <- run_npd_1cmpt_iv(
      dat = dat,
      est.method = est.method,
      input.cl =  input.cl ,
      input.vd =  input.vd
    )

    message(black(
      paste0("Run one-compartment model with Michaelis–Menten elimination", strrep(".", 20))
    ))

    if (vmax_km_threshold==T){
    npd_1cmpt_mm_out <- run_npd_1cmpt_mm_iv(
      dat = dat,
      est.method = est.method,
      npdmm_inputcl =  input.cl,
      npdmm_inputvd =  input.vd
    )
    }

    if (vmax_km_threshold==F){
      npd_1cmpt_mm_out <- run_npd_1cmpt_mm_iv(
        dat = dat,
        est.method = est.method,
        npdmm_inputvmax =  input.vmax,
        npdmm_inputkm =  input.km,
        npdmm_inputvd =  input.vd,
        km_threshold=F
      )
    }

    message(black(
      paste0("Run two-compartment model with first-order elimination", strrep(".", 20))
    ))

    npd_2cmpt_out <- run_npd_2cmpt_iv(
      dat = dat,
      est.method = est.method,
      input.cl=input.cl,
      input.vc2cmpt=  input.vc2cmpt,
      input.vp2cmpt=  input.vp2cmpt,
      input.q2cmpt=  input.q2cmpt
    )

    message(black(
      paste0("Run three-compartment model with first-order elimination", strrep(".", 20))
    ))

    npd_3cmpt_out <- run_npd_3cmpt_iv(
      dat = dat,
      est.method = est.method,
      input.cl=  input.cl,
      input.vc3cmpt =  input.vc3cmpt,
      input.vp3cmpt =  input.vp3cmpt  ,
      input.vp23cmpt =  input.vp23cmpt,
      input.q3cmpt =     input.q3cmpt ,
      input.q23cmpt =   input.q23cmpt)
  }


  if (noniv_flag==1){

    input.ka = mean(base.ka.best)

    message(black(
      paste0("Run one-compartment model with first-order absorption and elimination", strrep(".", 20))
    ))

    npd_1cmpt_out <- run_npd_1cmpt_oral(
      dat = dat,
      est.method = est.method,
      input.ka =  input.ka,
      input.cl = input.cl,
      input.vd =  input.vd
    )

    message(black(
      paste0("Run one-compartment model with first-order absorption and Michaelis–Menten elimination", strrep(".", 20))
    ))

    if (vmax_km_threshold==T){
    npd_1cmpt_mm_out <- run_npd_1cmpt_mm_oral(
      dat = Oral_1CPT,
      est.method = "nls",
      input.ka =  input.ka,
      input.cl =  input.cl,
      input.vd = input.vd,
      km_threshold = T
    )
    }

    if (vmax_km_threshold==F){
    npd_1cmpt_mm_out <- run_npd_1cmpt_mm_oral(
      dat = Oral_1CPT,
      est.method = "nls",
      input.ka =  input.ka,
      input.vmax =  input.vmax,
      input.km = input.km,
      km_threshold = F
    )
    }

    message(black(
      paste0("Run two-compartment model with first-order absorption and elimination", strrep(".", 20))
    ))

    npd_2cmpt_out <- run_npd_2cmpt_oral(
      dat = dat,
      est.method = est.method,
      input.cl=input.cl,
      input.vc2cmpt=    input.vc2cmpt,
      input.vp2cmpt=   input.vp2cmpt,
      input.q2cmpt= input.q2cmpt,
    )

    message(black(
      paste0("Run three-compartment model with first-order absorption and elimination", strrep(".", 20))
    ))

    npd_3cmpt_out <- run_npd_3cmpt_oral(
      dat = dat,
      est.method = est.method,
      input.cl=input.cl,
      input.vc3cmpt =  input.vc3cmpt,
      input.vp3cmpt =   input.vp3cmpt ,
      input.vp23cmpt =  input.vp23cmpt,
      input.q3cmpt =    input.q3cmpt  ,
      input.q23cmpt =  input.q23cmpt)
  }

# summarise results
    npd.1cmpt_results <- npd_1cmpt_out$npd.1cmpt_results
    npd.1cmpt.APE <- npd_1cmpt_out$npd.1cmpt.APE
    npd.1cmpt.MAPE <- npd_1cmpt_out$npd.1cmpt.MAPE

    npd.1cmpt.mm_results <- npd_1cmpt_mm_out$npd.1cmpt.mm_results
    npd.1cmpt.mm.APE <- npd_1cmpt_mm_out$npd.1cmpt.mm.APE
    npd.1cmpt.mm.MAPE <- npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE

    npd.2cmpt_results <- npd_2cmpt_out$npd.2cmpt_results
    npd.2cmpt.APE <- npd_2cmpt_out$npd.2cmpt.APE
    npd.2cmpt.MAPE <- npd_2cmpt_out$npd.2cmpt.MAPE

    npd.3cmpt_results <- npd_3cmpt_out$npd.3cmpt_results
    npd.3cmpt.APE <- npd_3cmpt_out$npd.3cmpt.APE
    npd.3cmpt.MAPE <- npd_3cmpt_out$npd.3cmpt.MAPE


   # Output
   # Naive pooled data approach (compartmental analysis)"
    if (noniv_flag==0){
    npdcmpt.all.out <- data.frame(

      Model = c("One-compartment & first-order (absorption) and elimination",
                "One-compartment & first-order (absorption) and nonlinear elimination",
                "Two-compartment & first-order (absorption) and elimination",
                "Three-compartment with first-order (absorption) and elimination"),

      ka= c(NA,
            npd.1cmpt.mm_results$ka,
            npd.2cmpt_results$ka,
            npd.3cmpt_results$ka
      ),

      cl = c(npd.1cmpt_results$cl,
             NA,
             npd.2cmpt_results$cl,
             npd.3cmpt_results$cl),

      vmax = c(NA,
             npd.1cmpt.mm_results$vmax,
             NA,
             NA),

       km = c(NA,
              npd.1cmpt.mm_results$km,
               NA,
               NA),

      vc = c(npd.1cmpt_results$vd,
             npd.1cmpt.mm_results$vd,
             npd.2cmpt_results$vc,
             npd.3cmpt_results$vc),

      vp =  c(NA,
              NA,
              npd.2cmpt_results$vp,
              npd.3cmpt_results$vp),

      vp2 =  c(NA,
               NA,
               NA,
               npd.3cmpt_results$vp2),

      q =   c(NA,
              NA,
              npd.2cmpt_results$q,
              npd.3cmpt_results$q),

      q2 =  c(NA,
              NA,
              NA,
              npd.3cmpt_results$q2),

      simAPE = c(npd.1cmpt.APE, npd.1cmpt.mm.APE, npd.2cmpt.APE, npd.3cmpt.APE),
      simMAPE = c(npd.1cmpt.MAPE, npd.1cmpt.mm.MAPE, npd.2cmpt.MAPE, npd.3cmpt.MAPE),

      time.spent = c(
        npd.1cmpt_results$timespent,
        npd.1cmpt.mm_results$timespent,
        npd.2cmpt_results$timespent,
        npd.3cmpt_results$timespent
      )
    )
    }


    if (noniv_flag==1){
    npdcmpt.all.out <- data.frame(

        Model = c("One-compartment & first-order (absorption) and elimination",
                  "One-compartment & first-order (absorption) and nonlinear elimination",
                  "Two-compartment & first-order (absorption) and elimination",
                  "Three-compartment with first-order (absorption) and elimination"),

        ka= c(npd.1cmpt_results$ka,
              npd.1cmpt.mm_results$ka,
              npd.2cmpt_results$ka,
              npd.3cmpt_results$ka
              ),

        cl = c(npd.1cmpt_results$cl,
               NA,
               npd.2cmpt_results$cl,
               npd.3cmpt_results$cl),


        vmax = c(NA,
                 npd.1cmpt.mm_results$vmax,
                 NA,
                 NA),

        km = c(NA,
               npd.1cmpt.mm_results$km,
               NA,
               NA),

        vc = c(npd.1cmpt_results$vd,
               npd.1cmpt.mm_results$vd,
               npd.2cmpt_results$vc,
               npd.3cmpt_results$vc),

        vp =  c(NA,
                NA,
                npd.2cmpt_results$vp,
                npd.3cmpt_results$vp),

        vp2 =  c(NA,
                NA,
                NA,
                npd.3cmpt_results$vp2),

        q =   c(NA,
                 NA,
                 npd.2cmpt_results$q,
                 npd.3cmpt_results$q),

        q2 =  c(NA,
                 NA,
                 NA,
                 npd.3cmpt_results$q2),

        simAPE = c(npd.1cmpt.APE, npd.1cmpt.mm.APE, npd.2cmpt.APE, npd.3cmpt.APE),
        simMAPE = c(npd.1cmpt.MAPE, npd.1cmpt.mm.MAPE, npd.2cmpt.MAPE, npd.3cmpt.MAPE),

        time.spent = c(
          npd.1cmpt_results$timespent,
          npd.1cmpt.mm_results$timespent,
          npd.2cmpt_results$timespent,
          npd.3cmpt_results$timespent
        )
      )

}
    # colnames(all.out) <-
    #   c(
    #     "Method",
    #     "Calculated Ka",
    #     "Calculated CL",
    #     "Calculated Vd",
    #     "Absolute Prediction Error (APE)",
    #     "Mean absolute prediction error (MAPE)",
    #     "Time spent"
    #   )

  }

  ######################## Finally selection########################################

   # ka,cl,vd.
    f_init_ka <- base.ka.best[1]
    f_init_cl <- base.cl.best[1]
    f_init_vd <- base.vd.best[1]

    if (selection.criteria=="MAPE"){
    sel.method.ka.cl.vd <-
      all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`,
                                                                               na.rm = T),]$Method[1]

       if (length(as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                  min(all.out$`Mean absolute prediction error (MAPE)`,
                                      na.rm = T),]$`Calculated CL`)) > 1) {
      #check the total volume of distribution of two compartment
      total.vd <-
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vc[1] +
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vp[1]


      f_init_ka <-
        as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                             min(all.out$`Mean absolute prediction error (MAPE)` ,
                                 na.rm = T) &
                             round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Ka`)[1]

      f_init_cl <-
        as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                             min(all.out$`Mean absolute prediction error (MAPE)` ,
                                 na.rm = T) &
                             round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated CL`)[1]

      f_init_vd <-
        as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                             min(all.out$`Mean absolute prediction error (MAPE)` ,
                                 na.rm = T) &
                             round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Vd`)[1]

      sel.method.ka.cl.vd <-
        all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)` ,
                                                                       na.rm = T) &
                  round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$Method[1]
    }

    }

    if (selection.criteria=="APE"){
      sel.method.ka.cl.vd <-
        all.out[all.out$`Absolute Prediction Error (APE)`== min(all.out$`Absolute Prediction Error (APE)`,
                                                                       na.rm = T),]$Method[1]

      if (length(as.numeric(all.out[all.out$`Absolute Prediction Error (APE)` ==
                                    min(all.out$`Absolute Prediction Error (APE)`,
                                        na.rm = T),]$`Calculated CL`)) > 1) {
        #check the total volume of distribution of two compartment
        total.vd <-
          sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE == min(sim.2cmpt.results.all$sim.2cmpt.APE, na.rm = T),]$vc[1] +
          sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE == min(sim.2cmpt.results.all$sim.2cmpt.APE, na.rm = T),]$vp[1]


        f_init_ka <-
          as.numeric(all.out[all.out$`Absolute Prediction Error (APE)` ==
                               min(all.out$`Absolute Prediction Error (APE)` ,
                                   na.rm = T) &
                               round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Ka`)[1]

        f_init_cl <-
          as.numeric(all.out[all.out$`Absolute Prediction Error (APE)`==
                               min(all.out$`Absolute Prediction Error (APE)`,
                                   na.rm = T) &
                               round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated CL`)[1]

        f_init_vd <-
          as.numeric(all.out[all.out$`Absolute Prediction Error (APE)` ==
                               min(all.out$`Absolute Prediction Error (APE)` ,
                                   na.rm = T) &
                               round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Vd`)[1]

        sel.method.ka.cl.vd <-
          all.out[all.out$`Absolute Prediction Error (APE)` == min(all.out$`Absolute Prediction Error (APE)` ,
                                                                         na.rm = T) &
                    round(all.out$`Calculated Vd`, 1) == round(total.vd, 1),]$Method[1]
      }

    }



    # vmax.km
    f_init_vmax <- recommended_vmax_init
    f_init_km <- recommended_km_init
    sel.method.vmax.km <- "Sensitivity analysis by simulation "

    # Multi-compartmental parameters
    f_init_vc2cmpt <-  recommended_vc2cmpt_init
    f_init_vp2cmpt <-  recommended_vp2cmpt_init
    f_init_q2cmpt <-   recommended_q2cmpt_init
    f_init_vc3cmpt <- recommended_vc3cmpt_init
    f_init_vp3cmpt <- recommended_vp3cmpt_init
    f_init_vp23cmpt <- recommended_vp23cmpt_init
    f_init_q3cmpt <- recommended_q3cmpt_init
    f_init_q23cmpt <-recommended_q23cmpt_init

    sel.method.multi <- "Sensitivity analysis by simulation "

    sel.method.ka<-"Wanger_nelson"

    if (sel.method.ka.cl.vd== "Hybrid approximate calculation" & sdflag==0 & fdobsflag==0 ){
        sel.method.ka<-"Wanger_nelson (median)"
    }

    if (sel.method.ka.cl.vd== "Graphic calculation"){
      sel.method.ka<-"Methods of residuals"
    }

    init.params.out.ka <- data.frame(method = sel.method.ka,
                                     vd = f_init_ka)

    init.params.out.cl <- data.frame(method = sel.method.ka.cl.vd,
                                     vd = f_init_cl)

    init.params.out.vd <- data.frame(method = sel.method.ka.cl.vd,
                                     vd = f_init_vd)

    init.params.out.vmax.km <-
      data.frame(method = sel.method.vmax.km,
                 vmax = f_init_vmax,
                 km =  f_init_km)

    init.params.out.vc.vp <- data.frame(
      method  =  sel.method.multi,
      vc2cmpt =  f_init_vc2cmpt,
      vp2cmpt =   f_init_vp2cmpt,
      q2cmpt =   f_init_q2cmpt,
      vc3cmpt =  f_init_vc3cmpt,
      vp3cmpt =  f_init_vp3cmpt,
      vp23cmpt =   f_init_vp23cmpt,
      q3cmpt =   f_init_q3cmpt,
      q23cmpt =   f_init_q23cmpt
    )


  colnames(init.params.out.ka) <- c("Method", "Ka")
  colnames(init.params.out.cl) <- c("Method", "CL")
  colnames(init.params.out.vd) <- c("Method", "Vd")
  colnames(init.params.out.vmax.km) <- c("Method", "Vmax", "Km")
  colnames(init.params.out.vc.vp) <-
    c("Method",
      "Vc2cmpt",
      "Vp2cmpt",
      "Q2cmpt",
      "Vc3cmpt",
      "Vp3cmpt",
      "Vp23cmpt",
      "Q3cmpt",
      "Q23cmpt"
      )

  init.params.out.all <- list(
    init.params.ka = init.params.out.ka,
    init.params.cl = init.params.out.cl,
    init.params.vd = init.params.out.vd,
    init.params.vmax.km = init.params.out.vmax.km,
    init.params.multi = init.params.out.vc.vp
  )

  Recommended_inits_df = data.frame(
    Parameters = c(
      "Ka",
      "CL",
      "Vd",
      "Vmax",
      "Km",
      "Vc(2CMPT)",
      "Vp(2CMPT)",
      "Q(2CMPT)",
      "Vc(3CMPT)",
      "Vp(3CMPT)",
      "Vp2(3CMPT)",
      "Q(3CMPT)",
      "Q2(3CMPT)"
    ),

    Methods = c(
      init.params.out.all$init.params.ka$Method,
      init.params.out.all$init.params.cl$Method,
      init.params.out.all$init.params.vd$Method,
      init.params.out.all$init.params.vmax.km$Method,
      init.params.out.all$init.params.vmax.km$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method,
      init.params.out.all$init.params.multi$Method

    ),
    Values = c(
      init.params.out.all$init.params.ka$Ka,
      init.params.out.all$init.params.cl$CL,
      init.params.out.all$init.params.vd$Vd,
      init.params.out.all$init.params.vmax.km$Vmax,
      init.params.out.all$init.params.vmax.km$Km,
      init.params.out.all$init.params.multi$Vc2cmpt,
      init.params.out.all$init.params.multi$Vp2cmpt,
      init.params.out.all$init.params.multi$Q2cmpt,
      init.params.out.all$init.params.multi$Vc3cmpt,
      init.params.out.all$init.params.multi$Vp3cmpt,
      init.params.out.all$init.params.multi$Vp23cmpt,
      init.params.out.all$init.params.multi$Q3cmpt,
      init.params.out.all$init.params.multi$Q23cmpt
    )
  )


  if (exists("all.out.vmax.km")) {
    colnames(all.out.vmax.km) <-
      c(
        "Method",
        "Estimated Vmax",
        "Estimated Km",
        "Absolute prediction error (APE)",
        "Mean absolute prediction error (MAPE)",
        "Time spent"
      )
  }
  if (exists("all.out.2cmpt")) {
    colnames(all.out.2cmpt) <-
      c(
        "Method",
        "Estimated CL",
        "Estimated Vp",
        "Estimated Q",
        "Absolute prediction error (APE)",
        "Mean absolute prediction error (MAPE)",
        "Time spent"
      )
  }
  if (exists("all.out.3cmpt")) {
    colnames(all.out.3cmpt) <-
      c(
        "Method",
        "Estimated CL",
        "Estimated Vp",
        "Estimated Vp2",
        "Absolute prediction error (APE)",
        "Estimated Q",
        "Estimated Q2",
        "Mean absolute prediction error (MAPE)",
        "Time spent"
      )
  }

  if (!exists("init.messages.multi")) {
    init.messages.multi <- NULL
  }
  if (!exists("init.messages.vmax.km")) {
    init.messages.vmax.km <- NULL
  }

  if (is.null(init.messages.vmax.km) &
      is.null(init.messages.vmax.km)) {
    init.messages <- NULL
  }

  if (is.null(init.messages.vmax.km) == F ||
      is.null(init.messages.vmax.km) == F) {
    init.messages <- c(init.messages.vmax.km, init.messages.multi)
  }

  if (!exists("npdcmpt.all.out")) {
    npdcmpt.all.out <-
      "No model fitting by naive pool data approach was conducted"
    npd.1cmpt_results=      "No model fitting by naive pool data approach was conducted"
    npd.1cmpt.mm_results=   "No model fitting by naive pool data approach was conducted"
    npd.2cmpt_results=      "No model fitting by naive pool data approach was conducted"
    npd.3cmpt_results=      "No model fitting by naive pool data approach was conducted"
    npd_1cmpt_out = "No model fitting by naive pool data approach was conducted"
    npd_1cmpt_mm_out = "No model fitting by naive pool data approach was conducted"
    npd_2cmpt_out = "No model fitting by naive pool data approach was conducted"
    npd_3cmpt_out= "No model fitting by naive pool data approach was conducted"

  }

  colnames(all.out) <-
    c(
      "Method",
      "Calculated Ka",
      "Calculated CL",
      "Calculated Vd",
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )
  colnames(sim.vmax.km.results.all) <-
    c(
      "Simulated Vmax",
      "Simulated Km",
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )
  colnames(sim.2cmpt.results.all) <-
    c(
      "Simulated Vc",
      "Simulated Vp",
      "Simulated Q",
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )
  colnames(sim.3cmpt.results.all) <-
    c(
      "Simulated Vc",
      "Simulated Vp",
      "Simulated Vp2",
      "Simulated Q",
      "Simulated Q2",
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )

  init.history <- list(
    base.out = all.out,
    sim.vmax.km = sim.vmax.km.results.all,
    sim.2cmpt = sim.2cmpt.results.all,
    sim.3cmpt = sim.3cmpt.results.all,
    npd_1cmpt_out =npd_1cmpt_out,
    npd_1cmpt_mm_out = npd_1cmpt_mm_out ,
    npd_2cmpt_out=npd_2cmpt_out,
    npd_3cmpt_out= npd_3cmpt_out
  )

  params.descriptions <- c(
    "CL: clearance",
    "Vd: volumn of distribution",
    "Vmax: maximum metobolic rate",
    "Km: Michaelis constant",
    "Vc: volume of distribution of the central compartment",
    "Vp: volume of distribution of the peripheral compartment",
    "Vp: volume of distribution of the second peripheral compartment",
    "Q: inter-compartmental clearance",
    "Q2: inter-compartmental clearance between central and second peripheral compartment"
  )

  output_env <- new.env()
  output_env$Datainfo <- Datainfo
  output_env$Recommended_initial_estimates <- Recommended_inits_df
  output_env$npdcmpt.all.out<-npdcmpt.all.out
  output_env$Run.history <- init.history
  output_env$Parameter.descriptions <- params.descriptions

  print_env_output(output_env)

  return(  output_env
    # list(
    #   Datainfo = Datainfo,
    #   Recommended_initial_estimates = Recommended_inits_df,
    #   Message =  init.messages.vmax.km,
    #   Run.history = init.history,
    #   Parameter.descriptions = params.descriptions
    # )
  )


} # end of function


#' Print environment summary for initial parameter estimation
#'
#' Prints a summary of initial parameter estimation stored in an environment.
#' It displays data information, recommended initial estimates, and parameter descriptions.
#' The function is designed to be used with results stored in a custom environment.
#'
#' @param env An environment containing the initial parameter estimation results.
#' The environment should include:
#'   \itemize{
#'     \item \code{Datainfo}: A string summarizing the data information (e.g., number of subjects, observations, etc.).
#'     \item \code{Recommended_initial_estimates}: A data frame containing recommended initial estimates for parameters.
#'     \item \code{Parameter_descriptions}: A character vector of descriptions for each parameter.
#'   }
#'
#' @return Prints a formatted summary to the console.
#' @examples
#' \dontrun{
#' # Assuming the environment 'env' contains the required fields:
#' print_env_output(env)
#' }
#' @export
#'
print_env_output <- function(env) {
  cat("===============Initial Parameter Estimation Summary ===============\n")
  cat("Data information:\n")
  print(env$Datainfo)

  cat("\nRecommended initial estimates by using non-model fitting methods:\n")
  print(head(env$Recommended_initial_estimates, 13))

  cat("\nRecommended initial estimates by naive pooled data compartmental analysis:\n")
  print(head(env$npdcmpt.all.out, 13))

  cat("\nParameter descriptions:\n")
  print(env$Parameter.descriptions)

  cat("\n=============== End of Summary ===============\n")
}

