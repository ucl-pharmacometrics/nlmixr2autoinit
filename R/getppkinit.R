#' Get initial estimates for a population pharmacokinetic modelling
#'
#' Calculates initial pharmacokinetic parameters using multiple methods including simplified calculation, non-compartmental analysis, graphical methods, and parameter estimation with naive pooled data approach if specified.
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param runnpd An integer indicating whether to run the naive pooled data approach (default is 0).
#' @param getinit.settings A list or data frame containing calculation settings (optional). It should include 'half_life', 'nlastpoints', 'nbins', and 'est.method'.
#' @return A list containing data information, initial parameter estimates, messages, and run history.
#' @importFrom dplyr %>% mutate filter select
#' @import nlmixr2
#' @importFrom tidyr fill
#' @import crayon
#' @examples
#' getppkinit(dat = Oral_1CPT,runnpd = 0)
#' getppkinit(dat = Bolus_1CPT,runnpd = 0)
#' getppkinit(dat = Bolus_1CPT,runnpd = 1)
#' getppkinit(dat = Bolus_1CPT,runnpd = 0, getinit.settings=c(trap.rule.method =2))
#' getppkinit(dat = Oral_1CPT,runnpd = 0)
#'
#' @export

getppkinit <- function(dat,
                       runnpd,
                       getinit.settings) {
  # do not run parameter estimation of model fitting
  if (missing(runnpd)) {
    runnpd = 0
  }
 # tool setting
  getinit.settings0 <- data.frame(
    half_life = NA,
    nlastpoints = 4,
    trap.rule.method=1,
    nbins = 8,
    est.method = "nls"
  )

  # If getinit.settings is missing, use the default settings
  if (missing(getinit.settings)) {
    getinit.settings <- getinit.settings0
  }

  else {
    getinit.settings <-
      as.data.frame(t(as.data.frame(getinit.settings)))
    # Update the values in getinit.settings0 (if they exist in getinit.settings)
    common_cols <-
      intersect(names(getinit.settings), names(getinit.settings0))
    getinit.settings0[common_cols] <- getinit.settings[common_cols]

  }

  message(black(
    paste0( "Settings of running nlmixr2autoinit",
            "\n")))

  message(blue(
      paste0(
            "Provided reference half-life------------------------------------------------------------------------------- ",
            getinit.settings0[1],
            "\n",
            "Number of plasma samples selected for linear regression on terminal phase slope---------------------------- ",
            getinit.settings0[2],
            "\n",
            "Trapezoidal rule method----------------------------------------------------------------------------------- ",
            getinit.settings0[3],
            "\n",
            "Number of bins during the naive pooled median data processing--------------------------------------------- ",
            getinit.settings0[4],
            "\n",
            "Estimated method for naive pooled approach data---------------------------------------------------------- ",
            getinit.settings0[5]
    )
  ))

  # Extract the final values
  half_life <- as.numeric(getinit.settings0$half_life)
  nlastpoints <- as.numeric(getinit.settings0$nlastpoints)
  trap.rule.method <- as.numeric(getinit.settings0$trap.rule.method)
  nbins <- as.numeric(getinit.settings0$nbins)
  est.method <- getinit.settings0$est.method
  ############################## Data information summary#######################

  # Calculate tad and dose number
  column_names <- toupper(colnames(dat))
  colnames(dat) <- toupper(colnames(dat))
  if ("ADDL" %in% column_names) {
    dat <- nmpkconvert(dat)
  }

  # Check whether infusion case, currently
  infusion_flag <- 0
  infusion_flag_c <- "N"

  if ("RATE" %in% column_names) {
    infusion_flag <- 1
    infusion_flag_c <- "Y"
  }

  if ("DUR" %in% column_names) {
    if (!"RATE" %in% column_names) {
      dat$RATE <- 0
      dat[dat$DUR > 0,]$RATE <-
        dat[dat$DUR > 0,]$AMT / dat[dat$DUR > 0,]$DUR
    }
  }

  dat <- calculate_tad(dat, infusion_flag)

  # check whether non intravenous case
  noniv_flag <- 0
  noniv_flag_c <- "N"

  if ("CMT" %in% column_names) {
    if (length(unique(dat$CMT)) > 1) {
      noniv_flag <- 1
      noniv_flag_c <- "Y"
      message(red(
        paste0(
          "Non-intravenous case(s) detected, Oral administration was assumed---------------------------------------- "
        )
      ))
    }
  }

  sdflag <- 0
  sdflag_c <- "N"
  # check whether single dose case
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
      ", Is infusion? ",
      infusion_flag_c,
      ", Is single dose? ",
      sdflag_c,
      ", First-dose data available? ",
      fdobsflag_c,

      ", is Oral case? ",
      noniv_flag_c
    )

  ########################Half-life estimated ##############################
  message(black(
    paste0("Performed linear regression on the terminal phase of pooled dose-normalized data, estimated half-life: ",
           half_life)))

   half_life<-half_life_estimated(dat = dat,
                                  sdflag = sdflag,
                                  fdobsflag = fdobsflag,
                                  nlastpoints = nlastpoints,
                                  nbins=nbins)

  ##############Rapid calculation############################################
  simpcal.APE <- Inf
  simpcal.MAPE <- Inf

  # iv case
  if (noniv_flag ==0){
    # half_life is estimated

    message(black(
      paste0("Run quick calculation with estimated half-life: ",
             half_life)))


    simpcal.results <- run_simpcal_iv(
    dat = dat,
    infusion_flag = infusion_flag,
    sdflag = sdflag,
    fdobsflag = fdobsflag,
    half_life = half_life)

    message(black(
      paste0("Run quick calculation without estimated half-life: ")))

    # use most commonly used dose-interval to replace half life
    simpcal.results.2 <- run_simpcal_iv(
      dat = dat,
      infusion_flag = infusion_flag,
      sdflag = sdflag,
      fdobsflag = fdobsflag,
      half_life = NA
  )

  simpcal.out <- simpcal.results$simpcal.results
  simpcal.out.2 <- simpcal.results.2$simpcal.results

 # Evaluate the provided initial values of cl and vd by their goodness of fit
  if (!is.na(simpcal.out$cl) & !is.na(simpcal.out$vd)) {
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
  if (!is.na(simpcal.out.2$cl) & !is.na(simpcal.out.2$vd)) {
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

  # oral case

    # message(black(
    #   paste0("Extravascular administration was observed; no rapid calculations were performed, and other processes continued.")))

  if (noniv_flag ==1){
    # half_life is estimated

    message(black(
      paste0("Run quick calculation with estimated half-life: ",
             half_life)))


    simpcal.results <- run_simpcal_iv(
      dat = dat,
      infusion_flag = infusion_flag,
      sdflag = sdflag,
      fdobsflag = 0, # do not run vd part
      half_life = half_life)

    message(black(
      paste0("Run quick calculation without estimated half-life: ")))

    # use most commonly used dose-interval to replace half life
    simpcal.results.2 <- run_simpcal_iv(
      dat = dat,
      infusion_flag = infusion_flag,
      sdflag = sdflag,
      fdobsflag = 0,
      half_life = NA
    )
    # obtained clearance
    simpcal.out <- simpcal.results$simpcal.results
    simpcal.out.2 <- simpcal.results.2$simpcal.results

  }

#################Non-compartmental analysis ##################################
 nca.results <- run_nca.normalised(
    dat = dat,
    nlastpoints = nlastpoints,
    trap.rule.method=trap.rule.method,
    nbins = nbins,
    fdobsflag = fdobsflag,
    sdflag=sdflag
  )

# run extra ka part
  if ( noniv_flag==1 ){
    # dat$DVnor <- dat$DV / dat$dose
    # dat_fd <- dat[dat$dose_number == 1, ]
    # datpooled_fd <- pk.time.binning(testdat = dat_fd,
    #                                 nbins = nbins)
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
    positive_ka_values <-    ka_values[   ka_values > 0 &  !is.na(ka_values)]
    ka_median <- round(median(positive_ka_values),2)

  }

  if (noniv_flag==0){
  # All pooled
  nca.APE <- Inf
  nca.MAPE <- Inf
  nca.results_all <- nca.results$nca.all.results
  nca.results_fd <- nca.results$nca.fd.results
  nca.results_efd <- nca.results$nca.efd.results

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

  nca_fd.APE <- Inf
  nca_fd.MAPE <- Inf
  if (!is.na(nca.results_fd$cl) & !is.na(nca.results_fd$vd)) {
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
  if (!is.na(nca.results_efd$cl) & !is.na(nca.results_efd$vd)) {
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
    nca.APE <- Inf
    nca.MAPE <- Inf
    nca.results_all <- nca.results$nca.all.results
    nca.results_fd <- nca.results$nca.fd.results
    nca.results_efd <- nca.results$nca.efd.results

    if (nca.results_all$cl > 0 & nca.results_all$vd > 0) {
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

    nca_fd.APE <- Inf
    nca_fd.MAPE <- Inf
    if (!is.na(nca.results_fd$cl) & !is.na(nca.results_fd$vd)) {
      if (nca.results_fd$cl > 0 & nca.results_fd$vd > 0) {
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
    if (!is.na(nca.results_efd$cl) & !is.na(nca.results_efd$vd)) {
      if (nca.results_efd$cl > 0 & nca.results_efd$vd > 0) {
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
if ( noniv_flag==0 ){
  graph.results_fd <- run_graphcal(
    dat = dat,
    fdobsflag =  fdobsflag,
    noniv_flag =  noniv_flag,
    nbins = nbins,
    nlastpoints = nlastpoints
  )

  graph_fd.APE <- Inf
  graph_fd.MAPE <- Inf
  if (!is.na(graph.results_fd$cl) & !is.na(graph.results_fd$vd)) {
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
if ( noniv_flag==1 ){
    graph.results_fd <- run_graphcal(
      dat = dat,
      fdobsflag =  fdobsflag,
      noniv_flag =  noniv_flag,
      nbins = nbins,
      nlastpoints = nlastpoints
    )

    graph_fd.APE <- Inf
    graph_fd.MAPE <- Inf
    if (!is.na(graph.results_fd$cl) & !is.na(graph.results_fd$vd)) {
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

  ############Hybrid simplified calculation##########################
  # volume of distribution was calculated by slope
  # Find an appropriate slope
  # if (nca_fd.MAPE < nca.MAPE) {
  #   usedslope <- nca.results_fd$slope
  # }
  #
  # if (nca_fd.MAPE >= nca.MAPE) {
  #   usedslope <- nca.results_all$slope
  # }

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
  if (!is.na(simpcal.out$cl)) {
    hybrid_cl <- simpcal.out$cl
    hybrid_vd <-  signif(-simpcal.out$cl / usedslope, 3)
  }

#iv case
 if ( noniv_flag==0 ){
  hybrid.APE <- Inf
  hybrid.MAPE <- Inf
  if (!is.na(hybrid_vd) & !is.na(simpcal.out$cl)) {
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
    if (!is.na(hybrid_vd) & !is.na(simpcal.out$cl)) {
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
      "Simplified calculation",
      "Graphic calculation",
      "NCA (only first dose)",
      "NCA (data exclude first-dose part)",
      "NCA (all pooled)",
      "Hybrid simplified calculation"
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
    ka=c(NA,graph.results_fd$ka,ka_method_1_fd,ka_median)
    all.out <- data.frame(
      method = c(
        "Simplified calculation",
        "Graphic calculation",
        "NCA (only first dose)",
        "Hybrid simplified calculation"
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
    round(all.out$`Mean absolute prediction error (MAPE)`, 0)
  # second check with APE if there are two rows
  # base.cl.best= as.numeric(min_mape_rows[min_mape_rows$`Mean absolute prediction error (MAPE)` ==
  #                                min(min_mape_rows$`Mean absolute prediction error (MAPE)`),]$`Calculated CL`)
  # base.vd.best= as.numeric(min_mape_rows[min_mape_rows$`Mean absolute prediction error (MAPE)` ==
  #                                min(min_mape_rows$`Mean absolute prediction error (MAPE)`),]$`Calculated Vd`)

  base.ka.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated Ka`

  base.cl.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated CL`

  base.vd.best <-
    all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`),]$`Calculated Vd`


message(black(
    paste0("Base parameter estimation finished. Estimated ka :",base.ka.best, ", estimated CL : ", base.cl.best, ", estimated Vd : ", base.vd.best )))

##############Naive pooled data approach for one-compartment model###############
# Start from one compartment to estimate cl and Vd.
  if (runnpd == 1) {
    if (noniv_flag==0){
    npd_1cmpt_out <- run_npd_1cmpt_iv(
      dat = dat,
      est.method = est.method,
      input.cl = mean(base.cl.best),
      input.vd = mean(base.vd.best)
    )
    }
    if (noniv_flag==1){
      npd_1cmpt_out <- run_npd_1cmpt_oral(
        dat = dat,
        est.method = est.method,
        input.ka = mean(base.ka.best),
        input.cl = mean(base.cl.best),
        input.vd = mean(base.vd.best)
      )
    }

    npd.1cmpt_results <- npd_1cmpt_out$npd.1cmpt_results
    npd.1cmpt.APE <- npd_1cmpt_out$npd.1cmpt.APE
    npd.1cmpt.MAPE <- npd_1cmpt_out$npd.1cmpt.MAPE

    all.out[(nrow(all.out) + 1),] <-
      c(
        "Naive pooled approach",
        npd.1cmpt_results$ka,
        npd.1cmpt_results$cl,
        npd.1cmpt_results$vd,
        npd.1cmpt.APE,
        npd.1cmpt.MAPE,
        npd.1cmpt_results$timespent
      )
    all.out$`Absolute Prediction Error (APE)` <-
      round(as.numeric(all.out$`Absolute Prediction Error (APE)`, 0))

    all.out$`Mean absolute prediction error (MAPE)` <-
      round(as.numeric(all.out$`Mean absolute prediction error (MAPE)`, 0))
  }


  ##############################Vmax and Km estimation##########################
  all.out.part <-
    all.out[all.out$Method != "Naive pooled approach",]
  all.out.part2 <- all.out

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

  #calculated cmax (c0) based on volume of distribution
  calc.cmax <- max.dose / min(base.vd.best)[1]
  fcmax <- max(c(mean.pop.cmax, calc.cmax))[1]
  if (fcmax == mean.pop.cmax) {
    cmax.message = "Cmax (mean maximum concontration) was obtained from mean observed cmax of analysis dataset"
  }

  if (fcmax == calc.cmax) {
    cmax.message = "Cmax was calculated from maximum dose divided by volume of distribution"
  }

  linear_minkm <-
    fcmax * 4 # if km>>4cmax, it nearly fall into the linear range

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


###############################wait############################################
  if (runnpd == 1) {
    # initial settings for cl and vd
    # use mean value if there are two options
    npdmm_inputcl <-
      mean(as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated CL`), na.rm = T)
    npdmm_inputvd <-
      mean(as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated Vd`), na.rm = T)

    npd_1cmpt_mm_out <- run_npd_1cmpt_mm_iv(
      dat = dat,
      est.method = est.method,
      npdmm_inputcl = npdmm_inputcl,
      npdmm_inputvd = npdmm_inputvd
    )

    npd.1cmpt.mm_results <- npd_1cmpt_mm_out$npd.1cmpt.mm_results
    npd.1cmpt.mm_results.2 <-
      npd_1cmpt_mm_out$npd.1cmpt.mm_results.2
    npd.1cmpt.mm.APE <- npd_1cmpt_mm_out$npd.1cmpt.mm.APE
    npd.1cmpt.mm.MAPE <- npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE
    npd.1cmpt.mm.APE.2 <- npd_1cmpt_mm_out$npd.1cmpt.mm.APE.2
    npd.1cmpt.mm.MAPE.2 <- npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE.2

    all.out.vmax.km <- data.frame(
      method = c(
        "Naive pooled approach (input_inits=1)",
        "Naive pooled approach (input_inits=threshold)"
      ),
      init.vmax = c(npd.1cmpt.mm_results$vmax, npd.1cmpt.mm_results.2$vmax),
      init.km = c(npd.1cmpt.mm_results$km, npd.1cmpt.mm_results.2$km),
      init.vd = c(npd.1cmpt.mm_results$vd, npd.1cmpt.mm_results.2$vd),
      modelAPE = c(npd.1cmpt.mm.APE, npd.1cmpt.mm.APE.2),
      modelMAPE = c(npd.1cmpt.mm.MAPE, npd.1cmpt.mm.MAPE.2),
      time.spent = c(
        npd.1cmpt.mm_results$timespent,
        npd.1cmpt.mm_results.2$timespent
      )
    )

    all.out.vmax.km$modelAPE <-
      round(as.numeric(all.out.vmax.km$modelAPE, 1))
    all.out.vmax.km$modelMAPE <-
      round(as.numeric(all.out.vmax.km$modelMAPE, 1))

  }

  ###########Multi-Compartmental Model Parameter Analysis#######################
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

  recommended_vc2cmpt_init <-
    sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE),]$vc[1]
  recommended_vp2cmpt_init <-
    sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE),]$vp[1]

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
                                              estka = base.vd.best[besti],
                                              noniv_flag = 1)
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results.all.i)
    rownames(sim.3cmpt.results.all) <-
      seq(1, nrow(sim.3cmpt.results.all), 1)
  }
 }

  recommended_vc3cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vc[1]
  recommended_vp3cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vp[1]
  recommended_vp23cmpt_init <-
    sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE),]$vp2[1]


  # Model fitting by naive pooled data approach
  if (runnpd == 1) {
    input.cl <- mean(base.cl.best, na.rm = T)
    input.vc2cmpt<- recommended_vc2cmpt_init
    input.vp2cmpt <- recommended_vp2cmpt_init
    input.q2cmpt <- input.cl

    npd_2cmpt_out <- run_npd_2cmpt_iv(
      dat = dat,
      est.method = est.method,
      input.cl=input.cl,
      input.vc2cmpt=  input.vc2cmpt,
      input.vp2cmpt= input.vp2cmpt,
      input.q2cmpt= input.q2cmpt,
    )
    npd.2cmpt_results <- npd_2cmpt_out$npd.2cmpt_results
    npd.2cmpt.APE <- npd_2cmpt_out$npd.2cmpt.APE
    npd.2cmpt.MAPE <- npd_2cmpt_out$npd.2cmpt.MAPE

    all.out.2cmpt <- data.frame(
      method = c("Naive pooled approach"),
      init.vc = c(npd.2cmpt_results$vc),
      init.vp = c(npd.2cmpt_results$vp),
      init.q = c(npd.2cmpt_results$q),
      modelAPE = c(npd.2cmpt.APE),
      modelMAPE = c(npd.2cmpt.MAPE),
      time.spent = c(npd.2cmpt_results$timespent)
    )

    all.out.2cmpt$modelAPE <-
      round(as.numeric(all.out.2cmpt$modelAPE, 1))

    all.out.2cmpt$modelMAPE <-
      round(as.numeric(all.out.2cmpt$modelMAPE, 1))

    input.vc3cmpt =   recommended_vc3cmpt_init
    input.vp3cmpt =   recommended_vp3cmpt_init
    input.vp23cmpt =  recommended_vp23cmpt_init
    input.q3cmpt =   input.cl

    npd_3cmpt_out <- run_npd_3cmpt_iv(dat = dat,
                                      est.method = est.method,
                                      input.cl=input.cl,
                                      input.vc3cmpt =  input.vc3cmpt,
                                      input.vp3cmpt =  input.vp3cmpt ,
                                      input.vp23cmpt = input.vp23cmpt,
                                      input.q3cmpt =   input.q3cmpt ,
                                      input.q23cmpt =  input.q23cmpt )

    npd.3cmpt_results <- npd_3cmpt_out$npd.3cmpt_results
    npd.3cmpt.APE <- npd_3cmpt_out$npd.3cmpt.APE
    npd.3cmpt.MAPE <- npd_3cmpt_out$npd.3cmpt.MAPE

    all.out.3cmpt <- data.frame(
      method = c("Naive pooled approach"),
      init.vc = c(npd.3cmpt_results$vc),
      init.vp = c(npd.3cmpt_results$vp),
      init.vp2 = c(npd.3cmpt_results$vp2),
      init.q = c(npd.3cmpt_results$q),
      init.q2 = c(npd.3cmpt_results$q2),
      modelAPE = c(npd.3cmpt.APE),
      modelMAPE = c(npd.3cmpt.MAPE),
      time.spent = c(npd.3cmpt_results$timespent)
    )

    all.out.3cmpt$modelAPE <-
      round(as.numeric(all.out.3cmpt$modelAPE, 1))

    all.out.3cmpt$modelMAPE <-
      round(as.numeric(all.out.3cmpt$modelMAPE, 1))
  }


  # Remove these temporary global variables run before.
  # List of variables to remove
  vars_to_remove <-
    c(
      "input.add",
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

  ######################## Finally selection########################################

  if (runnpd == 0) {
    # ka,cl,vd.
    f_init_ka <-
      as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                min(all.out.part$`Mean absolute prediction error (MAPE)`,
                                    na.rm = T),]$`Calculated Ka`)[1]
    f_init_cl <-
      as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                min(all.out.part$`Mean absolute prediction error (MAPE)`,
                                    na.rm = T),]$`Calculated CL`)[1]
    f_init_vd <-
      as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                min(all.out.part$`Mean absolute prediction error (MAPE)`,
                                    na.rm = T),]$`Calculated Vd`)[1]

    sel.method.ka.cl.vd <-
      all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` == min(all.out.part$`Mean absolute prediction error (MAPE)`,
                                                                               na.rm = T),]$Method[1]

    if (length(as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                       min(all.out.part$`Mean absolute prediction error (MAPE)`,
                                           na.rm = T),]$`Calculated CL`)) > 1) {
      #check the total volume of distribution of two compartment
      total.vd <-
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vc[1] +
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vp[1]


      f_init_ka <-
        as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                  min(all.out.part$`Mean absolute prediction error (MAPE)` ,
                                      na.rm = T) &
                                  round(all.out.part$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Ka`)[1]

      f_init_cl <-
        as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                  min(all.out.part$`Mean absolute prediction error (MAPE)` ,
                                      na.rm = T) &
                                  round(all.out.part$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated CL`)[1]

      f_init_vd <-
        as.numeric(all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` ==
                                  min(all.out.part$`Mean absolute prediction error (MAPE)` ,
                                      na.rm = T) &
                                  round(all.out.part$`Calculated Vd`, 1) == round(total.vd, 1),]$`Calculated Vd`)[1]

      sel.method.ka.cl.vd <-
        all.out.part[all.out.part$`Mean absolute prediction error (MAPE)` == min(all.out.part$`Mean absolute prediction error (MAPE)` ,
                                                                                 na.rm = T) &
                       round(all.out.part$`Calculated Vd`, 1) == round(total.vd, 1),]$Method[1]
    }


    # vmax.km
    f_init_vmax <-
      as.numeric(sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE ==
                                           min(sim.vmax.km.results.all$sim.mm.MAPE, na.rm = T),]$vmax)[1]
    f_init_km <-
      as.numeric(sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE ==
                                           min(sim.vmax.km.results.all$sim.mm.MAPE, na.rm = T),]$km)[1]
    sel.method.vmax.km <- "Sensitivity analysis by simulation "

    # multi-compartmental parameters
    f_init_vc2cmpt <-
      sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vc[1]
    f_init_vp2cmpt <-
      sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vp[1]

    f_init_vc3cmpt <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vc[1]
    f_init_vp3cmpt <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vp[1]
    f_init_vp23cmpt <-
      sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vp2[1]

    sel.method.multi <- "Sensitivity analysis by simulation "

    sel.method.ka<-"Wanger_nelson"

    if (sel.method.ka.cl.vd== "Hybrid simplified calculation"){
        sel.method.ka<-"Wanger_nelson_ (median)"
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
      vc3cmpt =  f_init_vc3cmpt,
      vp3cmpt =  f_init_vp3cmpt,
      vp23cmpt =   f_init_vp23cmpt
    )
  }

  # If naive pooled data approach is used
  if (runnpd == 1) {
    # cl selection

    f_init_cl <-
      as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                           min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated CL`)

    if (length(f_init_cl) == 1) {
      init.params.out.cl <-
        data.frame(method = all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                      min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$Method[1],
                   cl = f_init_cl)
    }

    if (length(f_init_cl) > 1) {
      f_init_cl <-
        mean(as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                  min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated CL`), na.rm = T)

      init.params.out.cl <-
        data.frame(
          method = paste0(
            "Mean values of method 1 (",
            all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                      min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$Method[1],
            ") and method 2",
            all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$Method[2],
            ")"
          ),
          cl = f_init_cl
        )

    }



    # Vd selection, identify whether vd should be selected from nonlinear model
    min.npdmmMAPE <- min(all.out.vmax.km$modelMAPE)[1]
    min.1cmptMAPE <-
      min(all.out.part$`Mean absolute prediction error (MAPE)`)[1]

    # Nnlinear model performs better, check values of vmax, km
    if (min.npdmmMAPE < min.1cmptMAPE) {
      # if fall into the linear range
      if (all.out.vmax.km[all.out.vmax.km$modelMAPE == min.npdmmMAPE,]$init.km[1] > linear_minkm) {
        f_init_vd <-
          as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                               min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated Vd`)

        if (length(f_init_vd) == 1) {
          sel.method <-
            all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                      min(all.out$`Mean absolute prediction error (MAPE)`,
                          na.rm = T),]$Method[1]

          init.params.out.vd <-
            data.frame(method = all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                          min(all.out$`Mean absolute prediction error (MAPE)`,
                                              na.rm = T),]$Method[1],
                       vd = f_init_vd)
        }

        if (length(f_init_vd) > 1) {
          f_init_vd <-
            mean(as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                                      min(all.out$`Mean absolute prediction error (MAPE)`,
                                          na.rm = T),]$`Calculated Vd`), na.rm = T)

          init.params.out.vd <-
            data.frame(
              method = paste0(
                "Mean values of method 1 (",
                all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                          min(all.out$`Mean absolute prediction error (MAPE)`,
                              na.rm = T),]$Method[1],
                ") and method 2",
                all.out[all.out$`Mean absolute prediction error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`,
                                                                               na.rm = T),]$Method[2],
                ")"
              ),
              vd = f_init_vd
            )
        }

      }

      else{
        f_init_vd <-
          all.out.vmax.km[all.out.vmax.km$modelMAPE == min.npdmmMAPE,]$init.vd[1]
        init.params.out.vd <-
          data.frame(method = "naive pooled data approach",
                     vd = f_init_vd)
      }
    }

    if (min.npdmmMAPE >= min.1cmptMAPE) {
      f_init_vd <-
        as.numeric(all.out[all.out$`Mean absolute prediction error (MAPE)` ==
                             min(all.out$`Mean absolute prediction error (MAPE)`, na.rm = T),]$`Calculated Vd`)[1]

      init.params.out.vd <-
        data.frame(method = "naive pooled data approach",
                   vd = f_init_vd)
    }


    # Check with simulation test
    npd_mm.MAPEmin <-
      min(npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE,
          npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE.2)
    sim.vmax.km.min <- min(sim.vmax.km.results.all$sim.mm.MAPE)

    # Vmax,km
    linear_minkm <- fcmax * 4
    linear_maxvmax <- linear_minkm * f_init_cl

    if (npd_mm.MAPEmin < sim.vmax.km.min) {
      if (all.out.vmax.km[all.out.vmax.km$modelMAPE == min.npdmmMAPE,]$init.km[1] > linear_minkm) {
        sel.method <- "Theshold value (km=4cmax,vmax=km*cl)"
        # vmax and km was recommended by simulation test
        f_init_vmax <- linear_maxvmax
        f_init_km <- linear_minkm

        init.params.out.vmax.km <- data.frame(method = sel.method,
                                              vmax = f_init_vmax,
                                              km =  f_init_km)

        init.messages.vmax.km <-
          c(
            "Estimated Km exceeds 4 times of the Cmax, it closely fall within the linear range. The threshold value (Km = 4*Cmax) is used.",
            cmax.message
          )
      }

      else{
        sel.method <- "naive pooled data approach"
        f_init_vmax <-
          all.out.vmax.km[all.out.vmax.km$modelMAPE == min(all.out.vmax.km$modelMAPE),]$init.vmax[1]
        f_init_km <-
          all.out.vmax.km[all.out.vmax.km$modelMAPE == min(all.out.vmax.km$modelMAPE),]$init.km[1]
        init.params.out.vmax.km <- data.frame(method = sel.method,
                                              vmax = f_init_vmax,
                                              km =  f_init_km)
      }
    }
    if (npd_mm.MAPEmin >= sim.vmax.km.min) {
      sel.method <- "Sensitivity analysis by simulation "
      # vmax and km was recommended by simulation test
      f_init_vmax <-
        sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE == sim.vmax.km.min,]$vmax[1]
      f_init_km <-
        sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE == sim.vmax.km.min,]$km[1]

      init.params.out.vmax.km <- data.frame(method = sel.method,
                                            vmax = f_init_vmax,
                                            km =  f_init_km)

      init.messages.vmax.km <-
        "Model with estimated Vmax and Km by naive pooled data approach has higher mean absolute prediction error (MAPE) than simulation test, the latter one was used for initial estimate recommendation "
    }


    # vc,vp,vp2
    # if two compartment does not fit well.

    if (npd_2cmpt_out$npd.2cmpt.MAPE >= min(sim.2cmpt.results.all$sim.2cmpt.MAPE)) {
      # still use the simulation test results
      f_init_vc2cmpt <-
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vc[1]
      f_init_vp2cmpt <-
        sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE == min(sim.2cmpt.results.all$sim.2cmpt.MAPE, na.rm = T),]$vp[1]

      f_init_vc3cmpt <-
        sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vc[1]
      f_init_vp3cmpt <-
        sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vp[1]
      f_init_vp23cmpt <-
        sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE == min(sim.3cmpt.results.all$sim.3cmpt.MAPE, na.rm = T),]$vp2[1]

      sel.method <- "Sensitivity analysis by simulation "

      init.params.out.vc.vp <- data.frame(
        method  =  sel.method,
        vc2cmpt =  f_init_vc2cmpt,
        vp2cmpt =   f_init_vp2cmpt,
        vc3cmpt =  f_init_vc3cmpt,
        vp3cmpt =  f_init_vp3cmpt,
        vp23cmpt =   f_init_vp23cmpt
      )

      init.messages.multi <-
        "Model with estimated mult-compartmental volume of distribution (vc,vpi) by naive pooled data approach shows higher mean absolute prediction error (MAPE) than simulation test, the latter one was used for initial estimate recommendation"

    }

    # if estimated parameters by npd is better than simulation test
    if (npd_2cmpt_out$npd.2cmpt.MAPE < min(sim.2cmpt.results.all$sim.2cmpt.MAPE)[1]) {
      f_init_vc2cmpt <- npd_2cmpt_out$npd.2cmpt_results$vc
      f_init_vp2cmpt <- npd_2cmpt_out$npd.2cmpt_results$vp
      f_init_vc3cmpt <- npd_3cmpt_out$npd.3cmpt_results$vc
      f_init_vp3cmpt <- npd_3cmpt_out$npd.3cmpt_results$vp
      f_init_vp23cmpt <- npd_3cmpt_out$npd.3cmpt_results$vp2


      init.params.out.vc.vp <- data.frame(
        method  =  "naive pooled data approach",
        vc2cmpt =  f_init_vc2cmpt,
        vp2cmpt =   f_init_vp2cmpt,
        vc3cmpt =  f_init_vc3cmpt,
        vp3cmpt =  f_init_vp3cmpt,
        vp23cmpt =   f_init_vp23cmpt
      )

    }
  }

  colnames(init.params.out.ka) <- c("Method", "Ka")
  colnames(init.params.out.cl) <- c("Method", "CL")
  colnames(init.params.out.vd) <- c("Method", "Vd")
  colnames(init.params.out.vmax.km) <- c("Method", "Vmax", "Km")
  colnames(init.params.out.vc.vp) <-
    c("Method",
      "Vc2cmpt",
      "Vp2cmpt",
      "Vc3cmpt",
      "Vp3cmpt",
      "Vp23cmpt")

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
      "Vc(3CMPT)",
      "Vp(3CMPT)",
      "Vp2(3CMPT)"
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
      init.params.out.all$init.params.multi$Vc3cmpt,
      init.params.out.all$init.params.multi$Vp3cmpt,
      init.params.out.all$init.params.multi$Vp23cmpt
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

  if (!exists("all.out.vmax.km")) {
    all.out.vmax.km <-
      "No model fitting by naive pool data approach was conducted"
  }

  if (!exists("all.out.2cmpt")) {
    all.out.2cmpt <-
      "No model fitting by naive pool data approach was conducted"
  }

  if (!exists("all.out.3cmpt")) {
    all.out.3cmpt <-
      "No model fitting by naive pool data approach was conducted"
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
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )
  colnames(sim.3cmpt.results.all) <-
    c(
      "Simulated Vc",
      "Simulated Vp",
      "Simulated Vp2",
      "Absolute prediction error (APE)",
      "Mean absolute prediction error (MAPE)",
      "Time spent"
    )

  init.history <- list(
    base.out = all.out,
    sim.vmax.km = sim.vmax.km.results.all,
    sim.2cmpt = sim.2cmpt.results.all,
    sim.3cmpt = sim.3cmpt.results.all,
    npd.out.vmax.km = all.out.vmax.km,
    npd.out.2cmpt = all.out.2cmpt,
    npd.out.3cmpt = all.out.3cmpt
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
  return(
    list(
      Datainfo = Datainfo,
      Recommended_initial_estimates = Recommended_inits_df,
      Message =  init.messages.vmax.km,
      Run.history = init.history,
      Parameter.descriptions = params.descriptions
    )
  )


} # end of function
