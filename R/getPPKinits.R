#' Get initial estimates for a population pharmacokinetic modelling
#'
#' Computes initial values of pharmacokinetic parameters using integrated pipeline which includes single-point method, non-compartmental analysis, graphical methods, simulation-based analysis and parameter estimation with naive pooled data approach compartmental analysis if specified.
#' @param dat A data frame containing the pharmacokinetic data.
#' @param run.option Integer value indicating which methods to use. It has three possible values:
#' \itemize{
#'   \item \code{0}: Default pipeline calculation, not involving parameter estimation or model fitting.
#'   \item \code{1}: Applies both the pipeline calculation methods and naive pooled data compartmental analysis.
#'   \item \code{2}: Uses only naive pooled data compartmental analysis.
#' }
#' @param getinit.settings A list or data frame containing calculation settings (optional). The following settings can be provided:
#' \itemize{
#'   \item \code{half_life}: Numeric value for the half-life of the drug (default is \code{NA}). If not provided, it will be estimated based on the data.
#'   \item \code{nlastpoints}: Numeric value specifying the number of last data points used for linear regression to obtain the slope in the terminal phase. (default is \code{4}).
#'   \item \code{trapezoidal.rule}: Numeric value indicating the trapezoidal rule method to use (default is \code{0}).
#'   \code{1} refers to the linear trapezoidal method, while \code{2} refers to the linear-up/log-down trapezoidal method.
#'   \item \code{nbins}: Numeric value specifying the number of time windows for quantile-based partitioning of the time variable in the dataset when performing naive pooling of data (default is \code{8}).
#'   \item \code{est.method}: Character string indicating the estimation method for naive pooled data compartmental analysis to use (default is \code{"nls"}) or other methods ((e.g., "nls", "nlm", "nlminb","foce", "focei") depending on the analysis.
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
#' @importFrom knitr kable
#' @examples
#' inits.out<-getPPKinits(dat = Infusion_1CPT,run.option = 2,getinitsControl = initsControl(est.method = "nls"))
#' inits.out
#'
#'
#'
#'
#' @export
#'

getPPKinits<- function(dat,
                       run.option=0,
                       getinitsControl=initsControl()) {

  # Get function settings
  half_life <-  getinitsControl$half_life
  nlastpoints <- getinitsControl$nlastpoints
  trapezoidal.rule <- getinitsControl$trapezoidal.rule
  nbins <-getinitsControl$nbins
  est.method <-getinitsControl$ est.method
  selection.criteria<-getinitsControl$ selection.criteria
  npdcmpt.inits.strategy<-getinitsControl$npdcmpt.inits.strategy

  function.params <- data.frame(
    Parameter = c(
      "Run option (0= only run pipeline, 1= only run NPD compartmental analysis, 3= run both)",
      "Half-life (User-defined half-life used as a reference for single-point calculation)",
      "Number of last points (Points used for linear regression during terminal elimination phase slope estimation)",
      "Number of bins (Number of time windows derived from quantile-based partitioning of the time variable in the dataset)",
      "Trapezoidal rule for AUC (1=linear,2=Linear-up log-down )",
      "Method used for naive pooled data compartmental analysis",
      "Selection criteria used for evaluating and selecting parameter values",
      "NPD compartmental analysis initial setting strategy (0 = set as 1, 1 = set from pipeline recommended values) "
    ),
    Value = c(
      run.option,
      half_life,
      nlastpoints,
      nbins,
      trapezoidal.rule,
      est.method,
      selection.criteria,
      npdcmpt.inits.strategy
    )
  )

  colnames(function.params)<-c("Settings","Values")
  # Display the initials setting
  message(magenta( paste(capture.output(knitr::kable(function.params, format = "simple")), collapse = "\n")))
  message(magenta(paste0("---------------------------------------------------------------------------------------------------------------------  -------")))
  ################## Data processing and information summary#################
  message(black(
    paste0("Processing data", strrep(".", 20))
  ))

  # fdat<-processData(dat = Bolus_1CPT)
  processData.out<-processData(dat = dat)
  dat<-  processData.out$dat
  Datainfo<-  processData.out$Datainfo

  # Reset ID
  dat$ID<-dat$ID+(dat$resetflag-1)*max(dat$ID)

  # prepare flag for analysis
    fdobsflag <-0
    mdobsflag<-0
    sdflag <-0
    bolus_flag<-0
    infusion_flag<-0
    oral_flag<-0

  if (nrow(dat[dat$dose_number==1 & dat$EVID==0 & dat$iiobs==0,])>0){
    fdobsflag = 1
  }

  if (nrow(dat[dat$dose_number>1 & dat$EVID==0,])>0){
    mdobsflag  =1
  }

   if (nrow(dat[dat$dose_number>1 & dat$EVID==0,])==0){
     sdflag <-1
   }

  if (unique(dat$routeobs)=="bolus"){
    bolus_flag<-1
    route<-"bolus"
  }

  if (unique(dat$routeobs)=="infusion"){
   infusion_flag<-1
   route<-"infusion"
  }

  if (unique(dat$routeobs)=="oral"){
    oral_flag<-1
    route<-"oral"
  }

  rawdat<-dat # keep raw dat for naive pooled data compartmental analysis

########################### Pipeline part ##################################
  if (run.option<2){
######################## Half-life estimated ################################

   half_life_out<-half_life_estimated(dat = dat,
                                  nlastpoints = nlastpoints,
                                  nbins=nbins,
                                  route=route)

   half_life<-  half_life_out$half_life_median

   message(black(
     paste0("Estimated half-life : ", half_life)))


  ####################Single point method ################################

  message(black(
     paste0("Run single-point method to calculate PK parameters",strrep(".", 20))))

 single.point.lst <- run_single_point(
    dat = dat,
    half_life = half_life)

 single.point.out <- single.point.lst$singlepoint.results
 dat<-single.point.lst$dat

################# Naive pooled Non-compartmental analysis ##################################
  message(black(
    paste0("Run non-compartmental analysis on naive pooling data", strrep(".", 20))
  ))

 nca.results <- run_nca.normalised(
    dat = dat,
    nlastpoints = nlastpoints,
    trapezoidal.rule = trapezoidal.rule,
    nbins = nbins,
    fdobsflag = fdobsflag,
    sdflag=sdflag,
    route=route
  )

 # Determine absorption rate of constant
 # Wanger nelson method (needs more than one data in the absorption phase)
    ka_method_1_fd=NA
    ka_method_1_out_fd=NA
    ka_method_1_efd=NA
    ka_method_1_out_efd=NA
    ka_method_1_all=NA
    ka_method_1_out_all=NA
    ka_wanger_nelson_result="No ka calculation using the Wagner-Nelson method was performed."

  if (oral_flag==1 & length(nca.results$datpooled_fd)==2 ){
     if (is.na(nca.results$nca.fd.results$cl)==F & is.na(nca.results$nca.fd.results$vd)==F){
     ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_fd$test.pool.normalised,
                      nlastpoints = nlastpoints,
                      nca.out = unlist(nca.results$nca.fd.results, use.names = FALSE))

     ka_method_1_fd <-signif(ka_wanger_nelson_result$ka,3)
     ka_method_1_out_fd<-ka_wanger_nelson_result$dat_out_wanger_nelson
     }
  }

  # Not applicable for multiple doses.
  # if (oral_flag==1 & length(nca.results$datpooled_efd)==2){
  #   if (is.na(nca.results$nca.efd.results$cl)==F & is.na(nca.results$nca.efd.results$vd)==F){
  #     ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_efd$test.pool.normalised,
  #                                               nlastpoints = nlastpoints,
  #                                               nca.out = unlist(nca.results$nca.efd.results, use.names = FALSE))
  #     ka_method_1_efd <-signif(ka_wanger_nelson_result$ka,3)
  #     ka_method_1_out_efd<-ka_wanger_nelson_result$dat_out_wanger_nelson
  #   }
  #   }

  # if (oral_flag==1 & length(nca.results$datpooled_all)==2){
  #     if (is.na(nca.results$nca.all.results$cl)==F & is.na(nca.results$nca.all.results$vd)==F){
  #     ka_wanger_nelson_result<-ka_wanger_nelson(dat = nca.results$datpooled_all$test.pool.normalised,
  #                                               nlastpoints = nlastpoints,
  #                                               nca.out = unlist(nca.results$nca.all.results, use.names = FALSE))
  #     ka_method_1_all <-signif(ka_wanger_nelson_result$ka,3)
  #     ka_method_1_out_all<-ka_wanger_nelson_result$dat_out_wanger_nelson
  #     }
  #  }

    # can be used for later hybrid method
    ka_values<-c(ka_method_1_fd, ka_method_1_efd, ka_method_1_all)
    # Remove negative numbers
    positive_ka_values <-  ka_values[   ka_values > 0 &  is.na(ka_values)==F]
    ka_median <- round(median(positive_ka_values),2)

 ############################ Graphic analysis#################################
   message(black(
     paste0("Run graphical analysis on naive pooling data only after the first dose", strrep(".", 20))
   ))

   # infusion case, the infusion rate should be very fast
   # if (route=="infusion"){
   #    duration_value_counts<- table(dat[dat$EVID %in% c(1,4) &!duplicated(dat$ID), ]$AMT
    # /dat[dat$EVID %in% c(1,4)&!duplicated(dat$ID), ]$RATE)
   #    most_commonly_used_infusion_duration<- as.numeric(names(which.max(duration_value_counts)))
   #  }
   graph_fd.APE <- NA
   graph_fd.MAE <- NA
   graph_fd.MAPE <- NA
   graph_fd.RMSE <- NA
   graph_fd.rRMSE <- NA

   graph.results_fd <- run_graphcal(
     dat = dat,
     route =   route ,
     nbins = nbins,
     nlastpoints = nlastpoints,
     fdobsflag=fdobsflag
   )


######### Predictive performance evaluation (one-compartmental parameter)#######
# default setting
   simpcal.APE <- NA
   simpcal.MAE <- NA
   simpcal.MAPE <- NA
   simpcal.RMSE <- NA
   simpcal.rRMSE <- NA

   nca.APE <- NA
   nca.MAE <- NA
   nca.MAPE <- NA
   nca.RMSE <- NA
   nca.rRMSE <- NA

   nca_fd.APE <- NA
   nca_fd.MAE <- NA
   nca_fd.MAPE <- NA
   nca_fd.RMSE <- NA
   nca_fd.rRMSE <- NA

   nca_efd.APE <- NA
   nca_efd.MAE <- NA
   nca_efd.MAPE <- NA
   nca_efd.RMSE <- NA
   nca_efd.rRMSE <- NA

   graph_fd.APE <- NA
   graph_fd.MAE <- NA
   graph_fd.MAPE <- NA
   graph_fd.RMSE <- NA
   graph_fd.rRMSE <- NA

#============================== for iv case===================================#
  if (bolus_flag==1 || infusion_flag==1){
    ka=  c(NA,NA,NA,NA,NA)
   # single-point method
    if (is.na(single.point.out$cl)==F & is.na(single.point.out$vd)==F) {
      simpcal_sim <-
        Fit_1cmpt_iv(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.cl = single.point.out$cl,
          input.vd = single.point.out$vd,
          input.add = 0
        )

      simpcal.APE <-  round(metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
      simpcal.MAE <-  round(metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
      simpcal.MAPE <- round( metrics.(pred.x = simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
      simpcal.RMSE <-  round(metrics.(pred.x = simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
      simpcal.rRMSE <- round( metrics.(pred.x = simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

      rm(simpcal_sim)
      gc()
    }


  # non-compartmental analysis
  nca.results_all <- nca.results$nca.all.results
  nca.results_fd <- nca.results$nca.fd.results
  nca.results_efd <- nca.results$nca.efd.results

  if (is.na(nca.results_all$cl)==F & is.na(nca.results_all$vd )==F) {
  if (nca.results_all$cl > 0 & nca.results_all$vd > 0) {
    nca_sim <- Fit_1cmpt_iv(
      data = dat[dat$EVID != 2,],
      est.method = "rxSolve",
      input.cl = nca.results_all$cl,
      input.vd = nca.results_all$vd,
      input.add = 0
    )
    nca.APE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
    nca.MAE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
    nca.MAPE <- round( metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
    nca.RMSE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
    nca.rRMSE <- round( metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

    rm(nca_sim)
    gc()
  }
  }

  if (is.na(nca.results_fd$cl)==F & is.na(nca.results_fd$vd)==F) {
    if (nca.results_fd$cl > 0 & nca.results_fd$vd > 0) {
      nca_fd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = nca.results_fd$cl,
        input.vd = nca.results_fd$vd,
        input.add = 0
      )

      nca_fd.APE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
      nca_fd.MAE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
      nca_fd.MAPE <- round( metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
      nca_fd.RMSE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
      nca_fd.rRMSE <- round( metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

      rm(nca_fd_sim)
      gc()


    }
  }

  if (is.na(nca.results_efd$cl)==F & is.na(nca.results_efd$vd)==F) {
    if (nca.results_efd$cl > 0 & nca.results_efd$vd > 0) {
      nca_efd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = nca.results_efd$cl,
        input.vd = nca.results_efd$vd,
        input.add = 0
      )

      nca_efd.APE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
      nca_efd.MAE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
      nca_efd.MAPE <- round( metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
      nca_efd.RMSE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
      nca_efd.rRMSE <- round( metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

      rm( nca_efd_sim )
      gc()


    }
  }

  # graphic analysis
  if (is.na(graph.results_fd$cl)==F & is.na(graph.results_fd$vd)==F) {
    if (graph.results_fd$cl > 0 & graph.results_fd$vd > 0) {
      graph_fd_sim <- Fit_1cmpt_iv(
        data = dat[dat$EVID != 2,],
        est.method = "rxSolve",
        input.cl = graph.results_fd$cl,
        input.vd = graph.results_fd$vd,
        input.add = 0
      )

      graph_fd.APE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
      graph_fd.MAE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
      graph_fd.MAPE <- round( metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
      graph_fd.RMSE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
      graph_fd.rRMSE <- round( metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

      rm( graph_fd_sim)
      gc()

    }
  }
}

#================================== for oral case================================#
  if (oral_flag==1){
    ka=  c(single.point.out$ka,graph.results_fd$ka,ka_method_1_fd,ka_method_1_efd,ka_method_1_all)
    # single-point method
    if (is.na(single.point.out$cl)==F & is.na(single.point.out$vd)==F & is.na(  single.point.out$ka)==F) {
      if (single.point.out$cl > 0 & single.point.out$vd > 0 & single.point.out$ka>0 ) {
        simpcal_sim <- Fit_1cmpt_oral(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.ka = single.point.out$ka,
          input.cl = single.point.out$cl,
          input.vd = single.point.out$vd,
          input.add = 0
        )

        simpcal.APE <-  round(metrics.(pred.x = simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
        simpcal.MAE <-  round(metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
        simpcal.MAPE <- round( metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
        simpcal.RMSE <-  round(metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
        simpcal.rRMSE <- round( metrics.(pred.x =  simpcal_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

        rm(simpcal_sim)
        gc()

      }
    }
    # non-compartmental analysis
    nca.results_all <- nca.results$nca.all.results
    nca.results_fd <- nca.results$nca.fd.results
    nca.results_efd <- nca.results$nca.efd.results

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

      nca.APE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
      nca.MAE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
      nca.MAPE <- round( metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
      nca.RMSE <-  round(metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
      nca.rRMSE <- round( metrics.(pred.x = nca_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

      rm(nca_sim)
      gc()

    }
    }

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
        nca_fd.APE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
        nca_fd.MAE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
        nca_fd.MAPE <- round( metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
        nca_fd.RMSE <-  round(metrics.(pred.x = nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
        nca_fd.rRMSE <- round( metrics.(pred.x =nca_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

        rm(nca_fd_sim)
        gc()

      }
    }

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
        nca_efd.APE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
        nca_efd.MAE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
        nca_efd.MAPE <- round( metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
        nca_efd.RMSE <-  round(metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
        nca_efd.rRMSE <- round( metrics.(pred.x = nca_efd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

        rm(nca_efd_sim)
        gc()
      }
    }

    # Graphic analysis
    if (is.na(graph.results_fd$cl)==F & is.na(graph.results_fd$vd)==F & is.na(graph.results_fd$ka)==F) {
      if (graph.results_fd$cl > 0 & graph.results_fd$vd > 0 & graph.results_fd$ka > 0) {
        graph_fd_sim <- Fit_1cmpt_oral(
          data = dat[dat$EVID != 2,],
          est.method = "rxSolve",
          input.ka = graph.results_fd$ka,
          input.cl = graph.results_fd$cl,
          input.vd = graph.results_fd$vd,
          input.add = 0
        )
        graph_fd.APE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[1],3)
        graph_fd.MAE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[2],3)
        graph_fd.MAPE <- round( metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[3],3)
        graph_fd.RMSE <-  round(metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[4],3)
        graph_fd.rRMSE <- round( metrics.(pred.x = graph_fd_sim$cp ,obs.y =dat[dat$EVID == 0,]$DV )[5],3)

        rm( graph_fd_sim)
        gc()

      }
    }

  }

##################### Summarise one-compartmental parameter results################


all.out <- data.frame(
    method = c(
      "Single-point method",
      "Graphic analysis",
      "NCA (only first dose)",
      "NCA (data exclude first-dose part)",
      "NCA (all pooled)"
    ),

    ka=ka,

    cl = c(
      single.point.out$cl,
      graph.results_fd$cl,
      nca.results_fd$cl,
      nca.results_efd$cl,
      nca.results_all$cl
    ),

    vd = c(
      single.point.out$vd,
      graph.results_fd$vd,
      nca.results_fd$vd,
      nca.results_efd$vd,
      nca.results_all$vd
    ),

    simAPE = c(simpcal.APE,
               graph_fd.APE,
               nca_fd.APE,
               nca_efd.APE,
               nca.APE),

    simMAE = c(
      simpcal.MAE,
      graph_fd.MAE,
      nca_fd.MAE,
      nca_efd.MAE,
      nca.MAE),

    simMAPE = c(simpcal.MAPE,
                graph_fd.MAPE,
                nca_fd.MAPE,
                nca_efd.MAPE,
                nca.MAPE),

    simRMSE = c(
      simpcal.RMSE,
      graph_fd.RMSE,
      nca_fd.RMSE,
      nca_efd.RMSE,
      nca.RMSE
    ),

    simrRMSE = c(simpcal.rRMSE,
                 graph_fd.rRMSE,
                 nca_fd.rRMSE,
                 nca_efd.rRMSE,
                 nca.rRMSE),

    time.spent = c(
      single.point.out$time.spent,
      graph.results_fd$time.spent,
      nca.results_fd$time.spent,
      nca.results_efd$time.spent,
      nca.results_all$time.spent
    )
  )

  colnames(all.out) <-
    c(
      "Method",
      "Calculated Ka",
      "Calculated CL",
      "Calculated Vd",
      "Absolute Predicted Error (APE)",
      "Mean Absolute Error (MAE)",
      "Mean Absolute Percentage Error (MAPE)",
      "Root Mean Squared Error (RMSE)",
      "Relative Root Mean Squared Error (rRMSE)",
      "Time spent"
    )


  if (selection.criteria=="APE"){
  base.best <-
    all.out[all.out$`Absolute Predicted Error (APE)` == min(all.out$`Absolute Predicted Error (APE)`,na.rm = T) & is.na(all.out$`Absolute Predicted Error (APE)`)==F,]
  }

  if (selection.criteria=="MAE"){
  base.best <-
    all.out[all.out$`Mean Absolute Error (MAE)` == min(all.out$`Mean Absolute Error (MAE)`,na.rm = T) & is.na(all.out$`Mean Absolute Error (MAE)`)==F ,]
  }

  if (selection.criteria=="MAPE"){
    base.best <-
      all.out[all.out$`Mean Absolute Percentage Error (MAPE)` == min(all.out$`Mean absolute prediction error (MAPE)`,na.rm = T) & is.na(all.out$`Mean absolute prediction error (MAPE)`)==F,]
  }

  if (selection.criteria=="RMSE"){
    base.best <-
      all.out[all.out$`Root Mean Squared Error (RMSE)` == min(all.out$`Root Mean Squared Error (RMSE)`,na.rm = T) & is.na(all.out$`Mean absolute prediction error (MAPE)`)==F ,]
  }

  if (selection.criteria=="rRMSE"){
    base.best <-
      all.out[all.out$`Relative Root Mean Squared Error (rRMSE)` == min(all.out$`Relative Root Mean Squared Error (rRMSE)`,na.rm = T),]
  }

# Only use the first record (if same, temporary setting)
  if (nrow(base.best)>1){
    base.best<-base.best[1,]
  }

  base.ka.best<-base.best$`Calculated Ka`
  base.cl.best <-base.best$`Calculated CL`
  base.vd.best <-base.best$`Calculated Vd`

message_text <- paste0("Base PK parameter analysis finished. Estimated ka: ", base.ka.best,
                       ", estimated CL: ", base.cl.best, ", estimated Vd: ", base.vd.best)

cat(message_text, "\n")


################# Simulation-based Vmax and Km analysis#######################

  message(black(
    paste0("Run simulation-based sensitivity analysis on nonlinear eliminiation kinetics PK parameters",strrep(".", 20))))
  # max.dose <-
  #     max(dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0, ]$dose)[1]

  # obtain cmax in the current dataset
  dat.obs <- dat[dat$EVID == 0,]
  pop.cmax <- aggregate(dat.obs$DV,
                        list(dat.obs$ID),
                        FUN = max,
                        na.rm = T)
  mean.pop.cmax <- mean(pop.cmax$x, na.rm = T)
  fcmax<-  mean.pop.cmax
 #  Make an assumption that when the parameter Km (Michaelis constant) is much greater than 4 * Cmax, drug metabolism largely enters the linear range.
  linear_minkm <- fcmax * 4 # if km>>4cmax, it nearly fall into the linear range

  sim.vmax.km.results.all <- NULL

  if (bolus_flag==1 || infusion_flag==1){
    sim.vmax.km.results.all.i <- sim_sens_vmax_km(
      dat = dat,
      estcmax =  fcmax,
      estcl = base.cl.best,
      estvd = base.vd.best
    )
    sim.vmax.km.results.all <-
      rbind(sim.vmax.km.results.all, sim.vmax.km.results.all.i)
    rownames(sim.vmax.km.results.all) <-
      seq(1, nrow(sim.vmax.km.results.all), 1)
   }

  if (oral_flag==1){
      sim.vmax.km.results.all.i <- sim_sens_vmax_km(
        dat = dat,
        estcmax =  fcmax,
        estcl = base.cl.best,
        estvd = base.vd.best,
        estka = base.ka.best,
        noniv_flag = 1
      )
      sim.vmax.km.results.all <-
        rbind(sim.vmax.km.results.all, sim.vmax.km.results.all.i)
      rownames(sim.vmax.km.results.all) <-
        seq(1, nrow(sim.vmax.km.results.all), 1)
    }

  # message(black(
  #   paste0("Nonlinear elimination parameter analysis finished. Estimated Vmax : ",   recommended_vmax_init, ", estimated km : ", recommended_km_init  )))

  ########### Simulation-based Multi-Compartmental Model Parameter Analysis#####
  message(black(
    paste0("Run simulation-based sensitivity analysis on multi-compartmental PK parameters",strrep(".", 20))))

  # Collect identified vc from single-point extra and base.best.vd

  # Two-compartment model simulation
  sim.2cmpt.results.all <- NULL

  approx.vc.value<-single.point.lst$approx.vc.out$approx.vc.value

  if (bolus_flag==1 || infusion_flag==1){
    sim.2cmpt.results.all.i <- sim_sens_2cmpt(dat = dat,
                                              estcl = base.cl.best,
                                              sim_vc_list = c( approx.vc.value, base.vd.best))
    sim.2cmpt.results.all <-
      rbind(sim.2cmpt.results.all, sim.2cmpt.results.all.i)
    rownames(sim.2cmpt.results.all) <-
      seq(1, nrow(sim.2cmpt.results.all), 1)
   }

  if (oral_flag==1){
      sim.2cmpt.results.all.i <- sim_sens_2cmpt(dat = dat,
                                                estcl = base.cl.best,
                                                sim_vc_list = c( approx.vc.value, base.vd.best),
                                                estka = base.ka.best,
                                                noniv_flag = 1)
      sim.2cmpt.results.all <-
        rbind(sim.2cmpt.results.all, sim.2cmpt.results.all.i)
      rownames(sim.2cmpt.results.all) <-
        seq(1, nrow(sim.2cmpt.results.all), 1)

  }


  # Three-compartment model simulation
  sim.3cmpt.results.all <- NULL

  if (bolus_flag==1 || infusion_flag==1){

    sim.3cmpt.results.all.i <- sim_sens_3cmpt(dat = dat,
                                              estcl = base.cl.best,
                                              sim_vc_list = c( approx.vc.value, base.vd.best))
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results.all.i)
    rownames(sim.3cmpt.results.all) <-
      seq(1, nrow(sim.3cmpt.results.all), 1)
  }


  if (oral_flag==1){

    sim.3cmpt.results.all.i <- sim_sens_3cmpt(dat = dat,
                                              estcl = base.cl.best,
                                              sim_vc_list = c( approx.vc.value, base.vd.best),
                                              estka = base.ka.best,
                                              noniv_flag = 1)
    sim.3cmpt.results.all <-
      rbind(sim.3cmpt.results.all, sim.3cmpt.results.all.i)
    rownames(sim.3cmpt.results.all) <-
      seq(1, nrow(sim.3cmpt.results.all), 1)
  }

################## Parameter Selection Selection############################
  if (selection.criteria=="APE"){
    recommended_mm<-
      sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.APE == min(sim.vmax.km.results.all$sim.mm.APE,na.rm = T) & is.na(sim.vmax.km.results.all$sim.mm.APE)==F,][1,]

    recommended.multi1<-sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.APE== min(sim.2cmpt.results.all$sim.2cmpt.APE,na.rm = T) & is.na(sim.2cmpt.results.all$sim.2cmpt.APE)==F,][1,]
    recommended.multi2<-sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.APE== min(sim.3cmpt.results.all$sim.3cmpt.APE,na.rm = T) & is.na(sim.3cmpt.results.all$sim.3cmpt.APE)==F,][1,]
  }

  if (selection.criteria=="MAE"){
    recommended_mm<-
      sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAE == min(sim.vmax.km.results.all$sim.mm.MAE,na.rm = T) & is.na(sim.vmax.km.results.all$sim.mm.MAE)==F,][1,]

    recommended.multi1<-sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAE== min(sim.2cmpt.results.all$sim.2cmpt.MAE,na.rm = T) & is.na(sim.2cmpt.results.all$sim.2cmpt.MAE)==F,][1,]
    recommended.multi2<-sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAE== min(sim.3cmpt.results.all$sim.3cmpt.MAE,na.rm = T) & is.na(sim.3cmpt.results.all$sim.3cmpt.MAE)==F,][1,]
  }

  if (selection.criteria=="MAPE"){
    recommended_mm<-
      sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.MAPE == min(sim.vmax.km.results.all$sim.mm.MAPE,na.rm = T) & is.na(sim.vmax.km.results.all$sim.mm.MAPE)==F,][1,]

    recommended.multi1<-sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.MAPE== min(sim.2cmpt.results.all$sim.2cmpt.MAPE,na.rm = T) & is.na(sim.2cmpt.results.all$sim.2cmpt.MAPE)==F,][1,]
    recommended.multi2<-sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.MAPE== min(sim.3cmpt.results.all$sim.3cmpt.MAPE,na.rm = T) & is.na(sim.3cmpt.results.all$sim.3cmpt.MAPE)==F,][1,]
  }

  if (selection.criteria=="RMSE"){
    recommended_mm<-
      sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.RMSE == min(sim.vmax.km.results.all$sim.mm.RMSE,na.rm = T) & is.na(sim.vmax.km.results.all$sim.mm.RMSE)==F,][1,]

    recommended.multi1<-sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.RMSE== min(sim.2cmpt.results.all$sim.2cmpt.RMSE,na.rm = T) & is.na(sim.2cmpt.results.all$sim.2cmpt.RMSE)==F,][1,]
    recommended.multi2<-sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.RMSE== min(sim.3cmpt.results.all$sim.3cmpt.RMSE,na.rm = T) & is.na(sim.3cmpt.results.all$sim.3cmpt.RMSE)==F,][1,]
  }

  if (selection.criteria=="rRMSE"){
    recommended_mm<-
      sim.vmax.km.results.all[sim.vmax.km.results.all$sim.mm.rRMSE == min(sim.vmax.km.results.all$sim.mm.rRMSE,na.rm = T) & is.na(sim.vmax.km.results.all$sim.mm.rRMSE)==F,][1,]

    recommended.multi1<-sim.2cmpt.results.all[sim.2cmpt.results.all$sim.2cmpt.rRMSE== min(sim.2cmpt.results.all$sim.2cmpt.rRMSE,na.rm = T) & is.na(sim.2cmpt.results.all$sim.2cmpt.rRMSE)==F,][1,]
    recommended.multi2<-sim.3cmpt.results.all[sim.3cmpt.results.all$sim.3cmpt.rRMSE== min(sim.3cmpt.results.all$sim.3cmpt.rRMSE,na.rm = T) & is.na(sim.3cmpt.results.all$sim.3cmpt.rRMSE)==F,][1,]
  }


  recommended_vmax_init<- recommended_mm$vmax
  recommended_km_init <-recommended_mm$km

  recommended_vc2cmpt_init <-recommended.multi1$vc
  recommended_vp2cmpt_init <-recommended.multi1$vp
  recommended_q2cmpt_init <-recommended.multi1$q

  recommended_vc3cmpt_init <-recommended.multi2$vc
  recommended_vp3cmpt_init <-recommended.multi2$vp
  recommended_vp23cmpt_init <-recommended.multi2$vp2
  recommended_q3cmpt_init <-recommended.multi2$q
  recommended_q23cmpt_init <-recommended.multi2$q2


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
  }

#################### Naive pooled data approach compartmental analysis #######
  if (run.option>0){

   dat<-rawdat  # recover raw data for compartmental analysis

    message(black(
      paste0("Run naive pooled data compartmental analysis ",strrep(".", 20))))

    if (npdcmpt.inits.strategy==0){
      message(red(
        paste0("Warning: there is no reference for intial estimates in current naive pooled data compartmental analysis, running maybe crashed down due to failure convergence")))

      # initial estimates = 1 (add a small random value to make the nls run)
      input.ka = round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.cl = round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vd = round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vc2cmpt=  round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vp2cmpt= round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vc3cmpt =  round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vp3cmpt = round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.vp23cmpt = round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.q2cmpt=    round(1 + runif(1, min = 0.01, max = 0.02),2)
      input.q3cmpt =   round(1+ runif(1, min = 0.01, max = 0.02),2)
      input.q23cmpt = round(1+ runif(1, min = 0.01, max = 0.02),2)
      input.vmax= round(1+ runif(1, min = 0.01, max = 0.02),2)
      input.km= round(1+ runif(1, min = 0.01, max = 0.02),2)
      vmax_km_threshold=F
    }

    if (npdcmpt.inits.strategy==1){
    input.ka = mean(base.ka.best)
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

     if (bolus_flag==1 || infusion_flag==1){

    dat[dat$CMT==1,]$CMT<-"centre"
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

     if (oral_flag==1){

    # change the input compartmental name consistent with model code.
    dat[dat$CMT==1,]$CMT="depot"
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
      est.method = est.method,
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

############ Summarise NPD compartmental analysis result#############################
    npd.1cmpt_results <- npd_1cmpt_out$npd.1cmpt_results
    npd.1cmpt.APE <- npd_1cmpt_out$npd.1cmpt.APE
    npd.1cmpt.MAE <- npd_1cmpt_out$npd.1cmpt.MAE
    npd.1cmpt.MAPE <- npd_1cmpt_out$npd.1cmpt.MAPE
    npd.1cmpt.RMSE <- npd_1cmpt_out$npd.1cmpt.RMSE
    npd.1cmpt.rRMSE <- npd_1cmpt_out$npd.1cmpt.rRMSE


    npd.1cmpt.mm_results <- npd_1cmpt_mm_out$npd.1cmpt.mm_results
    npd.1cmpt.mm.APE <- npd_1cmpt_mm_out$npd.1cmpt.mm.APE
    npd.1cmpt.mm.MAE <- npd_1cmpt_mm_out$npd.1cmpt.mm.MAE
    npd.1cmpt.mm.MAPE <- npd_1cmpt_mm_out$npd.1cmpt.mm.MAPE
    npd.1cmpt.mm.RMSE <- npd_1cmpt_mm_out$npd.1cmpt.mm.RMSE
    npd.1cmpt.mm.rRMSE <- npd_1cmpt_mm_out$npd.1cmpt.mm.rRMSE


    npd.2cmpt_results <- npd_2cmpt_out$npd.2cmpt_results
    npd.2cmpt.APE <- npd_2cmpt_out$npd.2cmpt.APE
    npd.2cmpt.MAE <- npd_2cmpt_out$npd.2cmpt.MAE
    npd.2cmpt.MAPE <- npd_2cmpt_out$npd.2cmpt.MAPE
    npd.2cmpt.RMSE <- npd_2cmpt_out$npd.2cmpt.RMSE
    npd.2cmpt.rRMSE <- npd_2cmpt_out$npd.2cmpt.rRMSE

    npd.3cmpt_results <- npd_3cmpt_out$npd.3cmpt_results
    npd.3cmpt.APE <- npd_3cmpt_out$npd.3cmpt.APE
    npd.3cmpt.MAE <- npd_3cmpt_out$npd.3cmpt.MAE
    npd.3cmpt.MAPE <- npd_3cmpt_out$npd.3cmpt.MAPE
    npd.3cmpt.RMSE <- npd_3cmpt_out$npd.3cmpt.RMSE
    npd.3cmpt.rRMSE <- npd_3cmpt_out$npd.3cmpt.rRMSE

   # Output
   # Naive pooled data approach (compartmental analysis)"
    if (bolus_flag==1 || infusion_flag==1){
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
      simMAE = c(npd.1cmpt.MAE, npd.1cmpt.mm.MAE, npd.2cmpt.MAE, npd.3cmpt.MAE),
      simMAPE = c(npd.1cmpt.MAPE, npd.1cmpt.mm.MAPE, npd.2cmpt.MAPE, npd.3cmpt.MAPE),
      simRMSE = c(npd.1cmpt.RMSE, npd.1cmpt.mm.RMSE, npd.2cmpt.RMSE, npd.3cmpt.RMSE),
      simrRMSE = c(npd.1cmpt.rRMSE, npd.1cmpt.mm.rRMSE, npd.2cmpt.rRMSE, npd.3cmpt.rRMSE),

      time.spent = c(
        npd.1cmpt_results$timespent,
        npd.1cmpt.mm_results$timespent,
        npd.2cmpt_results$timespent,
        npd.3cmpt_results$timespent
      )
    )
    }
    if (oral_flag==1){
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
        simMAE = c(npd.1cmpt.MAE, npd.1cmpt.mm.MAE, npd.2cmpt.MAE, npd.3cmpt.MAE),
        simMAPE = c(npd.1cmpt.MAPE, npd.1cmpt.mm.MAPE, npd.2cmpt.MAPE, npd.3cmpt.MAPE),
        simRMSE = c(npd.1cmpt.RMSE, npd.1cmpt.mm.RMSE, npd.2cmpt.RMSE, npd.3cmpt.RMSE),
        simrRMSE = c(npd.1cmpt.rRMSE, npd.1cmpt.mm.rRMSE, npd.2cmpt.rRMSE, npd.3cmpt.rRMSE),

        time.spent = c(
          npd.1cmpt_results$timespent,
          npd.1cmpt.mm_results$timespent,
          npd.2cmpt_results$timespent,
          npd.3cmpt_results$timespent
        )
      )
   }
}

  ######################## Finally selection########################################

  if (run.option<2){
   # Part 1. ka,cl,vd.
    f_init_ka <- base.best$`Calculated Ka`[1]
    f_init_cl <- base.best$`Calculated CL`[1]
    f_init_vd <- base.best$`Calculated Vd`[1]

    sel.method.ka.cl.vd <-base.best$Method[1]

    if (nrow(base.best)>1){
      #check the total volume of distribution of two compartment
      total.vd <- recommended.multi1$vc+recommended.multi1$vp

      f_init_cl <- base.best[base.best$`Calculated Vd`== round(total.vd,1), ]$`Calculated CL`
      f_init_vd <- base.best[base.best$`Calculated Vd`== round(total.vd,1), ]$`Calculated Ka`
      sel.method.ka.cl.vd <-base.best[base.best$`Calculated Vd`== round(total.vd,1), ]$Method
    }


    # vmax.km
    f_init_vmax <- recommended_vmax_init
    f_init_km <- recommended_km_init
    sel.method.vmax.km <- "Simulation-based analysis "

    # Multi-compartmental parameters
    f_init_vc2cmpt <-  recommended_vc2cmpt_init
    f_init_vp2cmpt <-  recommended_vp2cmpt_init
    f_init_q2cmpt <-   recommended_q2cmpt_init
    f_init_vc3cmpt <- recommended_vc3cmpt_init
    f_init_vp3cmpt <- recommended_vp3cmpt_init
    f_init_vp23cmpt <- recommended_vp23cmpt_init
    f_init_q3cmpt <- recommended_q3cmpt_init
    f_init_q23cmpt <-recommended_q23cmpt_init

    sel.method.multi <- "Simulation-based analysis "
    sel.method.ka<-"Wanger-nelson method"

    if (sel.method.ka.cl.vd== "Single-point method" & oral_flag==1 ){
      sel.method.ka<-"Single-point method"
    }

    if (sel.method.ka.cl.vd== "Graphic analysis"){
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



  # if (!exists("init.messages.multi")) {
  #   init.messages.multi <- NULL
  # }
  # if (!exists("init.messages.vmax.km")) {
  #   init.messages.vmax.km <- NULL
  # }
  #
  # if (is.null(init.messages.vmax.km) &
  #     is.null(init.messages.vmax.km)) {
  #   init.messages <- NULL
  # }
  #
  # if (is.null(init.messages.vmax.km) == F ||
  #     is.null(init.messages.vmax.km) == F) {
  #   init.messages <- c(init.messages.vmax.km, init.messages.multi)
  # }

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
      "Absolute Predicted Error (APE)",
      "Mean Absolute Error (MAE)",
      "Mean Absolute Percentage Error (MAPE)",
      "Root Mean Squared Error (RMSE)",
      "Relative Root Mean Squared Error (rRMSE)",
      "Time spent"
    )
  colnames(sim.vmax.km.results.all) <-
    c(
      "Simulated Vmax",
      "Simulated Km",
      "Absolute Predicted Error (APE)",
      "Mean Absolute Error (MAE)",
      "Mean Absolute Percentage Error (MAPE)",
      "Root Mean Squared Error (RMSE)",
      "Relative Root Mean Squared Error (rRMSE)",
      "Time spent"
    )
  colnames(sim.2cmpt.results.all) <-
    c(
      "Simulated Vc",
      "Simulated Vp",
      "Simulated Q",
      "Absolute Predicted Error (APE)",
      "Mean Absolute Error (MAE)",
      "Mean Absolute Percentage Error (MAPE)",
      "Root Mean Squared Error (RMSE)",
      "Relative Root Mean Squared Error (rRMSE)",
      "Time spent"
    )
  colnames(sim.3cmpt.results.all) <-
    c(
      "Simulated Vc",
      "Simulated Vp",
      "Simulated Vp2",
      "Simulated Q",
      "Simulated Q2",
      "Absolute Predicted Error (APE)",
      "Mean Absolute Error (MAE)",
      "Mean Absolute Percentage Error (MAPE)",
      "Root Mean Squared Error (RMSE)",
      "Relative Root Mean Squared Error (rRMSE)",
      "Time spent"
    )

  init.history <- list(
    base.out = all.out,
    single.point.lst=single.point.lst,
    nca.results=nca.results,
    ka_wanger_nelson_result=    ka_wanger_nelson_result,
    graph.results_fd=   graph.results_fd,
    sim.vmax.km = sim.vmax.km.results.all,
    sim.2cmpt = sim.2cmpt.results.all,
    sim.3cmpt = sim.3cmpt.results.all,
    npd_1cmpt_out =npd_1cmpt_out,
    npd_1cmpt_mm_out = npd_1cmpt_mm_out ,
    npd_2cmpt_out=npd_2cmpt_out,
    npd_3cmpt_out= npd_3cmpt_out
  )

  params.descriptions <- c(
    "Ka: absorption constant rate",
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
  }

  if (run.option==2){

    init.history <- list(
      npd_1cmpt_out =npd_1cmpt_out,
      npd_1cmpt_mm_out = npd_1cmpt_mm_out ,
      npd_2cmpt_out=npd_2cmpt_out,
      npd_3cmpt_out= npd_3cmpt_out
    )

    params.descriptions <- c(
      "Ka: absorption constant rate",
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
    output_env$npdcmpt.all.out<-npdcmpt.all.out
    output_env$Run.history <- init.history
    output_env$Parameter.descriptions <- params.descriptions
  }

  class(output_env) <- "getPPKinits"

  gc()

  return(output_env)

} # end of function


#' Print environment summary for initial parameter estimation
#'
#' Prints a summary of initial parameter estimation including Data information, recommended initial estimates by pipeline, naive pooled data compartmental analysis, run history and parameter descriptions stored in an environment by S3 method.
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
# Define a custom print method for the 'getPPKinits' by S3 method
print.getPPKinits <- function(env, ...) {
  cat("===============Initial Parameter Estimation Summary ===============\n")
  cat("Data information:\n")
  message(black(env$Datainfo))

  cat("\nRecommended initial estimates by using non-model fitting methods:\n")
  print(head(env$Recommended_initial_estimates, 13))

  cat("\nRecommended initial estimates by naive pooled data compartmental analysis:\n")
  print(head(env$npdcmpt.all.out, 13))

  cat("\nParameter descriptions:\n")
  print(env$Parameter.descriptions)

  cat("\n=============== End of Summary ===============\n")
}

# plot.getPPKinits
