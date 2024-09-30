#' Calculate the Absorption Rate Constant (ka) using the Wanger-Nelson Method
#'
#' Computes the absorption rate constant (ka) using the Wanger-Nelson method.
#' The function makes use of non-compartmental analysis (NCA) function part to calculate the area under the curve (AUC) and the elimination rate constant (ke).
#' Linear regression on the logarithm of the remaining fraction of drug is used to estimate the absorption rate constant (ka).
#'
#' @param dat A data frame with at least two columns: `DV` (drug concentration)
#' and `TIME` (time points of measurement).
#' @param nlastpoints An integer specifying the number of last points to be used for slope calculation
#' of the elimination rate constant (ke). However, this argument is not used in the current function logic but remains for future enhancements.
#'
#' @details
#' The function calculates the absorption rate constant (ka) as follows:
#' \itemize{
#'   \item The function uses the `nca.iv.normalised` method to calculate the area under the curve (AUC) and elimination rate constant (ke).
#'   \item AUC is calculated using trapezoidal integration for the observed data points.
#'   \item The fraction of drug absorbed is calculated based on ke and the cumulative AUC.
#'   \item The logarithm of the remaining fraction of the drug is then modeled using linear regression to determine the absorption rate constant (ka).
#' }
#'
#' This method assumes that the drug absorption follows first-order kinetics.
#'
#' @return A numeric value representing the estimated absorption rate constant (ka).
#'
#' @examples
#' # Example 1.  data frame with drug concentration (DV) and time (TIME)
#' data_example <- data.frame(
#'   TIME = c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10),
#'   DV = c(0, 5, 10, 30, 40, 45, 30, 20, 10, 5 )
#' )
#'
#' # Calculate ka using the Wanger-Nelson method with the last 4 data points
#' ka_wanger_nelson(data_example, 4)$ka
#'
#' # Example 2. data frame from nlmixr2data
#' dat<-Oral_1CPT[Oral_1CPT$ID==1 &Oral_1CPT$SD==1&Oral_1CPT$EVID==0,]
#' dat<-data.frame(TIME=dat$TIME,DV=dat$DV)
#' ka_wanger_nelson(dat=dat,nlastpoints=4)$ka
#' ka_wanger_nelson(dat=dat,nlastpoints=4)$dat_out_wanger_nelson
#'
#' @export
#'
ka_wanger_nelson<-function(dat,nlastpoints,nca.out){

  colnames(dat)[1]<-"TIME"
  colnames(dat)[2]<-"DV"
  x<-dat$TIME
  y<-dat$DV
  if (missing(nca.out)){
  nca.out<-nca.iv.normalised(dat = data.frame(time=x,conc=y),nlastpoints=nlastpoints)
  }
  auc_ <- nca.out[6]
  ke<- nca.out[8]
  auc_intervals<- diff(x) * (y[-1] + y[-length(y)])/2

  # use a triangle to calculate 0 - first sample
  dat$auc_intervals<-c( dat[1,]$DV * dat[1,]$TIME/2, auc_intervals)
  dat$auc_accumulate<-NA

  for (k in 1:nrow(dat)){
    if (k==1){
      dat[k,]$auc_accumulate<-dat[1,]$auc_intervals
    }
    else{
      dat[k,]$auc_accumulate<-sum(dat[dat$TIME<=dat[k,]$TIME, ]$auc_intervals)
    }
  }
  dat$frac_abs<-  (dat$DV+ ke*dat$auc_accumulate)/ (ke*auc_)

  # max_abs_time <- which(dat$frac_abs == max(dat$frac_abs ), arr.ind = T)
  # dat_absorb_phase <- dat[1:max_abs_time, ]
  # dat_absorb_phase$frac_remained<-1-  dat_absorb_phase$frac_abs
  dat$frac_remained<-1-  dat$frac_abs
  dat_absorb_phase<- dat[dat$frac_remained>0.1,]
  # dat_absorb_phase<- dat_absorb_phase[dat_absorb_phase$frac_remained>0.05,]
  abc2 <- lm(log(dat_absorb_phase$frac_remained) ~  dat_absorb_phase$TIME)
  slope_ka<-summary(abc2)[[4]][[2]]
  ka= -slope_ka
  return(list(ka=ka,dat_out_wanger_nelson=dat))
}


#' Calculate absorption rate constant (ka) using statistical moments
#'
#' Calculates the absorption rate constant (ka) for oral drug administration using statistical moments.
#' It calculates the area under the concentration-time curve (AUC) and the area under the moment curve (AUMC) based on non-compartmental analysis (NCA).
#' The absorption rate constant is derived from the mean absorption time (MAT), which is calculated using the difference between
#' the mean residence time (MRT) for oral and intravenous (IV) administration.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the concentrations corresponding to the time points.
#' @param nlastpoints An integer specifying the number of last points to be used for the slope calculation of the terminal elimination rate constant.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item It uses the `nca.iv.normalised` function to calculate the area under the curve (AUC), the area under the moment curve (AUMC), and the mean residence time (MRT) for IV administration.
#'   \item The MRT for oral administration is calculated as the ratio of AUMC to AUC.
#'   \item The mean absorption time (MAT) is calculated as the difference between the MRT for oral and IV administration.
#'   \item Finally, the absorption rate constant (ka) is estimated as the inverse of MAT.
#' }
#'
#' @return A numeric value representing the estimated absorption rate constant (ka).
#'
#' @examples
#' # Example data from Oral_1CPT dataset (first subject and single dose)
#' dat <- Oral_1CPT[Oral_1CPT$ID == 1 & Oral_1CPT$SD == 1 & Oral_1CPT$EVID == 0, ]
#' dat <- data.frame(TIME = dat$TIME, DV = dat$DV)
#'
#' # Calculate ka using statistical moments with the last 4 points for lambda_z
#' calculate_ka_statistical_moments(x = dat$TIME, y = dat$DV, nlastpoints = 4)
#'
#' @export
#'
calculate_ka_statistical_moments <- function(x, y, nlastpoints) {
  # Calculate AUC and AUMC for i.v. administration
  nca.out<-nca.iv.normalised(dat = data.frame(time=x,conc=y),nlastpoints=nlastpoints)
  auc_ <- nca.out[6]
  aumc_ <- nca.out[9]
  mrt_oral <- aumc_ / auc_
  # Calculate MRT for iv, Vss = CL * MRT
  mrt_iv <- nca.out[4]/log(2)
  # Calculate Mean Absorption Time (MAT)
  mat <- mrt_oral - mrt_iv
  # Estimate absorption rate constant (ka)
  ka <- 1 / mat
  return(ka)
}

