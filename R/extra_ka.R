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
  if (nrow(dat_absorb_phase)<2){
    ka=NA
  }

  if (nrow(dat_absorb_phase)>1){
  abc2 <- lm(log(dat_absorb_phase$frac_remained) ~  dat_absorb_phase$TIME)
  slope_ka<-summary(abc2)[[4]][[2]]
  ka= -slope_ka
  }
  return(list(ka=ka,dat_out_wanger_nelson=dat))
}



#' Calculate Absorption Rate Constant (ka) in a One-Compartment Model
#'
#' Calculates the absorption rate constant (\code{ka}) for a drug administered orally
#' using a one-compartment model with first-order absorption and elimination kinetics.
#' The calculation is based on known values of clearance (\code{cl}), elimination rate constant
#' (\code{ke}), time (\code{t}), observed concentration (\code{Ct}), bioavailability (\code{Fbio}),
#' and dose (\code{Dose}).
#'
#' @param cl Numeric. Clearance of the drug (in L/hr).
#' @param ke Numeric. Elimination rate constant (in 1/hr).
#' @param t Numeric. Time after drug administration (in hours) at which the concentration is measured.
#' @param Ct Numeric. Observed concentration of the drug at time \code{t} (in mg/L).
#' @param Fbio Numeric. Bioavailability fraction (default = 1, meaning 100% bioavailability).
#' @param Dose Numeric. Administered dose of the drug (in mg).
#'
#' @return A list containing the following components:
#' \item{ka}{The calculated absorption rate constant.}
#' \item{full_solution}{The full solution object returned by the \code{uniroot()} function, which includes additional details about the root-finding process.}
#'
#' @details
#' This function uses a one-compartment model with first-order absorption and elimination to
#' estimate the absorption rate constant (\code{ka}). The function solves the following equation
#' for \code{ka} numerically using \code{uniroot()}:
#'
#' \deqn{Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left( e^{-ke \cdot t} - e^{-ka \cdot t} \right)}
#'
#' The \code{uniroot()} function is used to find the value of \code{ka} that makes the difference
#' between the predicted and observed concentrations equal to zero. The reasonable range for
#' \code{ka} is set as [0.01, 100] for root-finding.
#'
#' @examples
#' # Example usage:
#' # Data from first point in simulated dataset Oral_1CPT
#' ka_result <- ka_calculation(cl = 4, ke = 0.057, t = 0.25, Ct = 204.8, Dose = 60000)
#' print(ka_result$ka)
#'
#' @export
#
ka_calculation <- function(cl,     # Clearance of the drug (L/hr)
                           ke,     # Elimination rate constant (1/hr)
                           t,      # Time (hr) after drug administration at which concentration is measured
                           Ct,     # Observed concentration of the drug at time 't' (mg/L)
                           Fbio=1, # Bioavailability fraction, default is 1 (100% bioavailability)
                           Dose) { # Administered dose of the drug (mg)

  # Define the equation to solve for the absorption rate constant (ka)
  ka.equation <- function(ka) {

    # Predicted concentration using the one-compartment model with first-order absorption and elimination
    predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
      (exp(-ke * t) - exp(-ka * t))

    # Return the difference between predicted concentration and observed concentration
    return(predicted_Ct - Ct)
  }

  # Use the uniroot function to numerically solve for ka (absorption rate constant)
  # The reasonable range for ka is provided as [lower, upper] for root finding
  solution <- uniroot(ka.equation, lower = 0.01, upper = 100)

  # Return the value of ka (solution$root) and full solution object
  return(list(ka = solution$root, full_solution = solution))
}






