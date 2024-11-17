#' Run non-compartmental analysis (NCA) for intravenous pharmacokinetic data
#'
#' Perform a non-compartmental analysis (NCA) for intravenous dosing data with normalised concentration with dose, calculating clearance, volume of distribution, slope, and half-life for all pooled data, as well as for only pooled first-dose data if available.
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data which include at least columns for ID, TIME, EVID, AMT, DV, and dose_number.
#' @param nlastpoints The number of last points to be used for slope calculation.
#' @param trapezoidal.rule A numeric value indicating the method used for calculating AUC using the trapezoidal rule.
#'        1 for the standard linear trapezoidal method.
#'        2 for the linear up and logarithmic down method, where the concentration increase is treated linearly,
#'        and the concentration decrease is treated logarithmically.
#' @param nbins The number of bins for time binning.
#' @param fdobsflag A flag indicating whether the dataset includes first-dose observations. Set to 1 if the data includes first-dose observations, otherwise set to 0.
#' @return A list containing the results of the NCA calculation for the first dose (`nca.fd.results`) and for all data (`nca.results`).
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @importFrom tidyr fill
#' @import nlmixr2
#' @examples
#' dat <- theo_sd
#' dat <- processData(dat)$dat
#' run_nca.normalised(dat, nlastpoints = 3, trapezoidal.rule=1,nbins = 8, fdobsflag = 1,sdflag=1)
#' run_nca.normalised(dat, nlastpoints = 3, trapezoidal.rule=2,nbins = 8, fdobsflag = 1,sdflag=1)
#' @export

run_nca.normalised <- function(dat,
                               nlastpoints=3,
                               trapezoidal.rule=1,
                               nbins=8,
                               fdobsflag=1,
                               sdflag=1,
                               route="bolus") {
  # default settings
  nca.fd.results <- data.frame(
    cl = NA,
    vd = NA,
    slope = NA,
    half_life = NA,
    auct=NA,
    auc0_inf=NA,
    C_last=NA,
    ke=NA,
    aumc_0_inf=NA,
    start.time = NA,
    time.spent = 0,
    nlastpoints=nlastpoints
  )

  nca.efd.results <- data.frame(
    cl = NA,
    vd = NA,
    slope = NA,
    half_life = NA,
    auct=NA,
    auc0_inf=NA,
    C_last=NA,
    ke=NA,
    aumc_0_inf=NA,
    start.time = NA,
    time.spent = 0,
    nlastpoints=nlastpoints
  )

  nca.all.results <- data.frame(
    cl = NA,
    vd = NA,
    slope = NA,
    half_life = NA,
    auct=NA,
    auc0_inf=NA,
    C_last=NA,
    ke=NA,
    aumc_0_inf=NA,
    start.time = NA,
    time.spent = 0,
    nlastpoints=nlastpoints
  )

  datpooled_fd<-NA
  datpooled_efd<-NA
  datpooled_all<-NA

  # Analyse data after the first dose
  if (fdobsflag==1){

    start.time <- Sys.time()

    dat_fd <- dat[dat$dose_number == 1 & dat$iiobs==0, ]
    datpooled_fd <- pk.time.binning(testdat = dat_fd,
                                    nbins = nbins)
    # clear nca.output
    nca.output<-NA
    nca.output <-
      nca.iv.normalised(dat = datpooled_fd$test.pool.normalised,
                        trapezoidal.rule= trapezoidal.rule,
                        nlastpoints = nlastpoints,
                        route=route)

    end.time <- Sys.time()
    time.spent <- round(difftime(end.time, start.time), 4)

    nca.fd.results <- data.frame(
      cl = signif(nca.output[1], 3),
      vd = signif(nca.output[2], 3),
      slope = signif(nca.output[3], 3),
      half_life = signif(nca.output[4], 3),
      auct = signif(nca.output[5], 3),
      auc0_inf = signif(nca.output[6], 3),
      C_last = signif(nca.output[7], 3),
      ke = signif(nca.output[8], 3),
      aumc_0_inf = signif(nca.output[9], 3),
      start.time = start.time,
      time.spent = time.spent,
      nlastpoints=nca.output[10]
    )
    }


  # Analyse data after the repeated dose

  if (sdflag == 0) {

    # Needs to introduce most commonly used dose-interval
    # Calculate the most commonly used dose interval
    dose_data <- dat[dat$EVID %in% c(1, 4) & dat$AMT > 0, ]
    dose_data <- dose_data[order(dose_data$ID, dose_data$TIME), ]
    dose_data$interval <-
      ave(
        dose_data$TIME,
        dose_data$ID,
        FUN = function(x)
          c(NA, diff(x))
      )
    dose_intervals <-
      round(dose_data$interval[!is.na(dose_data$interval)], 0)
     most_commonly_used_dose_interval <-
      as.numeric(names(sort(table(dose_intervals), decreasing = TRUE)[1]))

      start.time <- Sys.time()

      dat_efd1.obs <- dat[dat$dose_number != 1 & dat$EVID==0, ]
      dat_efd1.dose <- dat[dat$dose_number != 1 & dat$EVID==1, ]
      dat_efd1.obs <-dat_efd1.obs[ dat_efd1.obs$tad<=most_commonly_used_dose_interval*1.2, ]
      dat_efd2 <- dat[dat$dose_number==1 & dat$iiobs>0,]

      dat_efd<-rbind( dat_efd1.obs, dat_efd1.dose, dat_efd2)
      dat_efd<- dat_efd[with(dat_efd, order(ID, resetflag, TIME, -AMT)), ]

      datpooled_efd <- pk.time.binning(testdat = dat_efd,
                                       nbins = nbins)
      nca.output<-NA
      nca.output <-
        nca.iv.normalised(dat = datpooled_efd$test.pool.normalised,
                           trapezoidal.rule= trapezoidal.rule,
                          nlastpoints = nlastpoints,
                          ss = 1,route=route)

      end.time <- Sys.time()
      time.spent <- round(difftime(end.time, start.time), 4)

      nca.efd.results <- data.frame(
        cl = signif(nca.output[1], 3),
        vd = signif(nca.output[2], 3),
        slope = signif(nca.output[3], 3),
        half_life = signif(nca.output[4], 3),
        auct = signif(nca.output[5], 3),
        auc0_inf = signif(nca.output[6], 3),
        C_last = signif(nca.output[7], 3),
        ke = signif(nca.output[8], 3),
        aumc_0_inf = signif(nca.output[9], 3),
        start.time = start.time,
        time.spent = time.spent,
        nlastpoints=nca.output[10]
      )
    }

  # Analyse data after the first dose and repeated dose

  if (fdobsflag==1 & sdflag==0){

  start.time <- Sys.time()

  dat_all.obs <- dat[dat$EVID==0, ]
  dat_all.dose <- dat[dat$EVID==1, ]
  dat_all.obs <-dat_all.obs[ dat_all.obs$tad<=most_commonly_used_dose_interval*1.2, ]

  dat_all<-rbind( dat_all.obs, dat_all.dose)
  dat_all<- dat_all[with(dat_efd, order(ID, resetflag, TIME, -AMT)), ]

  datpooled_all <- pk.time.binning(testdat = dat_all,
                                   nbins = nbins)

  nca.output<-NA
  nca.output <-
    nca.iv.normalised(dat = datpooled_all$test.pool.normalised,
                       trapezoidal.rule= trapezoidal.rule,
                      nlastpoints = nlastpoints,
                      ss=1,route=route)

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  nca.all.results  <- data.frame(
    cl = signif(nca.output[1], 3),
    vd = signif(nca.output[2], 3),
    slope = signif(nca.output[3], 3),
    half_life = signif(nca.output[4], 3),
    auct = signif(nca.output[5], 3),
    auc0_inf = signif(nca.output[6], 3),
    C_last = signif(nca.output[7], 3),
    ke = signif(nca.output[8], 3),
    aumc_0_inf = signif(nca.output[9], 3),
    start.time = start.time,
    time.spent = time.spent,
    nlastpoints=nca.output[10]
  )
  }

  return(list(nca.fd.results = nca.fd.results,
              nca.efd.results = nca.efd.results,
              nca.all.results = nca.all.results,
              datpooled_fd=datpooled_fd,
              datpooled_efd=datpooled_efd,
              datpooled_all=datpooled_all
              ))

}


#' Non-compartmental analysis for intravenous pharmacokinetic data
#'
#' Perform non-compartmental analysis (NCA) on intravenous data to calculate pharmacokinetic parameters (clearance and volume of distribution) from the provided data.
#' @param dat A data frame with at least two columns: TIME and DV.
#' @param nlastpoints Number of last points to use for the linear regression on terminal slope (default is 4).
#' @param trapezoidal.rule A numeric value indicating the method used for calculating AUC using the trapezoidal rule.
#'        1 for the standard linear trapezoidal method.
#'        2 for the linear up and logarithmic down method, where the concentration increase is treated linearly,
#'        and the concentration decrease is treated logarithmically.
#' @param ss A flag indicating whether the analysis assumes steady state (TRUE or FALSE).
#' \itemize{
#'   \item \strong{FALSE (0)}: The model uses the extrapolated AUC (auc0_inf) to calculate clearance (CL) and volume of distribution (Vd).
#'   \item \strong{TRUE (1)}: If at steady state, the model uses only the AUC up to the last time point (auct) for calculating clearance (CL) and volume of distribution (Vd), instead of auc0_inf.
#' }
#' @return A named vector containing the calculated clearance (cl), volume of distribution (vd), slope, and half-life.
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(12, 8, 5, 3, 2, 1.5, 1 ))
#' nca.iv.normalised(dat, nlastpoints = 3)
#'
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(12, 8, 5, 3, 2, 1.5, 1 ))
#' nca.iv.normalised(dat, ss=1, nlastpoints = 3)
#'
#' @export

nca.iv.normalised <- function(dat,
                               trapezoidal.rule=1,
                              ss=0,
                              nlastpoints=3,
                              route=route) {

  cl=NA
  vd=NA
  slope=NA
  half_life=NA
  auct=NA
  auc0_inf=NA
  C_last=NA
  ke=NA
  aumc_0_inf=NA
  clnormalised=NA

  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  if ( trapezoidal.rule==1){
    auct <- trap.rule(dat$TIME, dat$DV)
  }
  if ( trapezoidal.rule==2){
    auct <- trap.rule_linear_up_log_down(dat$TIME, dat$DV)
  }

  # Identify the index of Tmax
  max_index <- which.max(dat$DV)

  if (route == "bolus" || route == "infusion") {
    temp1 <-
      dat[max_index:nrow(dat),] # Subset data points after Tmax
  }

  if (route == "oral") {
    temp1 <-
      dat[(max_index + 1):nrow(dat),] # Subset data points after Tmax
  }

  # Loop to adjust nlastpoints if there are insufficient points
  while (nlastpoints > 1 && nrow(temp1) < nlastpoints) {
    nlastpoints <- nlastpoints - 1
  }

  # check number of points available for slope linear regression
  if (nlastpoints<2){
    return(c(cl=cl,
             vd=vd,
             slope=slope,
             half_life=half_life,
             auct=auct,
             auc0_inf=auc0_inf,
             C_last=C_last,
             ke=ke,
             aumc_0_inf= aumc_0_inf,
             nlastpoints=nlastpoints))
  }

  # Select last 4/specified number for slope calculation
  temp1 <- tail(dat, n = nlastpoints)
  C_last <- tail(temp1$DV, 1)
  t_last <- tail(temp1$TIME, 1)

  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)
  slope <- summary(abc)[[4]][[2]]

  if (ss==1){
   clnormalised <- 1 / auct
   cl <- clnormalised
  }

  if (slope<0){
  ke <- -slope
  lambda_z<- ke
  half_life <- log(2) / ke
  auct_inf <- C_last / ke
  auc0_inf <- auct + auct_inf

  if (ss==0){
   clnormalised <- 1 / auc0_inf
   cl<-clnormalised
  }
  vdnormalised <- clnormalised / ke

  vd <- vdnormalised
  # AUMC calculation
  time<-dat$TIME
  concentration<-dat$DV
  moment_curve <- time * concentration
  # Calculate AUMC from 0 to t_last
  aumc0_t <- trap.rule(time, moment_curve)
  # Calculate AUMC from t_last to infinity
  aumc_tlast_to_inf <- (C_last * t_last) / lambda_z + C_last / (lambda_z^2)
  # Calculate AUMC from 0 to infinity
  aumc_0_inf<- aumc0_t + aumc_tlast_to_inf
  }

  return(c(cl=cl,
           vd=vd,
           slope=slope,
           half_life=half_life,
           auct=auct,
           auc0_inf=auc0_inf,
           C_last=C_last,
           ke=ke,
           aumc_0_inf= aumc_0_inf,
           nlastpoints=nlastpoints))
}




#' Calculate Area Under the Curve (AUC) using the trapezoidal rule
#'
#' This function computes the area under the curve (AUC) using the linear trapezoidal rule.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration at each time point.
#'
#' @return A numeric value representing the estimated AUC using the trapezoidal rule.
#'
#' @examples
#' x <- c(0.5, 1, 2, 4, 6, 8, 10)
#' y <- c(12, 8, 5, 3, 2, 1.5, 1)
#' trap.rule(x, y)
#'

trap.rule <- function(x, y) {
  # Calculate the first triangle (between the first two points)
  first_triangle_area <- x[1] * y[1] / 2

  # Calculate the trapezoidal areas for the middle points (excluding first and last segments)
  trapezoidal_area <- sum(diff(x) * (y[-1] + y[-length(y)])) / 2

  # Total area is the sum of the first triangle and middle trapezoids.
  total_area <- first_triangle_area + trapezoidal_area

  return(total_area)
}

#' Calculate AUC using linear up and log down trapezoidal rule
#'
#' This function computes the area under the curve (AUC) by using the trapezoidal rule
#' for phases where concentration is increasing, and a logarithmic rule for phases
#' where concentration is decreasing.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration at each time point.
#'
#' @return A numeric value representing the estimated AUC using the linear up/log down method.
#'
#' @examples
#' x <- c(0.5, 1, 2, 4, 6, 8, 10)
#' y <- c(12, 8, 5, 3, 2, 1.5, 1)
#' trap.rule_linear_up_log_down(x, y)
#'
#' @export
#
trap.rule_linear_up_log_down <- function(x, y) {

  # Calculate the first triangle (between the first two points)
  first_triangle_area <- x[1] * y[1] / 2

  delta_x <- diff(x)
  delta_y <- diff(y)

  linear_up_ <- delta_y > 0
  log_down_ <- delta_y < 0

  # For the linear up phase
  linear_up_auc <- (y[-length(y)][linear_up_] + y[-1][linear_up_]) / 2 * delta_x[linear_up_]

  # For the log down phase
  log_down_auc <- ((y[-length(y)][log_down_] - y[-1][log_down_]) /
                     log(y[-length(y)][log_down_] / y[-1][log_down_])) * delta_x[log_down_]

  # Sum the AUC from both the linear up and log down phases
  total_auc <- first_triangle_area + sum(linear_up_auc) + sum(log_down_auc)

  return(total_auc)
}















