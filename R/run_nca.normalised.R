#' Run non-compartmental analysis (NCA) for intravenous pharmacokinetic data
#'
#' Perform a non-compartmental analysis (NCA) for intravenous dosing data with normalised concentration with dose, calculating clearance, volume of distribution, slope, and half-life for all pooled data, as well as for only pooled first-dose data if available.
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data which include at least columns for ID, TIME, EVID, AMT, DV, and dose_number.
#' @param nlastpoints The number of last points to be used for slope calculation.
#' @param trap.rule.method A numeric value indicating the method used for calculating AUC using the trapezoidal rule.
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
#' dat <- Bolus_1CPT
#' dat <- nmpkconvert(dat)
#' dat <- calculate_tad(dat)
#' run_nca.normalised(dat, nlastpoints = 3, trap.rule.method=1,nbins = 8, fdobsflag = 1,sdflag=0)
#' run_nca.normalised(dat, nlastpoints = 3, trap.rule.method=2,nbins = 8, fdobsflag = 1,sdflag=0)
#' @export

run_nca.normalised <- function(dat,
                                  nlastpoints,
                                  trap.rule.method,
                                  nbins,
                                  fdobsflag,
                                  sdflag) {
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
    time.spent = 0
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
    time.spent = 0
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
    time.spent = 0
  )

  datpooled_fd<-NA
  datpooled_efd<-NA
  datpooled_all<-NA

  # Analyse data after the first dose
  if (fdobsflag == 1) {
    start.time <- Sys.time()
    dat$DVnor <- dat$DV / dat$dose

    dat_fd <- dat[dat$dose_number == 1, ]
    datpooled_fd <- pk.time.binning(testdat = dat_fd,
                                    nbins = nbins)

    nca.output <-
      nca.iv.normalised(dat = datpooled_fd$test.pool.normalised,
                        trap.rule.method=trap.rule.method,
                        nlastpoints = nlastpoints)

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
      time.spent = time.spent
    )
  }

  # Analyse data after the repeated dose

  if (sdflag == 0) {
      start.time <- Sys.time()
      dat$DVnor <- dat$DV / dat$dose

      dat_efd <- dat[dat$dose_number != 1, ]
      datpooled_efd <- pk.time.binning(testdat = dat_efd,
                                       nbins = nbins)

      nca.output <-
        nca.iv.normalised(dat = datpooled_efd$test.pool.normalised,
                          trap.rule.method=trap.rule.method,
                          nlastpoints = nlastpoints)

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
        time.spent = time.spent
      )
    }

  # Analyse data after the first dose and repeated dose

  if (fdobsflag==1 & sdflag==0){

  start.time <- Sys.time()
  dat$DVnor <- dat$DV / dat$dose
  datpooled_all <- pk.time.binning(testdat = dat,
                                   nbins = nbins)

  datpooled_all$test.pool.normalised
  nca.output <-
    nca.iv.normalised(dat = datpooled_all$test.pool.normalised,
                      trap.rule.method=trap.rule.method,
                      nlastpoints = nlastpoints)
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
    time.spent = time.spent
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
#' @return A named vector containing the calculated clearance (cl), volume of distribution (vd), slope, and half-life.
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(12, 8, 5, 3, 2, 1.5, 1 ))
#' calc.nca.iv.normalised(dat, nlastpoints = 4)
#' @export

nca.iv.normalised <- function(dat,
                              trap.rule.method,
                              nlastpoints) {

  if (missing(trap.rule.method)){
    trap.rule.method=1
  }
  if (missing(nlastpoints)) {
    nlastpoints <- 4
  }
  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  if (trap.rule.method==1){
    auct <- trap.rule(dat$TIME, dat$DV)
  }
  if (trap.rule.method==2){
    auct <- trap.rule_linear_up_log_down(dat$TIME, dat$DV)
  }

  # Select last 4/specified number for slope calculation
  temp1 <- tail(dat, n = nlastpoints)
  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)
  slope <- summary(abc)[[4]][[2]]
  ke <- -slope
  lambda_z<- ke
  half_life <- log(2) / ke
  C_last <- tail(temp1$DV, 1)
  t_last <- tail(temp1$TIME, 1)

  auct_inf <- C_last / ke
  auc0_inf <- auct + auct_inf
  clnormalised <- 1 / auc0_inf
  vdnormalised <- clnormalised / ke
  cl <- clnormalised
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

  return(c(cl, vd, slope, half_life, auct, auc0_inf, C_last, ke, aumc_0_inf))
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
  total_auc <- sum(linear_up_auc) + sum(log_down_auc)

  return(total_auc)
}















