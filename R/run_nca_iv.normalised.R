#' Run non-compartmental analysis (NCA) for intravenous pharmacokinetic data
#'
#' Perform a non-compartmental analysis (NCA) for intravenous dosing data with normalised concentration with dose, calculating clearance, volume of distribution, slope, and half-life for all pooled data, as well as for only pooled first-dose data if available.
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data which include at least columns for ID, TIME, EVID, AMT, DV, and dose_number.
#' @param nlastpoints The number of last points to be used for slope calculation.
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
#' dat$DVnor<-dat$DV/dat$dose
#' run_nca_iv.normalised(dat, nlastpoints = 3, nbins = 8, fdobsflag = 1)
#' @export

run_nca_iv.normalised <- function(dat,
                                  nlastpoints,
                                  nbins,
                                  fdobsflag) {
  # default settings
  nca.fd.results <- data.frame(
    cl = NA,
    vd = NA,
    slope = NA,
    half_life = NA,
    start.time = NA,
    time.spent = 0
  )
  
  # If there are multiple doses, consider pooling the data of the first dose together first
  if (fdobsflag == 1) {
    start.time <- Sys.time()
    dat$DVnor <- dat$DV / dat$dose
    
    dat_fd <- dat[dat$dose_number == 1, ]
    datpooled_fd <- pk.time.binning(testdat = dat_fd,
                                    nbins = nbins)
    
    nca.output <-
      nca.iv.normalised(dat = datpooled_fd$test.pool.normalised,
                        nlastpoints = nlastpoints)
    
    end.time <- Sys.time()
    time.spent <- round(difftime(end.time, start.time), 4)
    
    nca.fd.results <- data.frame(
      cl = signif(nca.output[1], 3),
      vd = signif(nca.output[2], 3),
      slope = signif(nca.output[3], 3),
      half_life = signif(nca.output[4], 3),
      start.time = start.time,
      time.spent = time.spent
    )
    
    
    # Ncacalplot(dat = dat,datpooled = datpooled_fd,ncasubfd = T,nlastpoints=nlastpoints)
  }
  
  start.time <- Sys.time()
  dat$DVnor <- dat$DV / dat$dose
  datpooled_all <- pk.time.binning(testdat = dat,
                                   nbins = nbins)
  
  datpooled_all$test.pool.normalised
  nca.output <-
    nca.iv.normalised(dat = datpooled_all$test.pool.normalised,
                      nlastpoints = nlastpoints)
  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  nca.results <- data.frame(
    cl = signif(nca.output[1], 3),
    vd = signif(nca.output[2], 3),
    slope = signif(nca.output[3], 3),
    half_life = signif(nca.output[3], 3),
    start.time = start.time,
    time.spent = time.spent
  )
  
  # Ncacalplot(dat = dat,datpooled = datpooled_all,ncasubfd = F,nlastpoints=nlastpoints)
  
  return(list(nca.fd.results = nca.fd.results,
              nca.results = nca.results))
  
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
                              nlastpoints) {
  if (missing(nlastpoints)) {
    nlastpoints <- 4
  }
  
  trap.rule <-
    function(x, y)
      sum(diff(x) * (y[-1] + y[-length(y)])) / 2
  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"
  
  auct <- trap.rule(dat$TIME, dat$DV)
  
  # Select last 4/specified number for slope calculation
  temp1 <- tail(dat, n = nlastpoints)
  
  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)
  
  slope <- summary(abc)[[4]][[2]]
  
  ke <- -slope
  
  half_life <- 0.693 / ke
  
  C_last <- tail(temp1$DV, 1)
  
  auct_inf <- C_last / ke
  
  auc0_inf <- auct + auct_inf
  
  clnormalised <- 1 / auc0_inf
  
  vdnormalised <- clnormalised / ke
  
  cl <- clnormalised
  
  vd <- vdnormalised
  
  return(c(cl, vd, slope, half_life))
}
