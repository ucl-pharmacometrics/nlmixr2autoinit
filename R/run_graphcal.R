#' Run graphical calculation for intravenous data
#'
#' Performs a graphical calculation of clearance and volume of distribution based on provided pharmacokinetic data with concentration normalised by dose
#'
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param fdobsflag A flag indicating whether the dataset includes first-dose observations. Set to 1 if the data includes first-dose observations, otherwise set to 0.
#' @param nbins The number of bins for time binning.
#' @param nlastpoints The number of last points to be used for slope linear regression.
#' @return A data frame with calculated clearance, volume of distribution, slope, and time spent for the calculation.
#' @importFrom dplyr %>% mutate if_else group_by ungroup
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPT
#' dat <- nmpkconvert(dat)
#' dat <- calculate_tad(dat)
#' run_graphcal_iv(dat, fdobsflag = 1, nbins = 8, nlastpoints = 4)
#' @export

run_graphcal_iv <- function(dat,
                         fdobsflag,
                         nbins,
                         nlastpoints) {
  graph.fd.results <- data.frame(
    cl = NA,
    vd = NA,
    slope = NA,
    time.spent = 0
  )

  ####################Calculation only based on pooled data of first dose#############
  if (fdobsflag == 1) {
    start.time <- Sys.time()
    dat$DVnor <- dat$DV / dat$dose
    dat_fd <- dat[dat$dose_number == 1,]
    datpooled_fd <- pk.time.binning(testdat = dat_fd,
                                    nbins = nbins)

    graph.fd.output <-
      graphcal_iv(dat = datpooled_fd$test.pool.normalised,
               nlastpoints = nlastpoints)

    end.time <- Sys.time()
    time.spent <- round(difftime(end.time, start.time), 4)

    graph.fd.results <- data.frame(
      cl = signif(graph.fd.output[1], 3),
      vd = signif(graph.fd.output[2], 3),
      slope = signif(graph.fd.output[3], 3),
      time.spent = time.spent
    )
  }

  return(graph.fd.results)

}


#' Graphical calculation of clearance and volume of distribution
#'
#' Graphic calculation for clearance and volume of distribution from the provided data.
#' @param dat A data frame with intravenous pharmacokinetic data, at least containing TIME and DV these two columns
#' @param nlastpoints Number of last points to use for the linear regression of terminal slope
#' @return A named vector containing the calculated clearance (cl), volume of distribution (vd), slope of terminal phase (slope), and concentration at the time 0 (C0).
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(12, 8, 5, 3, 2, 1.5, 1 ))
#' graphcal_iv(dat, nlastpoints = 3)
#' @export
#'

graphcal_iv <- function(dat,
                     nlastpoints) {
  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  temp1 <- tail(dat, n = nlastpoints)

  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)

  slope <- summary(abc)[[4]][[2]]

  Intercept <- summary(abc)[[4]][[1]]

  kel <- -slope

  C0 <- exp(Intercept)

  vd <- 1 / C0

  cl = kel * vd

  return(c(
    cl = cl,
    vd = vd,
    slope = slope,
    C0 = C0
  ))
}

#' Graphical calculation of clearance and volume of distribution for oral case
#'
#' Graphic calculation for clearance and volume of distribution from the provided data.
#' @param dat A data frame with intravenous pharmacokinetic data, at least containing TIME and DV these two columns
#' @param nlastpoints Number of last points to use for the linear regression of terminal slope
#' @return A named vector containing the calculated clearance (cl), volume of distribution (vd), slope of terminal phase (slope), and concentration at the time 0 (C0).
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(1, 2, 5, 3, 2, 1.5, 1 ))
#' graphcal_oral(dat, nlastpoints = 3)
#' @export
#'

graphcal_oral <- function(dat,
                          nlastpoints) {
  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  temp1 <- tail(dat, n = nlastpoints)
  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)
  slope <- summary(abc)[[4]][[2]]
  Intercept <- summary(abc)[[4]][[1]]
  kel <- -slope
  C0 <- exp(Intercept)
  vd <- 1 / C0
  cl = kel * vd
  # Identify Cmax point
  Cmax_point <- which(dat$DV == max(dat$DV), arr.ind = T)

  absorb_phase <- dat[1:Cmax_point, ]
  absorb_phase <- absorb_phase[absorb_phase$TIME > 0, ]

  # At least 2 points for ka slope regression
  # Cmax, and C0 are also accepted
  if (nrow(absorb_phase) < 2) {
    ka = NA
  }

  if (nrow(absorb_phase) > 1) {
    absorb_phase$IVconc <- -kel * absorb_phase$TIME + Intercept
    absorb_phase$residuals <-
      absorb_phase$IVconc - log(absorb_phase$DV)

    abslinear <- lm(absorb_phase$residuals ~ absorb_phase$TIME)
    slope_ka <- summary(abslinear)[[4]][[2]]
    ka =  -slope_ka

  }

  return(c(
    ka = ka,
    cl = cl,
    vd = vd,
    slope = slope,
    C0 = C0
  ))
}
