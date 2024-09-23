#' Estimate half-life from given pharmacokinetic Data
#'
#' Estimates the half-life based on pharmacokinetic data by linear regression on terminal phase,
#' considering different scenarios: data after the first dose, data after repeated doses,
#' and data including both first and repeated doses.
#'
#' @param dat A data frame containing the pharmacokinetic data. It should include columns
#'   for `DV` (dependent variable), `dose`, and `dose_number`.
#' @param sdflag A numeric flag indicating if the data includes repeated dose samples.
#'   `0` means repeated dose data is included, `1` otherwise.
#' @param fdobsflag A numeric flag indicating if the data includes only the first dose samples.
#'   `1` means only first dose data is included, `0` otherwise.
#' @param nlastpoints An integer specifying the number of last data points to be used for
#'   linear regression in the estimation.
#'
#' @return A numeric value representing the geometric mean of the estimated half-life
#'   values from the different data subsets (first dose, repeated doses, and all data).
#' @details The function splits the input data into subsets based on the dose number,
#'   performs time binning, and estimates the half-life using the slope of the log-linear
#'   regression on the last points of the normalized concentration-time data.
#'
#' The half-life is estimated separately for:
#' \itemize{
#'   \item Only first dose data
#'   \item Only repeated dose data
#'   \item All data combined
#' }
#' The final result is the geometric mean of the positive half-life estimates from these
#' subsets.
#' @examples
#' # Example usage:
#' dat <- Bolus_1CPT
#' dat <- nmpkconvert(dat)
#' dat <- calculate_tad(dat)
#' half_life_estimated(dat, sdflag = 0, fdobsflag = 1, nlastpoints = 4,nbins=8)
#'
#' @importFrom stats lm
#' @export
#'

half_life_estimated<-function(dat,
                              sdflag,
                              fdobsflag,
                              nlastpoints,
                              nbins){

# Split the data into three types
dat$DVnor <- dat$DV / dat$dose
half_life_fd<-NA
half_life_efd<-NA
half_life_all<-NA

if (fdobsflag==1){
  # Data with samples only after first dose
  dat_fd <- dat[dat$dose_number == 1, ]
  datpooled_fd <- pk.time.binning(testdat = dat_fd,
                                  nbins = nbins)
  # half_life_estimated
  temp1 <- tail( datpooled_fd$test.pool.normalised, n = nlastpoints)
  # linear regression for slope of log of DVs
  fslope <- lm(log(temp1$Conc) ~ temp1$Time)
  slope_fd <- summary(fslope)[[4]][[2]]
  ke_fd <-  - slope_fd
  half_life_fd <- 0.693 / ke_fd # change to log(2)
}

if (sdflag == 0){
  # Data with samples only after repeated dose
  dat_efd <- dat[dat$dose_number != 1, ]
  datpooled_efd <- pk.time.binning(testdat = dat_efd,
                                   nbins = nbins)
  # half_life_estimated
  temp1 <- tail( datpooled_efd$test.pool.normalised, n = nlastpoints)
  # linear regression for slope of log of DVs
  fslope <- lm(log(temp1$Conc) ~ temp1$Time)
  slope_efd <- summary(fslope)[[4]][[2]]
  ke_efd <-  - slope_efd
  half_life_efd <- 0.693 / ke_efd
}

# Data with samples after the first dose and repeated dose
if (fdobsflag==1 & sdflag==0){
  datpooled_all <- pk.time.binning(testdat = dat,
                                   nbins = nbins)
  # half_life_estimated
  temp1 <- tail( datpooled_all$test.pool.normalised, n = nlastpoints)
  # linear regression for slope of log of DVs
  fslope <- lm(log(temp1$Conc) ~ temp1$Time)
  slope_all <- summary(fslope)[[4]][[2]]
  ke_all <-  - slope_all
  half_life_all <- 0.693 / ke_all
}



half_life_values<-c( half_life_fd ,half_life_efd , half_life_all)

# Remove negative numbers
positive_values <-  half_life_values[ half_life_values > 0 &  !is.na( half_life_values)]

# Calculate geometric mean of the positive values
half_life_mean <- round(exp(mean(log(positive_values))),2)

return(half_life_mean)

}
