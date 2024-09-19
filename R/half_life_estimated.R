
half_life_estimated<-function(dat,
                              sdflag,
                              fdobsflag,
                              nlastpoints){

# split the data into three types
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
  half_life_fd <- 0.693 / ke_fd
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

# data with samples after the first dose and repeated dose
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
half_life_mean <- exp(mean(log(positive_values)))

return(half_life_mean)

}
