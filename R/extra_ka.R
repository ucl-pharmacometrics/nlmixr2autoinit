
ka_wanger_nelson<-function(dat,nlastpoints){
  # Select last 4/specified number for slope calculation
  temp1 <- tail(dat, n = nlastpoints)
  C_last <- tail(temp1$DV, 1)
  # linear regression for slope of log of DVs
  abc <- lm(log(temp1$DV) ~ temp1$TIME)
  slope <- summary(abc)[[4]][[2]]
  ke<- -slope
  x=dat$TIME
  y<-dat$DV
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
  dat$frac_abs<-  (dat$DV+ ke*dat$auc_accumulate)/ ( (ke*sum(dat$auc_intervals)) + C_last / ke)
  max_abs_time <- which(dat$frac_abs == max(dat$frac_abs ), arr.ind = T)
  dat_absorb_phase <- dat[1:max_abs_time, ]
  dat_absorb_phase$frac_remained<-1-  dat_absorb_phase$frac_abs
  abc2 <- lm(log(dat_absorb_phase$frac_remained) ~  dat_absorb_phase$TIME)
  slope_ka<-summary(abc2)[[4]][[2]]
  ka= -slope_ka
  return(ka)
}

dat<-Oral_1CPT[Oral_1CPT$ID==1 &Oral_1CPT$SD==1&Oral_1CPT$EVID==0,]
dat<-data.frame(TIME=dat$TIME,DV=dat$DV)

ka_wanger_nelson(dat=dat,nlastpoints=4)


# trap.rule <- function(x, y){
#   sum(diff(x) * (y[-1] + y[-length(y)])) / 2
# }



# Function to calculate AUMC using the trapezoidal rule
calculate_aumc <- function(x, y, nlastpoints) {
  time<-x
  concentration<-y
  moment_curve <- time * concentration
  aumc0_t <- trap.rule(time, moment_curve)
  # Calculate AUMC from t_last to infinity
  # Extract t_last and C_last
  t_last <- tail(time, n = 1)
  C_last <- tail(concentration, n = 1)
  # Select last 4/specified number for slope calculation
  temp1 <- tail(dat, n = nlastpoints)
  # linear regression for slope of log of DVs
  abc<- lm(log(temp1$DV) ~ temp1$TIME)
  lambda_z <- - summary(abc)[[4]][[2]]
  half_life <- log(2) /  (lambda_z)
  # Calculate AUMC from t_last to infinity
  aumc_tlast_to_inf <- (C_last * t_last) / lambda_z + C_last / (lambda_z^2)
  aumc_total<- aumc0_t + aumc_tlast_to_inf
  return(c(aumc_total,half_life))
}


# Main function to calculate ka using the statistical moments method
calculate_ka_statistical_moments <- function(x, y, nlastpoints) {
  # Calculate AUC and AUMC for i.v. administration
  auc_ <- nca.iv.normalised(dat = data.frame(time=x,conc=y))[6]
  aumc_ <- calculate_aumc(x, y,nlastpoints)[1]
  mrt_oral <- aumc_ / auc_
  # Calculate MRT for i.v. by half life
  mrt_iv <-calculate_aumc(x, y,nlastpoints)[2]/log(2)
  # Calculate Mean Absorption Time (MAT)
  mat <- mrt_oral - mrt_iv

  # Estimate absorption rate constant (ka)
  ka <- 1 / mat

  return(ka)
}

dat<-Oral_1CPT[Oral_1CPT$ID==1 &Oral_1CPT$SD==1&Oral_1CPT$EVID==0,]
dat<-data.frame(TIME=dat$TIME,DV=dat$DV)

calculate_ka_statistical_moments(x =dat$TIME,y=dat$DV,nlastpoints = nlastpoints)


