
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
dat<-data.frame(TIME=dat$TIME,DV=dat$DV)
ka_wanger_nelson(dat=dat,nlastpoints=4)



