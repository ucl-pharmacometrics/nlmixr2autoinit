#' Estimate half-life from given pharmacokinetic Data
#'
#' Estimates the half-life based on pharmacokinetic data by linear regression on terminal phase,
#' considering different scenarios: data after the first dose, data after repeated doses,
#' and data including both first and repeated doses.
#'
#' @param fdat A data frame containing the pharmacokinetic data.
#' @param nlastpoints An integer specifying the number of last data points to be used for
#'   linear regression in the estimation.
#'
#' @return A named vector with the estimated half-life values from each data subset, including:
#'   \itemize{
#'     \item `half_life_median`: The median of all positive half-life estimates from the subsets.
#'     \item `half_life_fd`: The half-life estimate from first dose data, if available.
#'     \item `half_life_md`: The half-life estimate from repeated dose data, if available.
#'     \item `half_life_all`: The half-life estimate from the complete dataset.
#'   }
#' @details
#' The function performs the following steps:
#'   1. Splits the input data into subsets based on dosing intervals.
#'   2. Applies time binning to normalize the concentration-time data in each subset.
#'   3. Performs linear regression on the log-transformed concentration of the last `nlastpoints`
#'      to calculate the elimination rate constant (`ke`) and the half-life (`half_life`).
#'   4. Estimates the half-life separately for:
#'      \itemize{
#'        \item First dose interval data
#'        \item Repeated dose interval data
#'        \item All data combined
#'      }
#'   5. Returns the median of the positive half-life estimates as `half_life_median`.
#' @examples
#' # Example usage:
#' dat <- Bolus_1CPT
#' dat<- processData(dat)$dat
#' half_life_estimated(fdat)
#'
#' @importFrom stats lm
#' @export
#'


half_life_estimated<-function(dat,
                              nlastpoints=3,
                              nbins=10,
                              route="bolus"
                              ){

  message(black(
    paste0("Estimating half-life",strrep(".", 20))))

  #Preset dataset
  datpooled_all<-NA
  datpooled_fd<-NA
  datpooled_md<-NA

  #Predefine half life obtained from three peroids of data
  half_life_fd <-NA
  half_life_md <-NA
  half_life_all <-NA

  datpooled_all <- pk.time.binning(testdat = dat,
                                   nbins = nbins)

  half_life_all<-get_hf(datpooled_all$test.pool.normalised,
                        nlastpoints = nlastpoints,
                        route=route)


  if (nrow(dat[dat$dose_number==1 & dat$EVID==0 & dat$iiobs==0,])>0){
    datpooled_fd <- pk.time.binning(testdat = dat[dat$dose_number==1 & dat$iiobs==0,],
                                    nbins = nbins)

    half_life_fd<-get_hf(datpooled_fd$test.pool.normalised,nlastpoints = nlastpoints,route = route)
  }


  if (nrow(dat[dat$dose_number>1 & dat$EVID==0 & dat$iiobs==0,])>0){

    datpart1<- dat[dat$dose_number>1,]
    datpart2<- dat[dat$dose_number==1 & dat$iiobs>0,]

    testdat<-rbind(datpart1,datpart2)
    testdat<- testdat[with(testdat, order(ID, resetflag, TIME, -AMT)), ]

    datpooled_md <- pk.time.binning(testdat = testdat,
                                    nbins = nbins)

    half_life_md<-get_hf(datpooled_md$test.pool.normalised,
                         nlastpoints = nlastpoints,
                         route=route)
    }


half_life_values<-c( half_life_fd[1] ,half_life_md[1] , half_life_all[1])

# Remove negative numbers
positive_values <-  half_life_values[ half_life_values > 0 &  !is.na( half_life_values)]
# Calculate geometric mean of the positive values
# half_life_geom_mean <- round(exp(mean(log(positive_values))),2)
half_life_median <- round(median(positive_values),2)

return(list(half_life_median=half_life_median,
         half_life_fd=  half_life_fd,
         half_life_md=  half_life_md,
         half_life_all=  half_life_all))

}


#' Calculate half-life from the slope of log-transformed concentration data
#'
#' Estimates the half-life of a substance using the slope of log-transformed
#' concentration data (`Conc`) over time (`Time`). It performs a linear regression on
#' the last `nlastpoints` of normalized concentration data to determine the elimination rate
#' constant (`ke`) and half-life (`half_life`).
#'
#' @param testdat A data frame containing at least two columns:
#'   - `Conc`: Concentration values at each time point.
#'   - `Time`: Corresponding time points for each concentration measurement.
#' @param nlastpoints Integer. The number of last data points to use for regression analysis.
#'
#' @return A numeric value representing the estimated half-life (`half_life_`) if the slope is negative.
#'   Returns `NA` if the slope is non-negative or if the dataset has fewer than `nlastpoints`.
#'
#' @details
#' The function performs the following steps:
#'   1. Identifies Tmax as the time point with the maximum concentration (`Conc`).
#'   2. For "bolus" and "infusion" routes, includes all points from Tmax onwards.
#'      For the "oral" route, only points after Tmax are included.
#'   3. Adjusts `nlastpoints` if there are fewer data points than specified. There are at least two points for linear regression.
#'   4. Performs a linear regression on `log(Conc)` against `Time`.
#'   5. If the slope from the regression is negative, calculates the elimination rate
#'      constant (`ke`) as the negative slope and the half-life as `log(2) / ke`.
#'
#' @examples
#' # Example usage:
#'
#'  # Test first ID in the case of Bolus_1CPT
#'  testdat= Bolus_1CPT[Bolus_1CPT$SD==1 & Bolus_1CPT$ID==1 & Bolus_1CPT$EVID==0,]
#'  testdat=subset(testdat,select=c(TIME,DV))
#'  colnames(testdat)<-c("Time","Conc")
#'  get_hf(testdat = testdat, nlastpoints = 3)
#'
#'  # Test first ID in the case of Oral_1CPT (sparse data)
#'  dat<-Oral_1CPT[Oral_1CPT$ID==1  &Oral_1CPT$EVID==0,]
#'  datobs.sparse<-datobs[datobs$TIME %in% c(146,164,167.99),]
#'  testdat=subset(datobs.sparse,select=c(TIME,DV))
#'  colnames(testdat)<-c("Time","Conc")
#'  get_hf(testdat = testdat, nlastpoints = 3,route="oral")
#'
#'
#' @export

get_hf <- function(testdat,
                   nlastpoints = 3,
                   route = "bolus") {
  half_life_ <- NA
  # Identify the index of Tmax
  max_index <- which.max(testdat$Conc)

  temp1<-testdat

  if (route == "infusion") {
    temp1 <-
      testdat[max_index:nrow(testdat),] # Subset data points after Tmax
  }

  if (route == "oral") {
    if (max_index==nrow(testdat)){
      return(half_life_)
    }
    temp1 <-
      testdat[(max_index + 1):nrow(testdat),] # Subset data points after Tmax
  }

  # Loop to adjust nlastpoints if there are insufficient points
  while (nlastpoints > 1 && nrow(temp1) < nlastpoints) {
    nlastpoints <- nlastpoints - 1
  }

  # If fewer than 2 points are available after the loop, exit as calculation is invalid
  if (nlastpoints < 2) {
    # message("Insufficient data points for valid regression.")
    return(half_life_)
  }

  # Select the required last nlastpoints data points
  temp1 <- tail(temp1, n = nlastpoints)

  # linear regression for slope of log of DVs
  funcslope <- lm(log(temp1$Conc) ~ temp1$Time)
  slope_ <- summary(funcslope)[[4]][[2]]

  if (slope_ < 0) {
    ke_ <-  -slope_
    half_life_ <- log(2) / ke_ # change to log(2)
  }

  return(c(half_life_=half_life_,nlastpoints=nlastpoints))
}











