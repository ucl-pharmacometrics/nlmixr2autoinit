#' Perform single-point pharmacokinetic calculations
#'
#' Performs pharmacokinetic calculations using a single-point method and calculates clearance (CL) and volume of distribution (Vd) based on the provided pharmacokinetic data.
#'
#' @details
#' The function includes two main calculations:
#'
#' 1. **Clearance (CL):**
#'    - Calculated using steady-state data from multiple dosing regimens.
#'    - Steady-state concentrations are identified and used to estimate clearance based on the dose, observed steady-state concentrations, and dosing intervals.
#'
#' 2. **Volume of Distribution (Vd):**
#'    - Calculated for intravenous single-dose data.
#'    - The first sampling point (C0) after dose administration is used to estimate Vd with the formula `Dose / C0`.
#'
#' The function supports both bolus and short-term infusion routes. It processes input data to identify steady-state and first-dose conditions and calculates trimmed mean values for CL and Vd to reduce the impact of outliers.
#'
#' @param dat A data frame containing pharmacokinetic data, including information on doses, sampling times, and observed concentrations.
#' @param half_life The half-life of the drug, used for determining elimination phase parameters and steady-state conditions.
#' @return A list containing the following:
#'   - `singlepoint.results`: A data frame with the calculated clearance (CL) and volume of distribution (Vd), along with processing start and elapsed time.
#'   - `dat`: The processed dataset with additional columns for steady-state and other calculated flags.
#'   - `dat.ss.obs`: Subset of the data used for steady-state clearance calculation.
#'   - `dat.fd.obs`: Subset of the data used for single-dose volume of distribution calculation.
#' @importFrom dplyr %>% mutate if_else arrange group_by ungroup slice_min filter select
#' @import nlmixr2
#'
#' @examples
#' # Example 1 (iv case)
#' dat <- Bolus_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)[1]
#' run_single_point(fdat,half_life)
#'
#' # Example 2 (infusion case).
#'
#' dat <- Infusion_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)[1]
#' run_single_point(fdat,half_life)
#'
#' dat <- Oral_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median[1]
#' run_single_point(fdat,half_life)
#'
#'
#' @export

run_single_point <- function(dat,
                              half_life=NA) {
  start.time <- Sys.time()
  # Default
  trimmed_mean_ka <- NA
  trimmed_mean_cl <- NA
  trimmed_mean_vd <- NA
  ka_single_point.out<-NA
  dat.ss.obs <- NA
  dat.fd.obs <- NA
  cmax_by_group<-NA
  single_point.message<-NA

##############################Single Point Base#############################

 ##################################Calculate of clearance######################
  dat <- is_ss(df = dat,
                 half_life = half_life )


  dat.ss.obs <- dat[dat$SteadyState == TRUE & !is.na(dat$SteadyState), ]

  if (nrow(dat.ss.obs) > 2) {
      # If there are multiple points within the same dose interval,
      # only the minimum and maximum values are selected as the steady-state points for calculation
      dat.ss.obs <- dat.ss.obs %>%
        group_by(ID, dose_number) %>%
        mutate(
          max_value = max(DV),
          min_value = min(DV),
          avg_value = (max_value + min_value) / 2,
          max_interval = ifelse(DV == max(DV), TRUE, FALSE),
          min_interval = ifelse(DV == min(DV), TRUE, FALSE)
        ) %>%
        ungroup()

      # Extract rows to be marked in dat.ss.obs
      rows_to_mark <- dat.ss.obs %>%
        filter(max_interval == TRUE | min_interval == TRUE) %>%
        select(ID, dose_number, TIME)

      dat <- dat %>%
        rowwise() %>%
        mutate(SteadyState = ifelse(
          any(
            ID == rows_to_mark$ID &
              dose_number == rows_to_mark$dose_number &
              TIME == rows_to_mark$TIME
          ),
          TRUE,
          FALSE
        )) %>%
        ungroup()

      # Same for dat.ss.obs
      dat.ss.obs  <- dat.ss.obs  %>%
        rowwise() %>%
        mutate(SteadyState = ifelse(
          any(
            ID == rows_to_mark$ID &
              dose_number == rows_to_mark$dose_number &
              TIME == rows_to_mark$TIME
          ),
          TRUE,
          FALSE
        )) %>%
        ungroup()

      # Second selection， only select the max, min points
      dat.ss.obs <- dat.ss.obs[dat.ss.obs$SteadyState == T,]

      # Css type identification
      # Only applicable for oral case with very fast absorption, and with the majority of absorption occurring before Tmax, and only a short period is required to reach Tmax.

      dat.ss.obs  <-   dat.ss.obs %>%
        mutate(
          # Initialise 'Css_type' to "Css,avg" by default
          Css_type = "Css_avg",

          # Identify Cssmax and Cssmin based on 'tad', 'ii', and the condition max_value == min_value
          Css_type = case_when(
            max_value == min_value & tad <= 0.2 * recent_ii ~ "Css_max",            # tad within 20% of ii -> Css,min
            max_value == min_value & tad >= 0.8 * recent_ii & tad <= recent_ii ~ "Css_min", # tad within 80-100% of ii -> Css,max
            TRUE ~ Css_type                                                   # Default to Css,avg for all other cases
          )
        )


      dat.ss.obs$cl <- NA

      # Calculate the clearance by single point and Css_type

      for (i in 1:nrow(dat.ss.obs)) {
        # Only one points in the dose interval
        if (dat.ss.obs[i,]$Css_type=="Css_avg") {
          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / dat.ss.obs[i,]$avg_value /dat.ss.obs[i,]$recent_ii ,
              digits = 3
            )
        }


        if (dat.ss.obs[i,]$Css_type=="Css_max" & dat.ss.obs[i,]$routeobs=="bolus") {

          Css_min_i<-dat.ss.obs[i,]$max_value * exp(-((log(2)/as.numeric(half_life))*dat.ss.obs[i,]$recent_ii))
          Css_avg_i<-(dat.ss.obs[i,]$max_value+   Css_min_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )
        }

        if (dat.ss.obs[i,]$Css_type=="Css_max" & dat.ss.obs[i,]$routeobs=="infusion") {

          Css_min_i<-dat.ss.obs[i,]$max_value * exp(-((log(2)/as.numeric(half_life))*(dat.ss.obs[i,]$recent_ii-dat.ss.obs[i,]$duration_obs)))
          Css_avg_i<-(dat.ss.obs[i,]$max_value+   Css_min_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )
        }

        if (dat.ss.obs[i,]$Css_type=="Css_min" & dat.ss.obs[i,]$routeobs=="bolus") {

          Css_max_i<-dat.ss.obs[i,]$min_value / exp(-((log(2)/as.numeric(half_life))*dat.ss.obs[i,]$recent_ii))
          Css_avg_i<-(dat.ss.obs[i,]$min_value + Css_max_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )

        }


        if (dat.ss.obs[i,]$Css_type=="Css_min" & dat.ss.obs[i,]$routeobs=="infusion") {

          Css_max_i<-dat.ss.obs[i,]$min_value / exp(-((log(2)/as.numeric(half_life))*(dat.ss.obs[i,]$recent_ii-dat.ss.obs[i,]$duration_obs)))
          Css_avg_i<-(dat.ss.obs[i,]$min_value + Css_max_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )

        }

      } # end loop


      # Calculate geometric mean cl for each individual
      individual_mean_cl <- tryCatch(aggregate(cl ~ ID, data = dat.ss.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_cl <- tryCatch(trimmed_geom_mean(individual_mean_cl$cl, trim=0.05,na.rm = TRUE),error=function(e) {NA})

   }

  ########################### Calculate the Volume of distribution###############

  if (unique(dat[dat$EVID==1,]$route) == "bolus") {

    # remove ssflag=1 case
    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 & dat$tad < half_life * 0.2 & dat$iiobs==0,]) > 0) {

      dat$C_first_flag <- 0
      dat[dat$EVID == 0 &
            dat$dose_number == 1 & dat$tad < half_life * 0.2 & dat$iiobs==0,]$C_first_flag <- 1


    dat.fd.obs <- dat[dat$C_first_flag == 1, ]

    dat.fd.obs$vd <- signif(dat.fd.obs$dose / dat.fd.obs$DV, 3)

    # Calculate median vd for each individual
    individual_mean_vd <- aggregate(vd ~ ID, data = dat.fd.obs, FUN = trimmed_geom_mean)

    # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
    trimmed_mean_vd <-trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE)

    }
  }

  # Short-term infusion assumed
  if (unique(dat[dat$EVID==1,]$route) == "infusion") {
    # remove ssflag=1 case
    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 &
                 dat$tad < half_life * 0.2 &
                 # dat$TIME < dat$duration_obs & # allow sampling outside infusion duration
                 dat$iiobs==0,]) > 0) {

      dat$C_first_flag <- 0

      dat[dat$EVID == 0 &
            dat$dose_number == 1 &
            dat$tad < half_life * 0.2 &
            # dat$TIME < dat$duration_obs  &
            dat$iiobs==0,]$C_first_flag <- 1

      dat.fd.obs <- dat[dat$C_first_flag == 1, ]
      dat.fd.obs$vd <-
        signif(pmin(dat.fd.obs$TIME,dat.fd.obs$duration_obs) * dat.fd.obs$rateobs / dat.fd.obs$DV, 3)

      # Calculate median vd for each individual
      individual_mean_vd <- tryCatch( aggregate(vd ~ ID, data = dat.fd.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_vd <-tryCatch(trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE),error=function(e) {NA})

    }
  }

  ##############################Single Point Extra#############################
  # Both can be calculated in the base part
  if (is.na(trimmed_mean_cl)==F & is.na(trimmed_mean_vd)==F) {
    single_point.message<-"Clearance was calculated using steady-state data and volume of distribution was calculated based on data within the dosing interval following the first dose "
  }

 # if half_life is available, Vd is not available
  if (is.na(trimmed_mean_cl)==F & is.na(trimmed_mean_vd)==T & is.na(half_life)==F) {
    message(black(
      paste0("Insufficient single-dose (IV) data for Vd calculation; derived from clearance and estimated half-life instead.", strrep(".", 20))))

    individual_mean_vd <- tryCatch( aggregate(cl* half_life/log(2) ~ ID, data = dat.ss.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

    trimmed_mean_vd <-  signif(trimmed_mean_cl * half_life/log(2), 3)
    single_point.message<-"Clearance was calculated using steady-state data and volume of distribution was calculated based on clearance and estimated half_life"
  }

  # single-point method volume of distribution + half-life estimated
  if (is.na(trimmed_mean_cl)==T & is.na(trimmed_mean_vd)==F & is.na(half_life)==F) {
    message(black(
      paste0("Insufficient steady-state data for CL calculation; derived from Vd and estimated half-life instead.", strrep(".", 20))))

    trimmed_mean_cl <-  signif(trimmed_mean_vd * log(2) / half_life, 3)
    single_point.message<-"Vd was calculated based on data within the dosing interval following the first dose and clearance was calculated based on volume of distribution and estimated half_life"
  }

  if (is.na(trimmed_mean_cl)==T & is.na(trimmed_mean_vd)==T & is.na(half_life)==F) {
    message(black(
      paste0("Neither single-dose (IV) data nor steady-state data supports the calculation of clearance (CL) and volume of distribution (Vd). Vd will be estimated using Dose/Cmax, and CL will be derived based on the estimated half-life and Vd.", strrep(".", 20))))

    datobs<-dat[dat$EVID==0,]

    cmax_by_group <- datobs %>%
      group_by(ID, dose_number) %>%
      slice_max(order_by = DV, with_ties = FALSE) %>%
      ungroup()

    cmax_by_group<-cmax_by_group[cmax_by_group$tad<half_life*0.2,]

    if (unique(dat[dat$EVID==1,]$route) == "infusion") {
    cmax_by_group$vd <-
       signif(pmin(cmax_by_group$TIME,cmax_by_group$duration_obs) * cmax_by_group$rateobs / cmax_by_group$DV, 3)
    }

    cmax_by_group$vd <-
      signif(cmax_by_group$dose/ cmax_by_group$DV, 3)

    # Calculate median vd for each individual
    individual_mean_vd <- tryCatch( aggregate(vd ~ ID, data = cmax_by_group, FUN = trimmed_geom_mean),error=function(e) {NA})

    # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
    trimmed_mean_vd <-tryCatch(trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE),error=function(e) {NA})


    trimmed_mean_cl <-  signif(trimmed_mean_vd * log(2) / half_life, 3)
    single_point.message<-"Vd was calculated based on Cmax in each dose interval and clearance was calculated based on volume of distribution and estimated half_life"
  }

  # extra ka for oral case
  if (unique(dat[dat$EVID==1,]$route) == "oral") {

    if (is.na( trimmed_mean_cl)==F & is.na( trimmed_mean_vd)==F){

      datobs<-dat[dat$EVID==0,]

      cmax_by_group2 <- datobs %>%
        group_by(ID, dose_number) %>%
        mutate(Tmax = TIME[which.max(DV)]) %>%
        filter(TIME < Tmax | (n() == 1 & tad < 0.2*log(2)*trimmed_mean_cl/trimmed_mean_vd)) %>%
        slice_max(order_by = DV, with_ties = FALSE) %>%
        ungroup()%>%
        select(-Tmax)

      if (nrow(cmax_by_group2[cmax_by_group2$EVID==0  & cmax_by_group2$iiobs==0,])>0){
        ka_single_point.out<-run_ka_solution(df =cmax_by_group2,
                                        cl = trimmed_mean_cl,
                                        ke = trimmed_mean_cl/trimmed_mean_vd)

        trimmed_mean_ka<-  tryCatch(trimmed_geom_mean(ka_single_point.out$ka_calc_dat$ka_calcv, trim = 0.05, na.rm = TRUE),error=function(e) {NA})

      }
    }
  }

  end.time <- Sys.time()

  time.spent <- round(difftime(end.time, start.time), 4)

  # Only selected the key columns
  singlepoint.results <- data.frame(
    ka = signif( trimmed_mean_ka, 3),
    cl = signif( trimmed_mean_cl, 3),
    vd = signif( trimmed_mean_vd, 3),
    starttime = Sys.time(),
    time.spent = time.spent,
    single_point.message=single_point.message
  )


  return(
    list(
      singlepoint.results = singlepoint.results,
      dat = dat,
      single_point_ka_df=   ka_single_point.out,
      single_point_cl_df =  dat.ss.obs,
      single_point_vd_df =  dat.fd.obs,
      single_point_vd_cmax_df=cmax_by_group
    )
  )

}

#' Calculate the trimmed geometric mean
#'
#' Computes the trimmed geometric mean of a numeric vector.
#' A specified proportion of the smallest and largest values is removed before calculating the geometric mean.
#'
#' @param x A numeric vector containing the values for which the trimmed geometric mean is to be calculated.
#' @param trim A numeric value between 0 and 0.5 indicating the proportion of observations
#'   to be trimmed from each end of the vector. Default is 0 (no trimming).
#' @param na.rm Logical. Should missing values (NA) be removed before computation? Defaults to TRUE.
#' @return The trimmed geometric mean of the input vector as a numeric value.
#'   If no values remain after trimming, returns `NA`.
#'
trimmed_geom_mean <- function(x, trim = 0, na.rm = TRUE) {
  if (na.rm) x <- na.omit(x)  # Remove NA values if na.rm is TRUE

  # Sort the vector and apply trimming
  x <- sort(x)
  n <- length(x)
  lower <- floor(n * trim) + 1  # Index of first value to keep
  upper <- n - floor(n * trim)  # Index of last value to keep

  # Keep only trimmed range
  x_trimmed <- x[lower:upper]

  # Calculate the geometric mean of the trimmed values
  exp(mean(log(x_trimmed), na.rm = FALSE))
}
