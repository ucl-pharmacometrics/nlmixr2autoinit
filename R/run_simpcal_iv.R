#' Run approximate pharmacokinetic calculation for intravenous pharmacokinetic data
#'
#' Perform a simplified calculation for clearance and volume of distribution in intravenous cases.
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param route A character string specifying the route of administration (e.g, "bolus" and "infusion"). The default value is "bolus".
#' @param sdflag A flag indicating whether this is only single-dose data.
#' @param fdobsflag A flag indicating whether first dose observations are available.
#' @param half_life The half-life of the drug.
#' @return A list containing the results of the simplified calculation, including clearance, volume of distribution, and the processed dataset.
#' @importFrom dplyr %>% mutate if_else arrange group_by ungroup slice_min filter select
#' @import nlmixr2
#' @examples
#'
#' # Example 1 (iv case)
#' dat <- Bolus_1CPT
#' fdat <- processData(dat)
#' half_life<-half_life_estimated(fdat = fdat)[1]
#' run_simpcal_iv(fdat,half_life)
#'
#' # Example 2 (infusion case).
#'
#' dat <- Infusion_1CPT
#' fdat <- processData(dat)
#' half_life<-half_life_estimated(fdat = fdat)[1]
#' run_simpcal_iv(fdat,half_life)

#' @export

run_simpcal_iv <- function(dat,
                           half_life=NA) {

  start.time <- Sys.time()

  # Default
  median.simpcal.cl <- NA
  dat.ss.obs <- NA
  median.simpcal.vd <- NA
  dat.fd.obs <- NA

 ##################################Calculate of clearance######################

  dat <- is_ss(df = dat,
                 half_life = half_life )
  dat$duration_obs<-0

  if (nrow(dat[dat$rateobs!=0,]>0)){
  dat$duration_obs <- dat[dat$rateobs!=0,]$dose / dat[dat$rateobs!=0,]$rateobs
  }

  dat.ss.obs <- dat[dat$SteadyState == T,]

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

      # Calculate median cl for each individual
      individual_mean_cl <- aggregate(cl ~ ID, data = dat.ss.obs, FUN = mean)

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_cl <- mean(individual_mean_cl$cl, trim = 0.05, na.rm = TRUE)

   }

  ########################### Calculate the Volume of distribution###############

  if (unique(dat[dat$EVID==1,]$route) == "bolus") {

    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 & dat$tad < half_life * 0.2,]) > 0) {

      dat$C_first_flag <- 0
      dat[dat$EVID == 0 &
            dat$dose_number == 1 & dat$tad < half_life * 0.2,]$C_first_flag <- 1


    dat.fd.obs <- dat[dat$C_first_flag == 1, ]

    dat.fd.obs$vd <- signif(dat.fd.obs$dose / dat.fd.obs$DV, 3)

    # Calculate median vd for each individual
    individual_mean_vd <- aggregate(vd ~ ID, data = dat.fd.obs, FUN = mean)

    # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
    trimmed_mean_vd <- mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE)

    }
  }

  # Short-term infusion assumed
  if (unique(dat[dat$EVID==1,]$route) == "Infusion") {

    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 &
                 dat$tad < half_life * 0.2 & dat$TIME < dat$duration_obs,]) > 0) {

      dat$C_first_flag <- 0

      dat[dat$EVID == 0 &
            dat$dose_number == 1 &
            dat$tad < half_life * 0.2 &
            dat$TIME < dat$duration_obs,]$C_first_flag <- 1

      dat.fd.obs <- dat[dat$C_first_flag == 1, ]
      dat.fd.obs$vd <-
        signif(dat.fd.obs$TIME * dat.fd.obs$rateobs / dat.fd.obs$DV, 3)

      # Calculate median vd for each individual
      individual_mean_vd <- aggregate(vd ~ ID, data = dat.fd.obs, FUN = mean)

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_vd <- mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE)

    }
  }

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  # Only selected the key columns
  simpcal.results <- data.frame(
    cl = signif( trimmed_mean_cl, 3),
    vd = signif( trimmed_mean_vd, 3),
    starttime = Sys.time(),
    time.spent = time.spent
  )


  return(
    list(
      simpcal.results = simpcal.results,
      dat = dat,
      dat.ss.obs = dat.ss.obs,
      dat.fd.obs = dat.fd.obs
    )
  )

}
