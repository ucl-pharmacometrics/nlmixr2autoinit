#' Run simplified calculation for intravenous pharmacokinetic data
#'
#' Perform a simplified calculation for clearance and volume of distribution in intravenous cases.
#' @param dat A data frame containing the intravenous pharmacokinetic data.
#' @param infusion_flag A flag indicating whether the dataset includes infusion doses. Set to 1 if the data includes infusions, otherwise set to 0.
#' @param sdflag A flag indicating whether this is only single-dose data.
#' @param fdobsflag A flag indicating whether first dose observations are available.
#' @param half_life The half-life of the drug.
#' @return A list containing the results of the simplified calculation, including clearance, volume of distribution, and the processed dataset.
#' @importFrom dplyr %>% mutate if_else arrange group_by ungroup slice_min filter select
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPT
#' dat<- nmpkconvert(dat)
#' dat<- calculate_tad(dat)
#' run_simpcal_iv(dat, infusion_flag = 0, sdflag = 0, fdobsflag = 1)
#' @export

run_simpcal_iv <- function(dat,
                        infusion_flag,
                        sdflag,
                        fdobsflag,
                        half_life) {
  start.time <- Sys.time()
  if (missing(half_life)) {
    half_life <- NA
  }

  ##################################Calculate of clearance######################
  median.simpcal.cl <- NA
  dat.ss.obs <- NA
  median.simpcal.vd <- NA
  dat.fd.obs <- NA

  if (sdflag == 0) {
    # Calculate the most commonly used dose interval
    dose_data <- dat[dat$EVID %in% c(1, 4, 101) & dat$AMT > 0, ]
    dose_data <- dose_data[order(dose_data$ID, dose_data$TIME), ]
    dose_data$interval <-
      ave(
        dose_data$TIME,
        dose_data$ID,
        FUN = function(x)
          c(NA, diff(x))
      )
    dose_intervals <-
      round(dose_data$interval[!is.na(dose_data$interval)], 0)
    most_commonly_used_dose_interval <-
      as.numeric(names(sort(table(dose_intervals), decreasing = TRUE)[1]))


    if (exists("half_life")) {
      if (is.na(half_life)) {
        half_life = most_commonly_used_dose_interval
        cat(
          "Warning: Half-life was not available due to some reasons. Change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered."
        )
      }
    }

    if (is.null(half_life)) {
      half_life = most_commonly_used_dose_interval
      cat(
        "Warning: Half-life was not available due to some reasons. Change the condition for reaching steady state after multiple dosing to whether 5 doses have been administered."
      )
    }

    # Identify potential points at the steady state.
    dat <- is_ss(df = dat,
                 half_life = half_life ,
                 dose_interval = most_commonly_used_dose_interval)

    # First selection, identify points during the steady-state phase.
    dat.ss.obs <- dat[dat$SteadyState == T,]

    if (nrow(dat.ss.obs) > 0) {
      # If there are multiple points within the same dose interval, only the minimum and maximum values are selected as the steady-state points for calculation
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
      dat.ss.obs$cl <- NA
      # Calculate the individual clearance.
      for (i in 1:nrow(dat.ss.obs)) {
        # Only one points in the dose interval
        if (dat.ss.obs[i,]$max_value == dat.ss.obs[i,]$min_value) {
          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / dat.ss.obs[i,]$DV / most_commonly_used_dose_interval,
              digits = 3
            )
        }
        # Two points
        if (dat.ss.obs[i,]$max_value > dat.ss.obs[i,]$min_value) {
          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / dat.ss.obs[i,]$avg_value / most_commonly_used_dose_interval,
              digits = 3
            )
        }
      }

      median.simpcal.cl <- median(dat.ss.obs$cl)
    }

    # No selected points
    if (nrow(dat.ss.obs) == 0) {
      median.simpcal.cl <- NA
    }

  }

  ########################### Calculate the Volume of distribution###############
  if (fdobsflag == 1) {
    dat$fdflag <- 0
    fp_tad <-
      3 # currently, it target the first plasma points within three hours, could make it variable in the future if needed.
    dat[dat$EVID == 0 &
          dat$dose_number == 1 & dat$tad < fp_tad, ]$fdflag <- 1
    dat.fd.obs <- dat[dat$fdflag == 1,]
    dat.fd.obs <- dat.fd.obs %>%
      group_by(ID) %>%
      slice_min(order_by = TIME, n = 1) %>%
      ungroup()


    # Extract rows to be marked in dat.ss.obs
    dat <- dat %>%
      rowwise() %>%
      mutate(fdflag = ifelse(
        any(
          ID == dat.fd.obs$ID &
            dose_number == dat.fd.obs$dose_number &
            TIME == dat.fd.obs$TIME
        ),
        1,
        0
      )) %>%
      ungroup()

    # Bolus case calculation
    if (infusion_flag == 0) {
      dat.fd.obs$vd <- signif(dat.fd.obs$dose / dat.fd.obs$DV, 3)
    }

    if (infusion_flag == 1) {
      # If less than infusion time, then it needs to multiply the time
      dat.fd.obs <- dat.fd.obs %>%
        group_by(ID) %>%
        mutate(vd = if_else(
          TIME < dose / RATE,
          signif(dose * TIME / DV, 3),
          signif(dose / DV, 3)
        )) %>%
        ungroup()


    }
    # if no enough points for statistics, even no available points found
    if (nrow(dat.fd.obs) < 3) {
      median.simpcal.vd <- NA
    }
    if (nrow(dat.fd.obs) > 2) {
      median.simpcal.vd <- median(dat.fd.obs$vd)
    }

  }

  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  # Only selected the key columns
  simpcal.results <- data.frame(
    cl = signif(median.simpcal.cl, 3),
    vd = signif(median.simpcal.vd, 3),
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
