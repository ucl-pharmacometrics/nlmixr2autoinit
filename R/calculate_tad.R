#' Calculate time after dose (TAD) for a pharmacokinetic data
#'
#' Calculates the time after dose (TAD) for each observation in a pharmacokinetic dataset. The TAD is calculated as the time difference between each observation and the most recent dose administration. The function handles both bolus and infusion dosing scenarios, adjusting the calculations accordingly based on the infusion_flag parameter.
#' @param dat A data frame containing the pharmacokinetic data.
#' @return A data frame with the calculated tad (time after the last dose).
#' @import nlmixr2
#' @import dplyr
#' @importFrom tidyr fill
#' @examples
#'
#' Example 1:
#' dat <- Bolus_1CPT
#' calculate_tad(dat)
#'
#' Example 2:
#' dat <- Infusion_1CPT
#' dat<-calculate_tad(dat)
#' dat
#' @export

calculate_tad <- function(dat) {

  if (!"dose_number" %in% colnames(dat)){
    cat("Calculate dose number...............")
    dat <- mark_dose_number(dat)
  }

  if (!"RATE" %in% colnames(dat)){
    dat$RATE<-NA
  }

  if (!"route" %in% colnames(dat)){
     dat$route = NA
  }

  dat <- dat %>% mutate(tad = NA_real_)
  # Identify concentration rows
  conc_rows <- dat$EVID == 0
  # Filter the dose lines
  dat <- dat %>%
    arrange(ID, resetflag, TIME) %>%
    group_by(ID,resetflag) %>%
    mutate(
      last_dose_time = if_else(EVID %in% c(1, 101, 4), TIME, NA_real_),
      last_dose_number = if_else(EVID %in% c(1, 101, 4), dose_number, NA_integer_),
      last_dose = if_else(EVID %in% c(1, 101, 4), AMT, NA_integer_),
      last_rate = if_else(EVID %in% c(1, 101, 4), RATE, NA_real_),  # Collecting rate information
      last_route = if_else(EVID %in% c(1, 101, 4), route, NA_character_)  # Collecting route information
    ) %>%

    fill(last_dose_time, last_dose_number, last_dose, last_rate,  last_route , .direction = "down") %>%
    ungroup()

  # Calculate TAD for each concentration
    dat <- dat %>%
      mutate(
        tad = if_else(conc_rows, TIME - last_dose_time, NA_real_),
        # if conc row, fill by last_dose_number, else, dose_number
        dose_number = if_else(conc_rows, last_dose_number, dose_number),
        dose = if_else(conc_rows, last_dose, AMT),
        rateobs= if_else(conc_rows, last_rate, RATE),
        routeobs=if_else(conc_rows, last_route, route)
      ) %>%

      select(-last_dose_time, -last_dose_number,-last_dose, -last_rate, - last_route  )

  return(dat)
}


#' Mark dose number
#'
#' Mark the dose number in the dataset by assigning a sequential dose number identified by EVID to each dose event, process dose number within the same resetflag and ID group.
#' @param dat A data frame containing the pharmacokinetic data. The data frame must include columns for ID, TIME, EVID, and AMT.
#' @return A data frame with an additional column `dose_number` indicating the sequential number of each dose event for each subject.
#' @import nlmixr2
#' @importFrom dplyr %>% group_by arrange mutate if_else ungroup
#' @examples
#' dat <-Bolus_1CPT
#' mark_dose_number(dat)
#' @export
#'
mark_dose_number <- function(dat) {

  if (!"resetflag" %in% colnames(dat)) {
    dat$resetflag<- 1
  }

  # dat <- dat[with(dat, order(ID, resetflag,TIME, -AMT)), ]
  dat <- dat %>%
    group_by(ID, resetflag) %>%
    arrange(TIME) %>%
    mutate(dose_number = if_else(EVID %in% c(1, 101, 4),
                                 cumsum(EVID %in% c(1, 101, 4)),
                                 NA_integer_)) %>%
    ungroup()

  dat <- dat[with(dat, order(ID, resetflag,TIME, -AMT)), ]
  return(dat)
}
