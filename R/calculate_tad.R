#' Calculate time after dose for pharmacokinetic data
#'
#' Calculate time after dose (TAD) for pharmacokinetic observations.
#'
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @return A modified data frame with added columns:
#'   - tad: time after dose, calculated as the observation time minus the time
#'     of the most recent prior dose; set to NA for dosing records
#'   - iiobs: interdose interval inherited from the most recent dosing record
#'   - rateobs: infusion rate inherited from the most recent dosing record
#'   - routeobs (optional): route of administration inherited from the most
#'     recent dosing record, included only if route information is present
#'   - dose_number: sequential dose number, generated if not already present
#'
#' @details
#' The procedure identifies dosing events based on the event identifier (EVID)
#' and assigns each observation the attributes of the most recent prior dose.
#' The time after dose is then calculated for observation rows. If `dose_number`
#' column is not present in the input, it is automatically created for each subject.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' calculate_tad(Bolus_1CPT)
#' calculate_tad(Infusion_1CPT)
#' calculate_tad(Oral_1CPT)
#'
#' @export
#'
calculate_tad <- function(dat) {

  if ("DOSE" %in% colnames(dat)) {
    if ("raw_dose" %in% colnames(dat)) {
      dat$raw_dose <- NULL
    }
    dat <- dplyr::rename(dat, raw_dose = DOSE)
  }

  if ("dose" %in% colnames(dat)) {
    if ("raw_dose" %in% colnames(dat)) {
      dat$raw_dose <- NULL
    }
    dat <- dplyr::rename(dat, raw_dose = dose)
  }

  if (!"dose_number" %in% colnames(dat)) {
    cat("Calculate dose number...............")
    dat <- mark_dose_number(dat)
  }

  if (!"RATE" %in% colnames(dat)) {
    dat$RATE <- 0
  }

  if (!"II" %in% colnames(dat)) {
    dat$II <- 0
  }

  if (!"route" %in% colnames(dat)) {
    dat$route = NA
  }

  dat <- dat %>% dplyr::mutate(tad = NA_real_)
  # Identify concentration rows
  conc_rows <- dat$EVID == 0
  # Filter the dose lines
  dat <- dat %>%
    dplyr::arrange(ID, resetflag, TIME) %>%
    dplyr::group_by(ID, resetflag) %>%
    dplyr::mutate(
      last_dose_time = dplyr::if_else(EVID %in% c(1, 101, 4), TIME, NA_real_),
      last_dose_number = dplyr::if_else(EVID %in% c(1, 101, 4), dose_number, NA_integer_),
      last_dose = dplyr::if_else(EVID %in% c(1, 101, 4), AMT, NA_integer_),
      last_ii = dplyr::if_else(EVID %in% c(1, 101, 4), II, NA_real_),
      last_rate = dplyr::if_else(EVID %in% c(1, 101, 4), RATE, NA_real_),
      last_route = dplyr::if_else(EVID %in% c(1, 101, 4), route, NA_character_)
    ) %>%
    tidyr::fill(
      last_dose_time,
      last_dose_number,
      last_dose,
      last_ii,
      last_rate,
      last_route,
      .direction = "down"
    ) %>%
    dplyr::ungroup()

  # Calculate TAD for each concentration
  dat <- dat %>%
    dplyr::mutate(
      tad = dplyr::if_else(conc_rows, TIME - last_dose_time, NA_real_),
      # if conc row, fill by last_dose_number, else, dose_number
      dose_number = dplyr::if_else(conc_rows, last_dose_number, dose_number),
      dose = dplyr::if_else(conc_rows, last_dose, AMT),
      iiobs = dplyr::if_else(conc_rows, last_ii, II),
      rateobs = dplyr::if_else(conc_rows, last_rate, RATE),
      routeobs = dplyr::if_else(conc_rows, last_route, route)
    ) %>%

    dplyr::select(-last_dose_time,-last_dose_number,-last_dose,-last_ii,-last_rate,-last_route)

  # Remove 'route' and 'routeobs' if 'route' column is all NA
  if (all(is.na(dat$route))) {
    dat <- dat %>% dplyr::select(-route, -routeobs)
  }

  return(dat)
}


#' Mark dose number
#'
#' Assigns sequential dose numbers based on dosing events (EVID) within each subject.
#'
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @return A modified data frame with an added column named dose_number, indicating
#'   the sequential dose count within each subject and reset group.
#'
#' @examples
#' mark_dose_number(Bolus_1CPT)
#' mark_dose_number(Infusion_1CPT)
#' mark_dose_number(Oral_1CPT)
#'
#' @author Zhonghui Huang
#'
#' @export
#'
mark_dose_number <- function(dat) {
  if (!"resetflag" %in% colnames(dat)) {
    dat$resetflag <- 1
  }

  dat <- dat %>%
    dplyr::group_by(ID, resetflag) %>%
    dplyr::arrange(TIME, CMT) %>%
    dplyr::mutate(dose_number = dplyr::if_else(EVID %in% c(1, 101, 4),
                                               cumsum(EVID %in% c(1, 101, 4)),
                                               NA_integer_)) %>%
    dplyr::ungroup()

  dat <- dat[with(dat, order(ID, resetflag, TIME, CMT,-AMT)), ]
  return(dat)
}
