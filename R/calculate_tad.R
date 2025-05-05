#' Calculate time after dose (TAD) for pharmacokinetic aata
#'
#' Computes time after dose (TAD) for pharmacokinetic observations. For observation records,
#' propagates dosing information from the most recent dose administration.
#'
#' @param dat A pharmacokinetic dataset containing:
#'   - **Required columns**: `ID`, `TIME`, `EVID`, `AMT`
#'
#' @return A modified dataframe with added columns:
#'   - `tad`: Time after dose (calculated as TIME minus the time of the most recent prior dose) for observations (`EVID == 0`); `NA` for dosing records.
#'   - `iiobs`: The inter-dose interval inherited from the most recent dose, applied to observation rows.
#'   - `rateobs`: The infusion rate inherited from the most recent dose, applied to observation rows.
#'   - `routeobs` *(optional)*: The administration route inherited from the most recent dose, applied to observation rows. This column is only included if the `route` information is available.
#'   - `dose_number`: Sequential dose number, automatically generated if not already present.
#'
#' @details
#' ### Key Operations:
#' 1. **Propagation of Dose Information**:
#'    - For **observation rows** (`EVID == 0`):
#'      - `iiobs`, `rateobs`, and `routeobs` inherit values from the most recent prior dose
#'      - `tad` is calculated as the time since the last dose
#'    - For **dose rows** (`EVID %in% c(1, 101, 4)`):
#'      - Retains original values for `II`, `RATE`, and `route`
#'      - Serves as the reference point for subsequent `tad` calculations
#'
#' 2. **Dose Number Assignment**:
#'    - Automatically generates a sequential `dose_number` using `mark_dose_number()` if not already present
#'    - Preserves existing `dose_number` values if available
#'
#'
#' @export
#'
#' @examples
#' # Bolus dosing
#' calculate_tad(Bolus_1CPT)
#' calculate_tad(Infusion_1CPT)
#' calculate_tad(Oral_1CPT)
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

  dat <- dat %>% mutate(tad = NA_real_)
  # Identify concentration rows
  conc_rows <- dat$EVID == 0
  # Filter the dose lines
  dat <- dat %>%
    dplyr::arrange(ID, resetflag, TIME) %>%
    dplyr::group_by(ID, resetflag) %>%
    dplyr::mutate(
      last_dose_time = if_else(EVID %in% c(1, 101, 4), TIME, NA_real_),
      last_dose_number = if_else(EVID %in% c(1, 101, 4), dose_number, NA_integer_),
      last_dose = if_else(EVID %in% c(1, 101, 4), AMT, NA_integer_),
      last_ii = if_else(EVID %in% c(1, 101, 4), II, NA_real_),
      last_rate = if_else(EVID %in% c(1, 101, 4), RATE, NA_real_),
      last_route = if_else(EVID %in% c(1, 101, 4), route, NA_character_)
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
      tad = if_else(conc_rows, TIME - last_dose_time, NA_real_),
      # if conc row, fill by last_dose_number, else, dose_number
      dose_number = if_else(conc_rows, last_dose_number, dose_number),
      dose = if_else(conc_rows, last_dose, AMT),
      iiobs = if_else(conc_rows, last_ii, II),
      rateobs = if_else(conc_rows, last_rate, RATE),
      routeobs = if_else(conc_rows, last_route, route)
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
#' Marks each dosing event in the dataset with a sequential dose number,
#' based on rows identified by `EVID`. Dose numbering is performed within
#' groups defined by `ID` and `resetflag`, ensuring numbering restarts
#' appropriately for each group.
#'
#' @param dat PK dataset containing columns:
#'   - **Required**: `ID`, `TIME`, `EVID`, `AMT`
#'   - **Optional**: `resetflag` (default=1 if missing), `CMT`
#'
#' @return A modified dataframe with added `dose_number` column
#'
#' @examples
#' mark_dose_number(Bolus_1CPT)
#' mark_dose_number(Infusion_1CPT)
#' mark_dose_number(Oral_1CPT)
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
    dplyr::mutate(dose_number = if_else(EVID %in% c(1, 101, 4),
                                        cumsum(EVID %in% c(1, 101, 4)),
                                        NA_integer_)) %>%
    dplyr::ungroup()

  dat <- dat[with(dat, order(ID, resetflag, TIME, CMT,-AMT)), ]
  return(dat)
}
