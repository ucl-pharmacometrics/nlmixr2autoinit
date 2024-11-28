#' Determine steady state for pharmacokinetic observations
#'
#' Evaluates whether observations in a pharmacokinetic dataset have reached steady state,
#' based on user-defined criteria. It supports two options: using a fixed number of doses or using a
#' fixed number of half-lives to determine the required dose count for steady state.
#'
#' @param df A data frame containing pharmacokinetic data. It should include columns:
#'   - `ID`: Unique identifier for each subject.
#'   - `EVID`: Event ID, where 1 indicates a dose and 0 indicates an observation.
#'   - `SSflag`: Steady-state flag, set to 1 if steady state is manually indicated.
#'   - `TIME`: Observation or dosing time.
#'   - `AMT`: Dose amount.
#'   - `tad`: Time after dose.
#' @param ss_option Integer, either 1, 2, or 3, indicating the method to determine steady state:
#'   - `1`: Steady state is determined based on the smaller value between the fixed number of doses (`no.doses`)
#'          and the estimated number of doses based on half-life.
#'   - `2`: Steady state is determined solely based on a fixed number of doses (`no.doses`).
#'   - `3`: Steady state is determined by calculating the number of doses as the dosing duration
#'          (a fixed number of half-lives, `no.half_life`, multiplied by `half_life`) divided by the dose interval.
#' @param no.doses Integer. Number of doses required to reach steady state when `ss_option = 1`. Default is 5.
#' @param no.half_life Integer. Number of half-lives required to reach steady state when `ss_option = 2`. Default is 3.
#' @param half_life Numeric. The half-life of the drug, required if `ss_option = 2`.
#'
#' @return A data frame similar to `df` with two additional columns:
#'   - `SteadyState`: Logical, indicating whether each observation has reached steady state.
#'   - `recent_ii`: Numeric, indicating the recent dosing interval for each observation that has achieved steady state.
#'
#' @details
#' The function evaluates steady state for each subject (`ID`) individually. It provides two options:
#'   - If `ss_option = 1`, steady state is determined based on the smaller value between the number of fixed doses (`no.doses`)
#'     and the estimated number of doses required based on half-life.
#'   - If `ss_option = 2`, steady state is determined solely based on a fixed number of doses (`no.doses`).
#'   - If `ss_option = 3`, steady state is determined by calculating the number of doses required as the dosing duration
#'     (a fixed number of half-lives, `no.half_life`, multiplied by `half_life`) divided by the dose interval.
#'
#' For each observation, the following criteria are used to assess steady state:
#'   1. If `SSflag = 1`, the observation is immediately marked as steady state.
#'   2. If `SSflag = 0`, the function checks that:
#'      - The required number of doses has been administered before the observation.
#'      - The dosing intervals are consistent, within 50% of the median dose interval.
#'      - The dose amounts are consistent, within 25% of the median dose.
#'      - The observation falls within the last dosing interval, allowing a 25% variance.
#'
#' @examples
#' # Example usage
#' df <- data.frame(
#'  ID = rep(1, 14),
#'  EVID = c(1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0),
#'  SSflag = rep(0, 14),  # Initialize SSflag as 0
#'  TIME = 1:14,
#'  AMT = c(100, 0, 0, 100, 0, 0, 100, 100, 100, 100, 100, 0, 0, 0),
#'  tad = runif(14)
#'  )
#'
#' is_ss(df,ss_option=1)
#' is_ss(df,ss_option=2)
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup
#' @importFrom crayon black
#' @import purrr
#' @export

is_ss <- function(df,
                  ss_option = 1,
                  no.doses = 5,
                  no.half_life = 5,
                  half_life = NA) {
  if (missing(df)) {
    stop("Error, no dataset provided")
  }

  time_to_ss <- NA

  if (!is.na(half_life)) {
    # Calculate number of doses required to reach steady state
    time_to_ss <- no.half_life * half_life
  }

  # Consider SS=1 column

  df$SteadyState <- FALSE
  df$recent_ii <- 0
  # df$SS.method <- NA


  # Calculate `previous_doses` and `previous_amts` for each observation time
  df <- df %>%
    group_by(ID) %>%
    mutate(
      # Filter `TIME` and `AMT` for dosing events where `EVID` indicates a dose (e.g., EVID in c(4, 101, 1))
      previous_doses = purrr::map(TIME, ~ TIME[EVID %in% c(4, 101, 1) & TIME <= .]),
      previous_amts = purrr::map(TIME, ~ AMT[EVID %in% c(4, 101, 1) & TIME <= .])
    )

  # Step 2: Calculate required number of doses based on fixed number of half-lives
  df <- df %>%
    mutate(
      # Define doses_required_1 as a fixed number of doses based on `no.doses`
      doses_required_1 = no.doses
    ) %>%
    mutate(
      #  Calculate last_two_doses_interval for each observation
      last_two_doses_interval = purrr::map2_dbl(previous_doses, TIME, ~ {
        last_two_doses <- tail(.x, 2)
        if (length(last_two_doses) == 2) {
          # Calculate and return the interval between the last two doses
          diff(last_two_doses)
        } else {
          NA_real_ # Return NA if there are not enough doses
        }
      })
    ) %>% # Calculate doses_required_2 (half-lives)
    mutate(
      doses_required_2 = ifelse(
        !is.na(last_two_doses_interval) & !is.na(time_to_ss),
        pmax(ceiling(time_to_ss / last_two_doses_interval),3),  # Calculate doses needed to reach steady state
        NA_real_  # Assign NA if conditions are not met
      )
    )

  # Calculate Final doses_required Based on ss_option
  df <- df %>%
    mutate(
      doses_required = case_when(
        ss_option == 1 ~ pmin(doses_required_1, doses_required_2, na.rm = TRUE), # Take minimum of two options
        ss_option == 2 ~ doses_required_1, # Use fixed doses if `ss_option` is 2
        ss_option == 3 ~ ifelse(!is.na(doses_required_2), doses_required_2, doses_required_1) # Based on half-lives if available
      ),
      doses_required = pmax(doses_required, 3)  # Ensure minimum of 3 doses
    )


  # df <- df %>%
  #   mutate(
  #     doses_to_check = purrr::map2(previous_doses, doses_required, ~ tail(.x, .y)),
  #     amts_to_check = purrr::map2(previous_amts, doses_required, ~ tail(.x, .y)),
  #     intervals = purrr::map(doses_to_check, diff),
  #     dose_interval = purrr::map_dbl(intervals, median, na.rm = TRUE),
  #     is_continuous = purrr::map2_lgl(intervals, dose_interval, ~ all(abs(.x - .y) <= .y * 0.5)),
  #     is_same_dose = purrr::map_lgl(amts_to_check, ~ all(abs(.x - median(.x, na.rm = TRUE)) <= median(.x, na.rm = TRUE) * 1.25))
  #   )


  df <- df %>%
    mutate(
      dose_count_before_obs =  purrr::map_int(previous_doses, length),
      # Extract the last `doses_required` doses and amounts
      doses_to_check = purrr::map2(previous_doses, doses_required, ~ tail(.x, .y)),
      amts_to_check = purrr::map2(previous_amts, doses_required, ~ tail(.x, .y)),

      # Calculate intervals between doses and the median interval
      intervals = purrr::map(doses_to_check, diff),
      dose_interval = purrr::map_dbl(last_two_doses_interval, median, na.rm = TRUE),

      # Check if all intervals are within 1.5 times the median interval
      is_continuous = purrr::map2_lgl(intervals, dose_interval, ~ all(abs(.x - .y) <= .y * 0.25)),

      # Check if all dose amounts are within 1.25 times the median dose amount
      is_same_dose = purrr::map_lgl(amts_to_check, ~ all(abs(.x - median(.x, na.rm = TRUE)) <= median(.x, na.rm = TRUE) * 1.25)),

      # Check if each observation occurs within dose_interval * 1.25
      is_within_last_dose_interval = tad < dose_interval * 1.25
    )

  # Determine steady state and recent_ii
  df <- df %>%
    mutate(
      SteadyState = ifelse(
        (SSflag == 1 & EVID == 0) | ( (dose_count_before_obs >= doses_required) &  is_continuous & is_same_dose & is_within_last_dose_interval),
        TRUE, FALSE
      ),
      SS.method = case_when(
        SteadyState & doses_required == doses_required_1 ~ "A fixed number of doses",
        SteadyState & doses_required == doses_required_2 ~ "A fixed number of half-lives",
        TRUE ~ NA_character_
      ),
      recent_ii = ifelse(SteadyState, dose_interval, NA)
    ) %>%
    ungroup()


  return(df)
}
