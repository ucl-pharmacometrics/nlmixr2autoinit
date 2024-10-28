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
#' @param ss_option Integer, either 1 or 2, indicating the method to determine steady state:
#'   - `1`: Use a fixed number of doses specified by `no.doses`.
#'   - `2`: Use a fixed number of half-lives specified by `no.half_life` and `half_life`.
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
#'   - If `ss_option = 1`, steady state is determined by a fixed number of doses (`no.doses`).
#'   - If `ss_option = 2`, steady state is determined based on a fixed number of half-lives (`no.half_life`)
#'     multiplied by `half_life` to estimate the required dosing duration.
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
#' @export




is_ss <- function(df,
                  ss_option = 1,
                  no.doses = 5,
                  no.half_life = 3,
                  half_life = NA) {
  if (missing(df)) {
    stop("Error, no dataset provided")
  }

  if (is.na(half_life)) {
    ss_option = 1
    message(black("No half life available, change to ss_option to 1"))
  }

  if (ss_option == 1) {
    doses_required <- no.doses
    message(
      black(
        "A fixed number of doses ",
        "(",
        no.doses,
        ") is used as the criterion to determine if steady state has been achieved."
      )
    )
  }

  if (ss_option == 2) {
    # Calculate number of doses required to reach steady state
    time_to_ss <- no.half_life * half_life
    doses_required <-
      ceiling(time_to_ss / dose_interval) #  rounded up to the next whole number.

    message(
      black(
        "A fixed number of half-lives",
        "(",
        no.half_life,
        ") is used to assess whether the dose count meets the criterion for achieving steady state"
      )
    )
  }

  # Consider SS=1 column

  df$SteadyState <- FALSE
  df$recent_ii <- 0
  # df[df$SSflag == 1 & df$EVID == 0, ]$SteadyState = TRUE

  for (id in unique(df$ID)) {
    id_df <- df[df$ID == id,]
    id_obs_df <- df[df$ID == id & df$EVID == 0,]

    obs_times <- id_obs_df$TIME
    dose_times <- id_df[id_df$EVID %in% c(4, 101, 1),]$TIME
    dose_ii <- id_df[id_df$EVID %in% c(4, 101, 1),]$II
    dose_amts <- id_df[id_df$EVID %in% c(4, 101, 1),]$AMT

    # Determine steady state by evaluating each sampling point individually.
    for (obsi in obs_times) {
      # if SSflag=1
      if (df[df$ID == id &  df$TIME == obsi, ]$SSflag == 1) {
        df[df$ID == id &  df$TIME == obsi, ]$SteadyState = T
        df[df$ID == id &
             df$TIME == obsi, ]$recent_ii =  tail(dose_ii[dose_ii <= obsi] , 1) # last one
      }

      if (df[df$ID == id &  df$TIME == obsi, ]$SSflag == 0) {
        # Find the doses before the current observation time
        previous_doses <- dose_times[dose_times <= obsi]
        previous_amts <- dose_amts[dose_times <= obsi]

        # Check if there are enough doses before observation
        if (length(previous_doses) >= doses_required) {
          # Extract the last few doses, with the number of doses specified by doses_required.
          doses_to_check <- tail(previous_doses, doses_required)
          amts_to_check <- tail(previous_amts, doses_required)

          # Calculate intervals between consecutive dosing times
          intervals <- diff(doses_to_check)

          # # Identify the most common interval as the standard interval
          # dose_interval <- as.numeric(names(sort(table(intervals), decreasing = TRUE)[1]))

          # Calculate the median interval to reduce the impact of outliers
          dose_interval <- median(intervals)

          # Check if all intervals are within 1.5 times the median interval
          is_continuous <-
            all(abs(intervals - dose_interval) <= dose_interval * 0.5)

          # Check if all dose amounts are within 1.25 times the median amount
          is_same_dose <-
            all(abs(amts_to_check - median(amts_to_check)) <= median(amts_to_check) * 1.25) # 20% amt difference

          # Check  if the observation falls within the last dose interval.
          is_within_last_dose_interval <-
            df[df$ID == id &
                 df$TIME == obsi &
                 df$EVID == 0, ]$tad < dose_interval * 1.25# allow 1.25 difference

          if (is_continuous &
              is_same_dose & is_within_last_dose_interval) {
            # Mark the corresponding observation as steady state
            df[df$ID == id &
                 df$TIME == obsi &  df$EVID == 0, ]$SteadyState <- TRUE
            df[df$ID == id &
                 df$TIME == obsi &  df$EVID == 0, ]$recent_ii <- median(intervals)
          }

        } # end previous_doses if

      } # end SSflag if

    } # end obs loop

  } # end id loop

  return(df)
}
