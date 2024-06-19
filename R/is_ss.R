#' Identify steady-state observations
#'
#' Checks if a sufficient number of doses have been administered before each observation point to determine if the observation point has reached a steady state, and calculates the time required to reach steady state based on the provided half-life and dose interval, and marks observations as steady state if they meet the criteria.
#' @param df A data frame containing the pharmacokinetic data. The data frame must include columns for ID, TIME, EVID, AMT, and tad (time after dose).
#' @param half_life The half-life of the drug.
#' @param dose_interval The dosing interval.
#' @return A data frame with an additional column `SteadyState` indicating whether each observation is at steady state.
#' @importFrom dplyr %>% mutate if_else arrange group_by ungroup
#' @examples
#' df <- data.frame(ID = rep(1, 10), EVID = c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1), TIME = 1:10, AMT = c(100, 0, 0, 100, 0, 0, 100, 0, 0, 100), tad = runif(10))
#' is_ss(df, half_life = 4, dose_interval = 8)
#' @export

is_ss <- function(df,
                  half_life,
                  dose_interval) {
  if (missing(df)) {
    stop("Error, no dataset provided")
  }

  if (missing(half_life)) {
    stop("Error, no half life provided")
  }

  if (missing(dose_interval)) {
    stop("Error, no dose_interval provided")
  }

  nhalf <- 5 # commonly 5 half life needed to reach steady-state

  for (nhalfloop in seq(nhalf, 3, -1)) {
    # Calculate time to reach steady state, 5 half life was set as default
    time_to_ss <- nhalfloop * half_life
    # Calculate number of doses required to reach steady state
    doses_required <- ceiling(time_to_ss / dose_interval)
    df$SteadyState <- FALSE

    for (id in unique(df$ID)) {
      id_df <- df[df$ID == id, ]
      id_obs_df <- df[df$ID == id & df$EVID == 0, ]

      obs_times <- id_obs_df$TIME
      dose_times <- id_df[id_df$EVID %in% c(4, 101, 1), ]$TIME
      dose_amts <- id_df[id_df$EVID %in% c(4, 101, 1), ]$AMT

      for (obsi in obs_times) {
        # Find the doses before the current observation time
        previous_doses <- dose_times[dose_times <= obsi]
        previous_amts <- dose_amts[dose_times <= obsi]

        # Check if there are enough doses before observation
        if (length(previous_doses) >= doses_required) {
          # Extract the relevant dose times for checking
          doses_to_check <- tail(previous_doses, doses_required)
          amts_to_check <- tail(previous_amts, doses_required)

          # Check that there are no missed doses or dose interruptions, but allow for variations in the dose interval.
          # Ensure that tad is less than the dose interval, and confirm that the observation falls within a dose interval.
          # Verify that the amount of each dose is the same.

          if (all(diff(doses_to_check) < dose_interval * 1.5) &
              all(amts_to_check == amts_to_check[1]) &
              df[df$ID == id &
                 df$TIME == obsi &
                 df$EVID == 0,]$tad < dose_interval) {
            # Mark the corresponding observation as steady state
            df[df$ID == id &
                 df$TIME == obsi &
                 df$EVID == 0,]$SteadyState <- TRUE
          }
        }
      }

    }

    # Minimum points at the steady state needed for statistics calculation
    if (nrow(df[df$EVID == 0 & df$SteadyState == T, ]) > 2) {
      if (nhalfloop != nhalf) {
        print(
          paste0(
            "Warning: required number of half life was decreased to ",
            nhalfloop,
            ", due to failure to find enough points based on current half life, ",
            half_life
          )
        )
      }
      break
    }

    if (nhalfloop == 3) {
      print(paste0("Warning: No steady-state concentration point found"))
    }

  }
  return(df)
}
