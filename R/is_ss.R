#' Internal control builder for steady-state evaluation
#'
#' Constructs a list of control parameters used by is_ss() to determine
#' pharmacokinetic steady state.
#'
#' @param ss_method Character string specifying the method used to determine
#'   steady state. One of:
#'   \itemize{
#'     \item "combined" (default): uses the smaller of the dose-based estimate
#'           (no.doses) and the half-life-based estimate (no.half_lives)
#'     \item "fixed_doses": considers steady state reached after no.doses
#'           administrations
#'     \item "half_life_based": uses the drug half-life and dosing interval to
#'           estimate the required number of doses
#'   }
#'
#' @param no.doses Integer indicating the number of doses assumed necessary to
#'   reach steady state when using the "fixed_doses" method or as part of the
#'   "combined" method. Default is 5.
#'
#' @param no.half_lives Integer indicating the number of half-lives required to
#'   reach steady state when using the "half_life_based" or "combined" method.
#'   Default is 5.
#'
#' @param allowed_interval_variation Numeric value specifying the acceptable
#'   fractional variation in dose interval. For example, 0.25 allows plus or
#'   minus 25 percent variation. Default is 0.25.
#'
#' @param allowed_dose_variation Numeric value specifying the acceptable
#'   fractional variation in dose amount. For example, 0.20 allows plus or minus
#'   20 percent variation. Default is 0.20.
#'
#' @param min_doses_required Integer specifying the minimum number of doses that
#'   must be administered regardless of method. Default is 3.
#'
#' @param tad_rounding Logical value. If TRUE (default), rounding is applied
#'   when comparing time after dose (tad) to dosing intervals to allow small
#'   numerical deviations.
#'
#' @return A named list containing the steady-state control parameters,
#'   typically passed as the ssctrl argument to is_ss().
#'
#' @seealso \link{is_ss}
#'
#' @examples
#' ss_control()
#' ss_control(ss_method = "fixed_doses", no.doses = 4)
#'
#' @export

ss_control <-
  function(ss_method = c("combined", "fixed_doses", "half_life_based"),
           no.doses = 5,
           no.half_lives = 5,
           allowed_interval_variation = 0.25,
           allowed_dose_variation = 0.20,
           min_doses_required = 3,
           tad_rounding = TRUE) {
    # Safe matching for ss_method
    ss_method <- tryCatch(
      match.arg(
        ss_method,
        choices = c("combined", "fixed_doses", "half_life_based")
      ),
      error = function(e) {
        stop(
          sprintf(
            "Invalid value for `%s`: '%s'. Must be one of: %s.",
            "ss_method",
            as.character(ss_method),
            paste(shQuote(
              c("combined", "fixed_doses", "half_life_based")
            ), collapse = ", ")
          ),
          call. = FALSE
        )
      }
    )
    list(
      ss_method = ss_method,
      no.doses = no.doses,
      no.half_lives = no.half_lives,
      allowed_interval_variation = allowed_interval_variation,
      allowed_dose_variation = allowed_dose_variation,
      min_doses_required = min_doses_required,
      tad_rounding = tad_rounding
    )
  }

#' Determine steady state for pharmacokinetic observations
#'
#' Evaluates whether pharmacokinetic observations have reached steady state
#' based on user-defined control settings. The classification can be based on
#' a fixed number of doses, the number of half-lives relative to the dosing
#' interval, or a combination of both criteria.
#'
#' @param df A data frame containing pharmacokinetic data. It should include
#'   columns for ID, EVID, SSflag, TIME, AMT, and tad.
#'
#' @param ssctrl A control list consistent with the structure returned by
#'   ss_control(). It specifies the method and thresholds for steady-state
#'   evaluation.
#'
#' @param half_life Numeric value representing the drug half-life. Required
#'   when the method in ss_control() is based on half-life or uses a combined
#'   approach.
#'
#' @details
#' The function determines steady state by examining each observation in
#' relation to prior dosing history. The required number of doses is calculated
#' based on the specified method in ss_control(). Observation times are
#' evaluated to confirm that dose interval and dose amount variability fall
#' within acceptable limits and that the time after dose is within the most
#' recent dosing interval. Observations manually marked as steady state using
#' SSflag are also recognized as steady state.
#'
#' @return A data frame with added columns indicating steady-state status,
#'   the dosing interval for steady-state observations, and the method used
#'   to classify steady state.
#'
#' @examples
#' dat <- pheno_sd
#' dat <- processData(dat)$dat
#' out <- is_ss(df = dat)
#' out[out$SteadyState == TRUE & !is.na(out$SteadyState),
#'     c("ID", "TIME", "DV", "EVID", "SteadyState")]
#'
#' @seealso ss_control()
#' @export


is_ss <- function(df,
                  ssctrl = ss_control(),
                  half_life = NA) {
  `%>%` <- magrittr::`%>%`

  # ---- Extract control parameters from ssctrl ----
  ss_method <- ssctrl$ss_method
  no.doses <- ssctrl$no.doses
  no.half_lives <- ssctrl$no.half_lives
  allowed_interval_variation <- ssctrl$allowed_interval_variation
  allowed_dose_variation <- ssctrl$allowed_dose_variation
  min_doses_required <- ssctrl$min_doses_required
  tad_rounding <- ssctrl$tad_rounding

  # ---- Defensive checks ----
  if (missing(df) ||
      !is.data.frame(df))
    stop("Input `df` must be a data.frame.")

  required_cols <- c("TIME", "EVID", "AMT", "ID", "tad", "SSflag")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  valid_methods <- c("combined", "fixed_doses", "half_life_based")
  ss_method <- tryCatch(
    match.arg(ss_method, choices = valid_methods),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "ss_method",
          as.character(ss_method),
          paste(shQuote(valid_methods), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # Compute time to steady state
  time_to_ss <-
    if (!is.na(half_life))
      no.half_lives * half_life
  else
    NA_real_

  # Calculate `dose_time_hist` and `dose_amt_hist` for each observation time
  df <- df %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(
      SteadyState = FALSE,
      recent_ii = 0,
      dose_time_hist = purrr::map(TIME, ~ TIME[EVID %in% c(4, 101, 1) &
                                                 TIME <= .]),
      dose_amt_hist = purrr::map(TIME, ~ AMT[EVID %in% c(4, 101, 1) &
                                               TIME <= .])
    )

  df <- df %>%
    dplyr::mutate(
      # Step 1: Fixed dose calculation (doses_req_fixed)
      doses_req_fixed = no.doses,

      # Step 2: Calculate dosing interval between last two doses
      last_two_doses_interval = purrr::map2_dbl(.x = dose_time_hist, .y = TIME,
                                                .f = ~ {
                                                  last_two_doses <- tail(.x, 2)
                                                  if (length(last_two_doses) == 2) {
                                                    diff(last_two_doses)  # Calculate time difference
                                                  } else {
                                                    NA_real_             # Return NA if <2 doses exist
                                                  }
                                                }),

      # Step 3: Half-life-based dose calculation (doses_req_hl)
      doses_req_hl = ifelse(
        !is.na(time_to_ss) & !is.na(last_two_doses_interval),
        pmax(ceiling(time_to_ss / last_two_doses_interval), 3),
        NA_real_
      )
    )

  # Calculate Final doses_required Based on ss_method
  df <- df %>%
    dplyr::mutate(
      doses_required = dplyr::case_when(
        ss_method == "combined"  ~ pmin(doses_req_fixed,
                                        doses_req_hl, na.rm = TRUE),
        ss_method == "fixed_doses"    ~ doses_req_fixed,
        ss_method == "half_life_based" ~ ifelse(!is.na(doses_req_hl),
                                                doses_req_hl, doses_req_fixed)
      ),
      doses_required = pmax(doses_required, min_doses_required)
    )

  df <- df %>%
    dplyr::mutate(
      dose_count_before_obs =  purrr::map_int(dose_time_hist, length),
      # Extract the last `doses_required` doses and amounts
      doses_to_check = purrr::map2(dose_time_hist,
                                   doses_required, ~ as.numeric(tail(.x, .y))),
      amts_to_check = purrr::map2(dose_amt_hist,
                                  doses_required, ~ as.numeric(tail(.x, .y))),

      # Calculate intervals between doses and the median interval
      intervals = purrr::map(doses_to_check, diff),
      dose_interval = purrr::map_dbl(last_two_doses_interval, median, na.rm = TRUE),

      # Check if all intervals are within ±25% (default) of the last interval
      is_continuous = purrr::map_lgl(intervals, ~ {
        if (length(.x) < 2)
          return(FALSE)
        last_val <- tail(.x, 1)
        all(abs(.x - last_val) <= last_val * (1 + allowed_interval_variation))
      }),

      # Check if all dose amounts are within ±20% (default) of the last dose amount
      is_same_dose = purrr::map_lgl(amts_to_check, ~ {
        if (length(.x) < 2)
          return(FALSE)
        last_val <- tail(.x, 1)
        all(abs(.x - last_val) <= last_val * (1 + allowed_dose_variation))
      }),

      is_within_last_dose_interval = if (tad_rounding) {
        round(tad) <= round(dose_interval)
      } else {
        tad <= dose_interval
      }
    )

  # Determine steady state and recent_ii
  df <- df %>%
    dplyr::mutate(
      SteadyState = ifelse(
        (SSflag == 1 & EVID == 0) |
          (
            (dose_count_before_obs >= doses_required) &
              is_continuous &
              is_same_dose & is_within_last_dose_interval
          ),
        TRUE,
        FALSE
      ),
      SS.method = dplyr::case_when(
        SteadyState &
          (SSflag == 1 &
             EVID == 0) ~ "Protocol-defined steady state",
        SteadyState &
          doses_required == doses_req_fixed ~
          paste0("Pre-specified fixed doses (ndoses = ", doses_required, ")"),
        SteadyState &
          doses_required == doses_req_hl ~
          paste0("Half-life-based steady state (ndoses = ", doses_required, ")"),
        TRUE ~ NA_character_
      ),

      recent_ii = dplyr::case_when(
        SteadyState & SSflag == 1 & EVID == 0 ~ iiobs,
        SteadyState ~ dose_interval,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::ungroup()

  # Remove internal list-columns to prevent large output size
  df <- df %>%
    dplyr::select(-dose_time_hist, -dose_amt_hist)

  return(df)
}

