#' Internal control builder for steady-state evaluation
#'
#' Constructs a list of control parameters used by \code{\link{is_ss}()} to
#' determine pharmacokinetic steady state.
#'
#' @param ss_method Character. Method used to determine steady state.
#'   One of:
#'   \itemize{
#'     \item \code{"combined"} (default): Uses the smaller of the dose-based
#'       estimate (\code{no.doses}) and the half-life-based estimate
#'       (\code{no.half_lives}).
#'     \item \code{"fixed_doses"}: Considers steady state reached after
#'       \code{no.doses} administrations.
#'     \item \code{"half_life_based"}: Uses the drug half-life and dosing
#'       interval to estimate the required number of doses.
#'   }
#' @param no.doses Integer. Number of doses assumed to reach steady state
#'   when \code{ss_method = "fixed_doses"} or part of the \code{"combined"} method.
#'   Default is 5.
#' @param no.half_lives Integer. Number of half-lives required to reach
#'   steady state when using \code{"half_life_based"} or \code{"combined"} methods.
#'   Default is 5.
#' @param allowed_interval_variation Numeric. Acceptable fractional variation in
#'   dose interval (e.g., 0.25 allows ±25% variation). Default is 0.25.
#' @param allowed_dose_variation Numeric. Acceptable fractional variation in
#'   dose amount (e.g., 0.20 allows ±20% variation). Default is 0.20.
#' @param min_doses_required Integer. Minimum number of doses that must be
#'   administered regardless of method. Default is 3.
#' @param tad_rounding Logical. If \code{TRUE} (default), applies rounding when
#'   comparing observation times (\code{tad}) to dosing intervals, allowing
#'   small numerical deviations (e.g., 24.3 ≈ 24).
#'
#' @return A named list containing the steady-state control parameters.
#'   The returned object is typically passed as the \code{ssctrl} argument
#'   to \code{\link{is_ss}()}.
#'
#' @seealso \code{\link{is_ss}}
#' @export
#' @examples
#' ss_control()
#' ss_control(ss_method = "fixed_doses", no.doses = 4)
#'
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
#' Evaluates whether observations in a pharmacokinetic dataset have reached steady state,
#' based on user-defined criteria. It supports three methods: using a fixed number of doses,
#' using a fixed number of half-lives, or combining both criteria.
#'
#' @param df A data frame containing pharmacokinetic data. It should include columns:
#'   - `ID`: Unique identifier for each subject.
#'   - `EVID`: Event ID, where 1 indicates a dose and 0 indicates an observation.
#'   - `SSflag`: Steady-state flag, set to 1 if steady state is manually indicated.
#'   - `TIME`: Observation or dosing time.
#'   - `AMT`: Dose amount.
#'   - `tad`: Time after the last dose.
#' @param ssctrl A list created by \code{\link{ss_control}()}, specifying control settings
#'   for steady-state evaluation. See \code{\link{ss_control}} for details on available options.
#' @param half_life Numeric. The drug's half-life. Required if `ss_method` in `ssctrl` is
#'   `"combined"` or `"half_life_based"`.
#'
#' @return A data frame similar to `df` with additional columns:
#'   - `SteadyState`: Logical, indicating whether each observation has reached steady state.
#'   - `recent_ii`: Numeric, the most recent dosing interval for observations marked as steady state.
#'   - `SS.method`: Character, description of the method used to identify steady state.
#'
#' @examples
#' # Example usage
#' dat <- pheno_sd
#' dat <- processData(dat)$dat
#' out<-is_ss(df = dat)
#' out[out$SteadyState == TRUE & !is.na(out$SteadyState),
#'            c("ID", "TIME", "DV", "EVID", "SteadyState")]

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

