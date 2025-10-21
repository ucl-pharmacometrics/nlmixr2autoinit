#' Compute overall residual variability from elimination phase
#'
#' Applies `getsigmas` to each individual and dose group after filtering
#' observation records (EVID == 0), and calculates trimmed mean estimates
#' of additive and proportional residual variability.
#'
#' @param df Full pharmacokinetic dataset containing at least the columns:
#'   EVID, ID, TIME, DV, and routeobs.
#' @param nlastpoints Number of terminal points used for elimination phase
#'   regression in each group (passed to getsigmas).
#' @param sigma_trim Trimming proportion used when calculating trimmed means
#'   of residual standard deviations. Default is 0.05.
#'
#' @details
#' The function groups the dataset by subject and dose occasion, applies
#' elimination-phase residual analysis using `getsigmas`, and summarizes the
#' individual residual standard deviations by their trimmed means. This
#' provides population-level estimates of additive and proportional residual
#' unexplained variability (RUV).
#'
#' @return A list containing:
#'   - summary: Named list with trimmed mean values of additive and
#'     proportional residual variability
#'   - full: Data frame with residual estimates for each individual-dose group
#'
#' @seealso \link{getsigmas}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_1CPT
#' dat <- processData(dat)$dat
#' getsigma(dat)
#' }
#'
#' @export

getsigma <- function(df,
                     nlastpoints = 3,
                     sigma_trim = 0.05) {
  # Filter to observation data
  obs_df <- df %>% dplyr::filter(EVID == 0)

  # Apply getsigmas to each group
  full_result <- obs_df %>%
    dplyr::group_by(ID, resetflag, dose_number) %>%
    dplyr::group_modify( ~ getsigmas(.x, nlastpoints = nlastpoints)) %>%
    dplyr::ungroup()

  # Compute trimmed means of residual SDs
  sigma_additive <- tryCatch(
    mean(
      full_result$residual_sd_additive,
      na.rm = TRUE,
      trim = sigma_trim
    ),
    error = function(e)
      NA_real_
  )

  sigma_proportional <- tryCatch(
    mean(
      full_result$residual_sd_proportional,
      na.rm = TRUE,
      trim = sigma_trim
    ),
    error = function(e)
      NA_real_
  )

  summary_result <- list(sigma_additive = sigma_additive,
                         sigma_proportional = sigma_proportional)

  return(list(summary = summary_result,
              full = full_result))
}


#' Estimate individual-level residual error from the elimination phase
#'
#' Performs log-linear regression on the elimination phase of a single individual's
#' or one group's pharmacokinetic concentration–time data to estimate additive and
#' proportional residual standard deviations.
#'
#' @param group_df A data frame for a single group (e.g., one subject or dose),
#'   containing columns: EVID (event ID), DV (observed concentration),
#'   TIME (time after dose), and routeobs (administration route).
#' @param nlastpoints Integer specifying the number of terminal data points used
#'   for regression.
#'
#' @details
#' Residuals are computed from individual-predicted concentrations (IPRED) and
#' observed concentrations (DV) using the following definitions:
#' \deqn{
#'   \sigma_{add} = \sqrt{Var(C_{obs} - C_{pred})}
#' }
#'
#' \deqn{
#'   \sigma_{prop} = \sqrt{Var\left(\frac{C_{obs}}{C_{pred}} - 1\right)}
#' }
#'
#' where \eqn{C_{obs}} is the observed concentration and \eqn{C_{pred}} is the
#' model-predicted concentration obtained by back-transformation of the
#' log-linear regression. The additive residual standard deviation
#' (\eqn{\sigma_{add}}) and proportional residual standard deviation
#' (\eqn{\sigma_{prop}}) are calculated per individual.
#'
#' @return A tibble with the following columns:
#'   - intercept: Intercept of the log-linear regression line
#'   - slope: Estimate of the terminal elimination rate constant
#'   - residual_sd_additive: Standard deviation of additive residuals
#'   - residual_sd_proportional: Standard deviation of proportional residuals
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \dontrun{
#' dat <- Bolus_1CPT
#' dat <- processData(dat)$dat
#' getsigmas(dat[dat$ID == 1 & dat$dose_number == 1 & dat$resetflag == 1 &
#'               dat$EVID == 0, ])
#' }
#'
#' @export

getsigmas <- function(group_df, nlastpoints = 3) {
  if (nrow(group_df) < (nlastpoints + 1)) {
    return(
      tibble::tibble(
        intercept = NA_real_,
        slope = NA_real_,
        residual_sd_additive = NA_real_,
        residual_sd_proportional = NA_real_
      )
    )
  }

  # Identify Tmax
  tmax_index <- which.max(group_df$DV)

  # Determine routeobs
  route <- unique(group_df$routeobs)
  if (length(route) != 1) {
    stop(
      "Multiple routeobs values found in one group (ID = ",
      unique(group_df$ID),
      ", dose_number = ",
      unique(group_df$dose_number),
      "). Check data integrity."
    )
  }

  # Set elimination start index based on routeobs
  if (route %in% c("bolus", "infusion")) {
    elim_start <- tmax_index + 1
  } else if (route == "oral") {
    elim_start <- tmax_index
  } else {
    stop("Unknown routeobs: ", route)
  }

  # Check if enough points after elim_start
  if ((nrow(group_df) - elim_start + 1) < nlastpoints) {
    return(
      tibble::tibble(
        intercept = NA_real_,
        slope = NA_real_,
        residual_sd_additive = NA_real_,
        residual_sd_proportional = NA_real_
      )
    )
  }

  # Select elimination phase data
  elim_df <-
    utils::tail(group_df[elim_start:nrow(group_df),], nlastpoints)

  if (any(elim_df$DV <= 0)) {
    return(
      tibble::tibble(
        intercept = NA_real_,
        slope = NA_real_,
        residual_sd_additive = NA_real_,
        residual_sd_proportional = NA_real_
      )
    )
  }

  # Perform regression
  model_obj <- lm(log(elim_df$DV) ~ elim_df$TIME)
  intercept <- stats::coef(model_obj)[1]
  slope <- stats::coef(model_obj)[2]
  ipred <- exp(intercept + slope * elim_df$TIME)

  residuals_additive <- ipred - elim_df$DV
  residuals_proportional <- (elim_df$DV / ipred) - 1
  sd_additive <- stats::sd(residuals_additive)
  sd_proportional <- stats::sd(residuals_proportional)

  tibble::tibble(
    intercept = intercept,
    slope = slope,
    residual_sd_additive = sd_additive,
    residual_sd_proportional = sd_proportional
  )
}
