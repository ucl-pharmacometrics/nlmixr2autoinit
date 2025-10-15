#' Compute overall residual variability (sigma) from grouped elimination phase analysis
#'
#' This function applies `getsigmas()` to each individual-dose group after filtering
#' to retain only observation records (EVID == 0), and returns trimmed mean estimates
#' of residual additive and proportional variability.
#'
#' @inheritParams getsigmas
#' @param df The full PK dataset with required columns: `EVID`, `ID`, `TIME`, `DV`, `routeobs`, etc.
#' @param sigma_trim Proportion of trimmed mean to calculate (default = 0.05).
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Named list of trimmed mean residual SDs}
#'   \item{full}{Tibble with full results for each group}
#' }
#'
#' @examples
#' \dontrun{
#' dat<- Bolus_1CPT
#' dat<-processData(dat)$dat
#' getsigma(dat)
#'}
#' @export
#'
getsigma <- function(df,
                     nlastpoints = 3,
                     sigma_trim=0.05) {

  # Filter to observation data
  obs_df <- df %>% dplyr::filter(EVID == 0)

  # Apply getsigmas to each group
  full_result <- obs_df %>%
    dplyr::group_by(ID, resetflag, dose_number) %>%
    dplyr::group_modify(~ getsigmas(.x, nlastpoints = nlastpoints)) %>%
    dplyr::ungroup()

  # Compute trimmed means of residual SDs
  sigma_additive <- tryCatch(
    mean(full_result$residual_sd_additive, na.rm = TRUE, trim = sigma_trim),
    error = function(e) NA_real_
  )

  sigma_proportional <- tryCatch(
    mean(full_result$residual_sd_proportional, na.rm = TRUE, trim = sigma_trim),
    error = function(e) NA_real_
  )

  summary_result <- list(
    sigma_additive = sigma_additive,
    sigma_proportional = sigma_proportional
  )

  return(list(
    summary = summary_result,
    full = full_result
  ))
}


#' Estimate individual-level residual error from elimination phase
#'
#' Performs log-linear regression on the elimination phase of a
#' single individual's (or one group’s) pharmacokinetic data to estimate additive
#' and proportional residual standard deviations.
#'
#' @param group_df A data frame for a single group (e.g., one subject/dose),
#'   containing columns `EVID`, `DV`, `TIME`, and `routeobs`.
#' @param nlastpoints Integer. Number of terminal points to include in regression.
#'
#' @details
#' The elimination phase is selected differently based on the `routeobs` value:
#' - For `bolus` and `infusion`, it starts *after* Tmax.
#' - For `oral`, it starts *at* Tmax.
#'
#' The function uses the last `nlastpoints` data points in the elimination phase.
#' Residuals are calculated as:
#' - Additive: \eqn{IPRED - DV}
#' - Proportional: \eqn{DV / IPRED - 1}
#'
#' Returns NA values if the data is insufficient or inconsistent.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{intercept}{Intercept of the regression line}
#'   \item{slope}{Slope of the regression line}
#'   \item{residual_sd_additive}{Standard deviation of additive residuals}
#'   \item{residual_sd_proportional}{Standard deviation of proportional residuals}
#' }
#'
#' @examples
#' \dontrun{
#' dat<- Bolus_1CPT
#' dat<-processData(dat)$dat
#' getsigmas(dat[dat$ID==1 & dat$dose_number==1 & dat$resetflag==1&dat$EVID==0,])
#' }
#'
#'
#' @export
#'

getsigmas <- function(group_df, nlastpoints = 3) {

  if (nrow(group_df) < (nlastpoints + 1)) {
    return(tibble::tibble(
      intercept = NA_real_,
      slope = NA_real_,
      residual_sd_additive = NA_real_,
      residual_sd_proportional = NA_real_
    ))
  }

  # Identify Tmax
  tmax_index <- which.max(group_df$DV)

  # Determine routeobs
  route <- unique(group_df$routeobs)
  if (length(route) != 1) {
    stop("Multiple routeobs values found in one group (ID = ",
         unique(group_df$ID),
         ", dose_number = ", unique(group_df$dose_number),
         "). Check data integrity.")
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
    return(tibble::tibble(
      intercept = NA_real_,
      slope = NA_real_,
      residual_sd_additive = NA_real_,
      residual_sd_proportional = NA_real_
    ))
  }

  # Select elimination phase data
  elim_df <- utils::tail(group_df[elim_start:nrow(group_df), ], nlastpoints)

  if (any(elim_df$DV <= 0)) {
    return(tibble::tibble(
      intercept = NA_real_,
      slope = NA_real_,
      residual_sd_additive = NA_real_,
      residual_sd_proportional = NA_real_
    ))
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




