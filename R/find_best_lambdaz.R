#' Find the best terminal elimination rate constant (Lambda Z)
#'
#' Identifies the optimal terminal phase for lambda-z estimation using a systematic
#' log-linear regression approach with adjusted R-squared optimization criteria.
#'
#' @param time Numeric vector of observation time points.
#' @param conc Numeric vector of concentration measurements corresponding to time points.
#' @param route Administration method specification:
#'   \itemize{
#'     \item "bolus" (default) - Excludes time of maximum concentration (Tmax) point
#'     \item "infusion" - Includes Tmax point in terminal phase evaluation
#'   }
#' @param adj_r_squared_threshold Minimum acceptable adjusted R-squared value for valid
#'   estimation (default = 0.7). Values below this threshold will generate warnings.
#' @param tolerance Threshold for considering adjusted R-squared values statistically
#'   equivalent (default = 1e-4). Used when selecting between fits with similar goodness-of-fit.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Lamdaz}: Estimated terminal elimination rate constant (λ_z), or NA if no valid fit
#'   \item \code{UsedPoints}: Number of data points used in the optimal fit
#'   \item \code{adj.r.squared}: Adjusted R-squared value of the optimal regression
#'   \item \code{message}: Character vector containing diagnostic messages/warnings
#' }
#'
#' @details
#' The algorithm implements the following decision logic:
#' \enumerate{
#'   \item Identifies the time of maximum observed concentration (Tmax)
#'   \item Defines candidate terminal phases starting from the last 3 measurable concentrations
#'   \item Iteratively evaluates longer time spans by including preceding data points
#'   \item For each candidate phase:
#'   \itemize{
#'     \item Performs log-concentration vs. time linear regression
#'     \item Requires negative regression slope (positive λ_z)
#'     \item Calculates adjusted R-squared metric
#'   }
#'   \item Selects the optimal phase based on:
#'   \itemize{
#'     \item Highest adjusted R-squared value
#'     \item When R-squared differences are < tolerance, selects the fit with more points
#'   }
#'   \item Validates final selection against R-squared threshold
#' }
#'
#' @examples
#' # Basic usage
#' time <- c(0.5, 1, 2, 4, 6, 8, 10)
#' conc <- c(12, 8, 5, 3, 2, 1.5, 1)
#' result <- find_best_lambdaz(time, conc)
#'
#' # With infusion route specification
#' result_infusion <- find_best_lambdaz(time, conc, route = "infusion")
#'
#' # Custom threshold settings
#' result_custom <- find_best_lambdaz(time, conc, adj_r_squared_threshold = 0.8, tolerance = 0.001)
#'
#' @export
find_best_lambdaz <- function(time,
                       conc,
                       route = "bolus",
                       adj_r_squared_threshold = 0.7,
                       nlastpoints = 3,
                       tolerance = 1e-4) {
  # Initialize variables
  best_lamdaz <- NA
  best_points <- NULL
  best_r2 <- -Inf
  last_best_msg <- NULL
  warn_msgs <- character(0)

  n_points <- length(time)
  tmax_idx <- which.max(conc)

  if (route == "bolus"){
    search_start <-  tmax_idx + 1
  } else {search_start <- tmax_idx}

  # search_start<-tmax_idx
  max_possible_points <- n_points - search_start + 1

  #----- Early Checks -----#
  if (n_points < nlastpoints) {
    return(
      list(
        lamdaz = NA,
        UsedPoints = NULL,
        adj.r.squared = NA,
        message = "ERROR: Insufficient data points",
        slopefit =NULL
      )
    )
  }

  if (max_possible_points < nlastpoints) {
    return(
      list(
        lamdaz = NA,
        UsedPoints = NULL,
        adj.r.squared = NA,
        message = "ERROR: Insufficient terminal phase points",
        slopefit =NULL
      )
    )
  }

  #----- Main Loop -----#
  for (point_count in nlastpoints:max_possible_points) {
    start_idx <- n_points - point_count + 1
    subset <- start_idx:n_points

    # Regression attempt
    fit <- try(lm(log(conc[subset]) ~ time[subset]), silent = TRUE)
    if (inherits(fit, "try-error")) {
      warn_msgs <-
        c(warn_msgs,
          paste(point_count, "points: regression failed"))
      next
    }

    # Extract parameters
    s <- summary(fit)
    current_r2 <- s$adj.r.squared
    current_slope <- coef(fit)[2]

    # Catch invalid R² or slope
    if (is.na(current_r2) || is.nan(current_r2) || is.na(current_slope) || is.nan(current_slope)) {
      warn_msgs <- c(warn_msgs, paste(point_count, "points: invalid regression (NaN or NA)"))
      next
    }

    # Slope validation
    if (current_slope >= 0) {
      warn_msgs <-
        c(warn_msgs,
          paste(point_count, "points: non-negative slope"))
      next
    }

    # Update criteria
    update <- FALSE
    if (current_r2 > best_r2 + tolerance) {
      update <- TRUE
      reason <- "higher R²"
    } else if (abs(current_r2 - best_r2) < tolerance &&
               point_count > best_points) {
      update <- TRUE
      reason <- "more points"
    }

    if (update) {
      best_lamdaz <- as.numeric(-current_slope)
      best_points <- point_count
      best_r2 <- current_r2
      best_fit <- fit
      last_best_msg <- sprintf(
        "Selected %d points (%s) R²=%.4f λz=%.4f",
        point_count,
        reason,
        current_r2,
        as.numeric(best_lamdaz)
      )
    }
  }

  #----- Final Messages -----#
  final_msgs <- character(0)
  if (!is.null(last_best_msg))
    final_msgs <- c(final_msgs, last_best_msg)
  final_msgs <- c(final_msgs, warn_msgs)

  if (is.na(best_lamdaz)) {
    final_msgs <-
      c(final_msgs, "ERROR: No valid elimination phase found")
  } else if (best_r2 < adj_r_squared_threshold) {
    final_msgs <- c(
      final_msgs,
      sprintf(
        "WARNING: R² %.4f < threshold %.2f",
        best_r2,
        adj_r_squared_threshold
      )
    )
  }

  list(
    lamdaz = best_lamdaz,
    UsedPoints = best_points,
    adj.r.squared = best_r2,
    message = if (length(final_msgs) > 0) paste(final_msgs, collapse = "\n"),
    slopefit = best_fit

  )


}
