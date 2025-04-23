#' Perform non-compartmental pharmacokinetic analysis
#'
#' Calculates key pharmacokinetic parameters using non-compartmental methods for both intravenous and oral administration data.
#'
#' @param x Numeric vector of observation times
#' @param y Numeric vector of drug concentration measurements
#' @param dose Administered dose (default = 1)
#' @param trapezoidal.rule Method for AUC calculation:
#' \itemize{
#'   \item 1 - Linear trapezoidal method (default)
#'   \item 2 - Linear-up/log-down method (linear for ascending concentrations, logarithmic for descending)
#' }
#' @param ss Steady-state flag:
#' \itemize{
#'   \item 0 - Use extrapolated AUC₀→∞ for clearance (default)
#'   \item 1 - Use observed AUC₀→last for clearance
#' }
#' @param nlastpoints Number of terminal points for slope estimation (default = 3)
#' @param route Administration route:
#' \itemize{
#'   \item "bolus" - Intravenous bolus (default)
#'   \item "oral" - Oral administration
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item cl - Clearance (CL)
#'   \item vz - Volume of distribution (Vz)
#'   \item half_life - Terminal half-life
#'   \item auct - AUC₀→last
#'   \item auc0_inf - AUC₀→∞
#'   \item C_last - Last measurable concentration
#'   \item lamdaz - Terminal slope (λ_z)
#'   \item aumc_0_inf - AUMC₀→∞
#'   \item used_points - Number of points used for terminal slope
#'   \item adj.r.squared - Adjusted R² for terminal slope
#'   \item messages - Concatenated warning/status messages
#' }
#'
#' @examples
#' # IV bolus example
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8),
#'                   DV = c(12, 8, 5, 3, 2, 1))
#' getnca(x = dat$TIME, y = dat$DV, dose = 1)
#'
#' # Oral administration example
#' oral_dat <- data.frame(TIME = c(0, 1, 2, 4, 6, 8),
#'                        DV = c(0, 9, 12, 8, 4, 2))
#' getnca(x = oral_dat$TIME, y = oral_dat$DV, route = "oral")
#'
#' @export

getnca <- function(x,
                   y,
                   dose = 1,
                   trapezoidal.rule = 1,
                   ss = 0,
                   nlastpoints = 3,
                   route = "bolus") {
  # Initialize parameters and messages
  cl <-
    vz <-
    lamdaz <-
    half_life <- auct <- auc0_inf <- C_last <- aumc_0_inf <- NA
  messages <- character(0)

  dat <- data.frame(TIME = x, DV = y)
  dat <- dat[order(dat$TIME),]

  #----- C0 Handling -----#
  if (!0 %in% dat$TIME) {
    if (route == "bolus") {
      # Back-extrapolation logic
      warn_flag <- FALSE
      if (nrow(dat) >= 2 && all(dat$DV[1:2] > 0)) {
        t1 <- dat$TIME[1]
        c1 <- dat$DV[1]
        t2 <- dat$TIME[2]
        c2 <- dat$DV[2]

        slope_est <- (log(c2) - log(c1)) / (t2 - t1)

        if (slope_est < 0 & c1 > c2) {
          # Valid extrapolation
          c0 <- exp(log(c1) - slope_est * t1)
          dat <- rbind(data.frame(TIME = 0, DV = c0), dat)
        } else {
          warn_flag <- TRUE
          msg <- paste(
            "Back-extrapolation failed:",
            ifelse(
              slope_est >= 0,
              "non-negative slope",
              "concentration increase detected"
            )
          )
          messages <- c(messages, msg)
        }
      } else {
        warn_flag <- TRUE
        msg <- paste(
          "Insufficient data for back-extrapolation:",
          ifelse(nrow(dat) < 2, "n < 2",
                 "zero in first two points")
        )
        messages <- c(messages, msg)
      }

      if (warn_flag) {
        # Fallback handling
        dat <- rbind(data.frame(TIME = 0, DV = dat$DV[1]), dat)
        messages <-
          c(messages, "Used first observed concentration as C0")
      }

      dat <- dat[order(dat$TIME),]

    } else if (route == "oral") {
      # Oral administration handling
      dat <- rbind(data.frame(TIME = 0, DV = 0), dat)
      dat <- dat[order(dat$TIME),]
    }
  }

  #----- AUC Calculation -----#
  if (trapezoidal.rule == 1) {
    auct <- trapezoidal_linear(dat$TIME, dat$DV)
  } else if (trapezoidal.rule == 2) {
    auct <- trapezoidal_linear_up_log_down(dat$TIME, dat$DV)
  } else {
    messages <- c(messages, "Invalid trapezoidal rule specified")
  }

  #----- Terminal Phase Analysis -----#
  slope_result <- best_slope(
    time = dat$TIME,
    conc = dat$DV,
    route = route,
    nlastpoints = nlastpoints
  )

  # Process slope results
  ke <- slope_result$lamdaz
  half_life <- ifelse(ke > 0, log(2) / ke, NA)
  C_last <- tail(dat$DV[dat$DV > 0], 1)
  t_last <- tail(dat$TIME[dat$DV > 0], 1)

  #----- PK Parameters -----#
  if (ke > 0) {
    auc_inf <- C_last / ke
    auc0_inf <- auct + auc_inf
    cl <- ifelse(ss == 0, dose / auc0_inf, dose / auct)
    vz <- cl / ke
  } else {
    messages <- c(messages, "Terminal slope calculation failed")
  }

  #----- AUMC Calculation -----#
  moment_curve <- dat$TIME * dat$DV
  aumc0_t <- trapezoidal_linear(dat$TIME, moment_curve)
  aumct_inf <-
    ifelse(ke > 0, (C_last * t_last) / ke + C_last / (ke ^ 2), NA)
  aumc_0_inf <- aumc0_t + aumct_inf

  #----- Message Consolidation -----#
  # Add slope calculation messages
  if (!is.null(slope_result$message)) {
    slope_msgs <- unlist(strsplit(slope_result$message, "\n"))
    messages <- c(messages, slope_msgs)
  }

  # Final message formatting
  if (length(messages) > 0) {
    messages <-
      paste("[Message]", seq_along(messages), messages, sep = ": ")
    messages <- paste(messages, collapse = "\n")
  } else {
    messages <- NA_character_
  }

  # Return structured results
  list(
    clobs = cl,
    vzobs = vz,
    half_life = half_life,
    auct = auct,
    auc0_inf = auc0_inf,
    C_last = C_last,
    lamdaz = slope_result$lamdaz,
    aumc_0_inf = aumc_0_inf,
    used_points = slope_result$UsedPoints,
    adj.r.squared = slope_result$adj.r.squared,
    messages = messages
  )
}

#' Calculate Area Under the Curve (AUC) using the trapezoidal linear rule
#'
#' Computes the area under the curve (AUC) using the linear trapezoidal rule.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration at each time point.
#'
#' @return A numeric value representing the estimated AUC using the trapezoidal rule.
#'
#' @examples
#' x <- c(0.5, 1, 2, 4, 6, 8, 10)
#' y <- c(12, 8, 5, 3, 2, 1.5, 1)
#' trapezoidal_linear(x, y)
#'

# Linear trapezoidal rule
trapezoidal_linear <- function(x, y) {
  n <- length(x)
  sum(diff(x) * (y[-n] + y[-1])) / 2
}

#' Calculate AUC using linear up and log down trapezoidal rule
#'
#' This function computes the area under the curve (AUC) by using the trapezoidal rule
#' for phases where concentration is increasing, and a logarithmic rule for phases
#' where concentration is decreasing.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration at each time point.
#'
#' @return A numeric value representing the estimated AUC using the linear up/log down method.
#'
#' @examples
#' x <- c(0.5, 1, 2, 4, 6, 8, 10)
#' y <- c(12, 8, 5, 3, 2, 1.5, 1)
#' trapezoidal_linear_up_log_down(x, y)
#'
#' @export
#
trapezoidal_linear_up_log_down <- function(x, y) {

  delta_x <- diff(x)
  delta_y <- diff(y)

  # Classify intervals
  linear_up_   <- delta_y > 0
  log_down_    <- delta_y < 0
  constant_    <- delta_y == 0

  # For the linear up phase
  linear_up_auc <- (y[-length(y)][linear_up_] + y[-1][linear_up_]) / 2 * delta_x[linear_up_]

  # For the log down phase
  log_down_auc <- ((y[-length(y)][log_down_] - y[-1][log_down_]) /
                     log(y[-length(y)][log_down_] / y[-1][log_down_])) * delta_x[log_down_]

  # Constant components
  constant_auc <- y[-length(y)][constant_] * delta_x[constant_]

  # Sum all components
  total_auc <- sum(linear_up_auc, na.rm = TRUE) +
    sum(log_down_auc, na.rm = TRUE) +
    sum(constant_auc, na.rm = TRUE)

  return(total_auc)
}



#' Estimate terminal elimination rate constant (Lambda Z)
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
#' result <- best_slope(time, conc)
#'
#' # With infusion route specification
#' result_infusion <- best_slope(time, conc, route = "infusion")
#'
#' # Custom threshold settings
#' result_custom <- best_slope(time, conc, adj_r_squared_threshold = 0.8, tolerance = 0.001)
#'
#' @export
best_slope <- function(time,
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
  search_start <- if (route == "bolus")
    tmax_idx + 1
  else
    tmax_idx
  max_possible_points <- n_points - search_start + 1

  #----- Early Checks -----#
  if (n_points < nlastpoints) {
    return(
      list(
        lamdaz = NA,
        UsedPoints = NULL,
        adj.r.squared = NA,
        message = "ERROR: Insufficient data points"
      )
    )
  }

  if (max_possible_points < nlastpoints) {
    return(
      list(
        lamdaz = NA,
        UsedPoints = NULL,
        adj.r.squared = NA,
        message = "ERROR: Insufficient terminal phase points"
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
    message = if (length(final_msgs) > 0)
      paste(final_msgs, collapse = "\n")
  )


}
