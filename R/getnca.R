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
#'   \item lambdaz - Terminal slope (λ_z)
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
  # Record start time
  start.time <- Sys.time()

  # Initialize parameters and messages
  cl <-
    vz <-
    lambdaz <-
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
  slope_result <- find_best_lambdaz(
    time = dat$TIME,
    conc = dat$DV,
    route = route,
    nlastpoints = nlastpoints
  )

  # Terminal Phase Analysis
  ke <- slope_result$lambdaz

  if (!is.na(ke) && ke > 0) {
    half_life <- log(2) / ke
    C_last <- tail(dat$DV[dat$DV > 0], 1)
    t_last <- tail(dat$TIME[dat$DV > 0], 1)

    auc_inf <- C_last / ke
    auc0_inf <- auct + auc_inf
    cl <- ifelse(ss == 0, dose / auc0_inf, dose / auct)
    vz <- cl / ke

  #----- AUMC Calculation -----#
    moment_curve <- dat$TIME * dat$DV
    aumc0_t <- trapezoidal_linear(dat$TIME, moment_curve)
    aumct_inf <- (C_last * t_last) / ke + C_last / (ke ^ 2)
    aumc_0_inf <- aumc0_t + aumct_inf

  } else {
    messages <- c(messages, "Terminal slope (lambdaz) is NA or non-positive; unable to calculate terminal parameters.")
    half_life <- NA
    auc0_inf <- NA
    cl <- NA
    vz <- NA
    aumc_0_inf <- NA
  }

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

  end.time <- Sys.time()
  time.spent <-
    round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

  # Return structured results
  list(
    clobs = cl,
    vzobs = vz,
    half_life = half_life,
    auct = auct,
    auc0_inf = auc0_inf,
    C_last = C_last,
    lambdaz = slope_result$lambdaz,
    aumc_0_inf = aumc_0_inf,
    used_points = slope_result$UsedPoints,
    adj.r.squared = slope_result$adj.r.squared,
    messages = messages,
    time.spent = time.spent
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



