#' Perform non-compartmental pharmacokinetic analysis
#'
#' Calculates key pharmacokinetic parameters using non-compartmental methods for both intravenous and oral administration data.
#'
#' @param x Numeric vector of observation times.
#' @param y Numeric vector of drug concentration measurements.
#' @param dose Administered dose (default = 1).
#' @param trapezoidal.rule Method for AUC calculation:
#' \itemize{
#'   \item \code{"linear"} - Linear trapezoidal method (default)
#'   \item \code{"linear_up_log_down"} - Linear-up/log-down method (linear for ascending concentrations, logarithmic for descending)
#' }
#' @param ss Steady-state flag:
#' \itemize{
#'   \item 0 - Use extrapolated \eqn{AUC_{0 \rightarrow \infty}} for clearance (default)
#'   \item 1 - Use observed \eqn{AUC_{0 \rightarrow \mathrm{last}}} for clearance
#' }
#' @param nlastpoints Number of terminal points for slope estimation (default = 3).
#' @param slope.method Method for estimating terminal slope (\eqn{\lambda_z}):
#' \itemize{
#'   \item \code{"bestfitforce"} - Force estimation using decreasing number of terminal points if best-fit fails (default)
#'   \item \code{"bestfit"} - Use automated best-fit selection based on adjusted R-squared
#' }
#' @param duration Infusion duration (required if \code{route = "infusion"}).
#' @param route Administration route:
#' \itemize{
#'   \item \code{"bolus"} - Intravenous bolus (default)
#'   \item \code{"oral"} - Oral administration
#'   \item \code{"infusion"} - Intravenous infusion
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item cl - Clearance (CL), calculated as Dose/AUC
#'   \item vz - volume of distribution (Vz), calculated as CL / lambdaz
#'   \item half_life - Terminal elimination half-life, computed as ln(2) / lambdaz
#'   \item auct - Area under the concentration–time curve from time 0 to last measurable concentration
#'   \item auc0_inf - AUC extrapolated to infinity
#'   \item C_last - Last non-zero measurable concentration
#'   \item lambdaz - Terminal elimination rate constant
#'   \item aumc_0_t - Area under the first moment curve from time 0 to last measurable concentration
#'   \item aumc_0_inf - AUMC extrapolated to infinity
#'   \item used_points - Number of time–concentration points used to estimate lambdaz
#'   \item adj.r.squared - Adjusted R-squared of the terminal phase regression
#'   \item messages - Warning or diagnostic messages returned during the calculation
#' }
#' @examples
#' \dontrun{
#' # IV bolus example
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8),
#'                   DV = c(12, 8, 5, 3, 2, 1))
#' getnca(x = dat$TIME, y = dat$DV, dose = 1)
#'
#' # IV infusion example
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8),
#'                   DV = c(2, 8, 5, 3, 2, 1))
#' getnca(x = dat$TIME, y = dat$DV, dose = 1, route = "infusion", duration = 1)
#'
#' # Oral administration example
#' dat <- data.frame(TIME = c(0, 1, 2, 4, 6, 8),
#'                   DV = c(0, 9, 12, 8, 4, 2))
#' getnca(x = dat$TIME, y = dat$DV, route = "oral")
#' }
#'
#' @export
#'
getnca <- function(x,
                   y,
                   dose = 1,
                   trapezoidal.rule = c("linear_up_log_down", "linear"),
                   ss = 0,
                   duration =NULL,
                   nlastpoints = 3,
                   slope.method = c("bestfitforce", "bestfit"),
                   route = c("bolus", "oral", "infusion")) {
  # Record start time
  start.time <- Sys.time()


  # safe check
  trapezoidal.rule <- tryCatch(
    match.arg(trapezoidal.rule, choices = c("linear_up_log_down", "linear")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "trapezoidal_rule",
          as.character(trapezoidal.rule),
          paste(shQuote(c(
            "linear_up_log_down", "linear"
          )), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  route <- tryCatch(
    match.arg(route, choices = c("bolus", "oral", "infusion")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "route",
          as.character(route),
          paste(shQuote(c(
            "bolus", "oral", "infusion"
          )), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  slope.method <- tryCatch(
    match.arg(slope.method, choices = c("bestfitforce","bestfit")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "slope.method",
          as.character(slope.method),
          paste(shQuote(c("bestfitforce","bestfit")), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # --- Defensive programming for infusion ---
  if (route == "infusion" && is.null(duration)) {
    stop("For infusion route, `duration` must be specified.", call. = FALSE)
  }

  # Initialize parameters
  cl         <- NA_real_  # Clearance
  vz         <- NA_real_  # Volume of distribution
  lambdaz    <- NA_real_  # Terminal slope
  half_life  <- NA_real_  # Half-life
  auct       <- NA_real_  # AUC_0–t
  auc0_inf   <- NA_real_  # AUC_0–inf
  C_last     <- NA_real_  # Last measurable conc
  aumc_0_t   <- NA_real_  # AUMC_0–t
  aumc_0_inf <- NA_real_  # AUMC_0–inf

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

    } else if (route == "oral" || route == "infusion") {
      # Oral administration handling
      dat <- rbind(data.frame(TIME = 0, DV = 0), dat)
      dat <- dat[order(dat$TIME),]
    }
  }

  #----- AUC Calculation -----#
  if (trapezoidal.rule == "linear") {
    auct <- trapezoidal_linear(dat$TIME, dat$DV)
    aumc_0_t <- trapezoidal_linear(dat$TIME, dat$DV, moment = T)
  } else if (trapezoidal.rule == "linear_up_log_down") {
    auct <- trapezoidal_linear_up_log_down(dat$TIME, dat$DV)
    aumc_0_t <-
      trapezoidal_linear_up_log_down(dat$TIME, dat$DV, moment = T)
  } else {
    messages <- c(messages, "Invalid trapezoidal rule specified")
  }

  if (slope.method == "bestfit") {
    slope_result <- find_best_lambdaz(
      time = dat$TIME,
      conc = dat$DV,
      route = route,
      duration = duration,
      nlastpoints = nlastpoints
    )
  } else if (slope.method == "bestfitforce") {
    slope_result <- force_find_lambdaz(
      time = dat$TIME,
      conc = dat$DV,
      route = route,
      duration = duration,
      nlastpoints = nlastpoints
    )
  }


  # Terminal Phase Analysis
  ke <- slope_result$lambdaz

  if (!is.na(ke) && ke > 0) {
    half_life <- log(2) / ke
    C_last <- utils::tail(dat$DV[dat$DV > 0], 1)
    t_last <- utils::tail(dat$TIME[dat$DV > 0], 1)

    auc_inf <- C_last / ke
    auc0_inf <- auct + auc_inf

    cl <- ifelse(ss == 0, dose / auc0_inf, dose / auct)
    vz <- cl / ke

    aumct_inf <- (C_last * t_last) / ke + C_last / (ke ^ 2)
    aumc_0_inf <- aumc_0_t + aumct_inf

  } else {
    messages <-
      c(
        messages,
        "Terminal slope (lambdaz) is NA or non-positive; unable to calculate terminal parameters."
      )
    half_life  <- NA_real_
    auc0_inf   <- NA_real_
    cl         <- NA_real_
    vz         <- NA_real_
    aumc_0_inf <- NA_real_

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
    aumc_0_t = aumc_0_t,
    aumc_0_inf = aumc_0_inf,
    used_points = slope_result$UsedPoints,
    adj.r.squared = slope_result$adj.r.squared,
    messages = messages,
    time.spent = time.spent
  )
}

#' Linear trapezoidal rule
#'
#' Computes the area under the curve (AUC) or the area under the moment curve (AUMC) using the linear trapezoidal rule.
#' If \code{moment = TRUE}, the function estimates AUMC by integrating \code{time * concentration}.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration values at each time point.
#' @param moment Logical. If \code{TRUE}, computes AUMC by integrating \code{t * C(t)} instead of just \code{C(t)}.
#'
#' @return A numeric value representing the estimated AUC or AUMC using the linear trapezoidal rule.
#'
#' @examples
#' x <- c(0.5, 1, 2, 4, 6, 8)
#' y <- c(12, 8, 5, 3, 2, 1)
#' trapezoidal_linear(x, y)                # AUC
#' trapezoidal_linear(x, y, moment = TRUE) # AUMC
#'
#' @export

# Linear trapezoidal rule
trapezoidal_linear <- function(x, y, moment = FALSE) {
  # If moment = TRUE, interpret y as concentration and compute t*C(t)
  y1 <- if (moment)
    x * y
  else
    y

  n <- length(x)
  sum(diff(x) * (y1[-n] + y1[-1])) / 2
}

#' Linear-up and log-down trapezoidal rule
#'
#' Computes the area under the curve (AUC) or the area under the moment curve (AUMC)
#' using a hybrid trapezoidal rule. The method uses linear interpolation for increasing
#' or constant concentration segments, and logarithmic interpolation for decreasing segments.
#'
#' If \code{moment = TRUE}, the function calculates the area under the moment curve (AUMC),
#' i.e., it integrates \code{t * C(t)} over time instead of just \code{C(t)}.
#'
#' @param x A numeric vector representing the time points.
#' @param y A numeric vector representing the corresponding concentration values at each time point.
#' @param moment Logical. If \code{TRUE}, computes AUMC by integrating \code{t * C(t)} instead of \code{C(t)}.
#'
#' @return A numeric value representing the estimated AUC or AUMC using the linear-up/log-down trapezoidal method.
#'
#' @examples
#' x <- c(0, 0.5, 1, 2, 4, 6, 8)
#' y <- c(0, 2, 8, 5, 3, 2, 1)
#' trapezoidal_linear_up_log_down(x, y)                # AUC
#' trapezoidal_linear_up_log_down(x, y, moment = TRUE) # AUMC
#'
#' @export

trapezoidal_linear_up_log_down <- function(x, y, moment = FALSE) {
  # delta
  delta_x <- diff(x)
  delta_y <- diff(y)

  # Classify intervals
  linear_up_ <- delta_y > 0
  log_down_  <- delta_y < 0
  constant_  <- delta_y == 0

  # Common values
  y_prev <- y[-length(y)]
  y_curr <- y[-1]
  x_prev <- x[-length(x)]
  x_curr <- x[-1]

  ### --- LINEAR UP --- ###
  if (moment) {
    # moment curve = t * C(t)
    linear_up_auc <- (x_prev[linear_up_] * y_prev[linear_up_] +
                        x_curr[linear_up_] * y_curr[linear_up_]) / 2 * delta_x[linear_up_]

    constant_auc <-
      x_prev[constant_] * y_prev[constant_] * delta_x[constant_]
  } else {
    linear_up_auc <-
      (y_prev[linear_up_] + y_curr[linear_up_]) / 2 * delta_x[linear_up_]
    constant_auc  <- y_prev[constant_] * delta_x[constant_]
  }

  ### --- LOG DOWN --- ###
  log_down_auc <- numeric(sum(log_down_))
  if (moment) {
    # Correct moment integral formula from NonCompart
    k <- (log(y_prev[log_down_]) - log(y_curr[log_down_])) / delta_x[log_down_]
    log_down_auc <- (
      x_prev[log_down_] * y_prev[log_down_] -
        x_curr[log_down_] * y_curr[log_down_]
    ) / k + (y_prev[log_down_] - y_curr[log_down_]) / (k^2)
  } else {
    # Standard log trapezoid AUC
    log_down_auc <- ((y_prev[log_down_] - y_curr[log_down_]) /
                       log(y_prev[log_down_] / y_curr[log_down_])) * delta_x[log_down_]
  }

  # Sum all
  total_auc <- sum(linear_up_auc, na.rm = TRUE) +
    sum(log_down_auc, na.rm = TRUE) +
    sum(constant_auc, na.rm = TRUE)

  return(total_auc)
}
