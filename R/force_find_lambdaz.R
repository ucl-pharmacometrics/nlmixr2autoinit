#' Forceful estimation of terminal slope
#'
#' Estimates the terminal elimination rate constant (lambda_z) of a pharmacokinetic
#' profile. The function first attempts to use the `find_best_lambdaz` method. If no
#' valid estimate is obtained, it falls back to a simplified log-linear regression
#' using progressively fewer data points to enforce a negative slope.
#'
#' @param time Numeric vector of time points.
#' @param conc Numeric vector of concentration values corresponding to time.
#' @param ... Additional arguments passed to find_best_lambdaz (e.g., nlastpoints).
#'
#' @details
#' This function implements a two-step strategy to ensure estimation of the
#' terminal elimination slope:
#'   - First, it applies `find_best_lambdaz` to automatically select the best
#'     fitting terminal phase segment based on adjusted R-squared optimization.
#'   - If `find_best_lambdaz` fails (e.g., limited data), the function forcibly
#'     fits simplified linear models using progressively fewer points (starting
#'     from n-1 down to 2) until a negative slope is identified. In fallback
#'     mode, adjusted R-squared is not considered.
#'
#' @return A list containing:
#'   - lambdaz: Estimated terminal elimination rate constant (1/time)
#'   - intercept: Intercept of the log-linear regression, used to extrapolate concentration at time zero
#'   - method: Method used (`find_best_lambdaz` or fallback regression)
#'   - UsedPoints: Number of time-concentration points used for estimation
#'   - adj.r.squared: Adjusted R-squared (available only when using `find_best_lambdaz`)
#'   - message: Diagnostic message summarizing the outcome
#'   - slopefit: Fitted linear model object
#'
#' @seealso \link{find_best_lambdaz}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' time <- c(0.5, 1, 2, 4, 6, 8, 10)
#' conc <- c(12, 8, 5, 3, 2, 1.5, 1)
#' force_find_lambdaz(time, conc)
#'
#' @export
#'
force_find_lambdaz <- function(time, conc, ...) {
  dots <- list(...)
  slope_args <-
    dots[names(dots) %in% names(formals(find_best_lambdaz))]

  result <-
    do.call(find_best_lambdaz, c(list(time = time, conc = conc),
                                 slope_args))

  lambdaz <- result$lambdaz
  method <- NA
  intercept <- NA
  used_points <- result$UsedPoints
  adj_r2 <- result$adj.r.squared
  message_text <- result$message
  slopefit <- result$slopefit

  if (!is.na(lambdaz)) {
    method <- "find_best_lambdaz"
    if (!is.null(slopefit)) {
      intercept <- summary(slopefit)[[4]][[1]]
    }
  } else {
    if (length(time) == 2) {
      fallback_points <- data.frame(time = time, conc = conc)
      fit <-
        try(lm(log(fallback_points$conc) ~ fallback_points$time), silent = TRUE)

      if (!inherits(fit, "try-error")) {
        coefs <- summary(fit)[[4]]
        slope_val <- coefs[2]
        intercept_val <- coefs[1]

        if (!is.na(slope_val) && slope_val < 0) {
          method <- "fallback_regression"
          lambdaz <- -slope_val
          intercept <- intercept_val
          used_points <- nrow(fallback_points)
          adj_r2 <- NA
          message_text <-
            "Fallback regression successful with 2 points."
          slopefit <- fit
        }
      }
    } else {
      for (k in seq(length(time) - 1, 2, by = -1)) {
        fallback_points <-
          utils::tail(data.frame(time = time, conc = conc), k)
        fit <-
          try(lm(log(fallback_points$conc) ~ fallback_points$time), silent = TRUE)
        if (inherits(fit, "try-error"))
          next

        coefs <- summary(fit)[[4]]
        slope_val <- coefs[2]
        intercept_val <- coefs[1]

        if (!is.na(slope_val) && slope_val < 0) {
          method <- "fallback_regression"
          lambdaz <- -slope_val
          intercept <- intercept_val
          used_points <- nrow(fallback_points)
          adj_r2 <- NA
          message_text <- "Fallback regression successful."
          slopefit <- fit
          break
        }
      }
    }
  }

  return(
    list(
      lambdaz = lambdaz,
      intercept = intercept,
      method = method,
      UsedPoints = used_points,
      adj.r.squared = adj_r2,
      message = message_text,
      slopefit = slopefit
    )
  )
}
