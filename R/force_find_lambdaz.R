#' Forceful estimation of terminal slope (Lambda_z)
#'
#' Estimates the terminal slope (lambda_z) of a pharmacokinetic profile
#' by first attempting to use the \code{find_best_lambdaz} method. If no valid lambda_z
#' is identified, it falls back to simplified linear regression (log concentration vs time)
#' using a decreasing number of points to forcibly find an acceptable estimate.
#'
#' @param time Numeric vector of time points.
#' @param conc Numeric vector of concentration measurements corresponding to \code{time}.
#' @param ... Additional arguments passed to \code{find_best_lambdaz}, such as \code{nlastpoints}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{lambdaz}}{Estimated terminal elimination rate constant (lambda_z, 1/time units)}
#'   \item{\code{intercept}}{Intercept of the log-linear regression, used to extrapolate C0}
#'   \item{\code{method}}{Method used: \code{\"find_best_lambdaz\"} or \code{\"fallback_regression\"}}
#'   \item{\code{UsedPoints}}{Time-concentration points used for terminal slope estimation}
#'   \item{\code{adj.r.squared}}{Adjusted R-squared (only available when using find_best_lambdaz, otherwise NA)}
#'   \item{\code{message}}{Text message indicating fit diagnostics}
#'   \item{\code{slopefit}}{The fitted lm() object from regression}
#' }
#'
#' @details
#' This function implements a two-step strategy to ensure estimation of the terminal elimination slope:
#' \enumerate{
#'   \item First, it applies \code{find_best_lambdaz} to automatically select the best fitting terminal phase segment
#'   based on adjusted R-squared optimization.
#'   \item If \code{find_best_lambdaz} fails (e.g., poor data quality or lack of clear terminal phase), the function
#'   forcibly fits simplified linear models using progressively fewer points (starting from \code{n-1} down to 2)
#'   until a negative slope is identified. In fallback mode, adjusted R-squared is not considered.
#' }
#'
#' Application scenarios include:
#' \itemize{
#'   \item Exploratory pharmacokinetic (PK) data analysis
#'   \item Analysis of datasets with limited or poor quality data
#'   \item Automated processing pipelines where manual curve fitting is impractical
#'   \item Cases where half-life estimation is required for reporting or modeling purposes
#' }
#'
#' @seealso \code{\link{find_best_lambdaz}}
#'
#' @examples
#' # Basic usage
#'
#' time <- c(0.5, 1, 2, 4, 6, 8, 10)
#' conc <- c(12, 8, 5, 3, 2, 1.5, 1)
#' find_best_lambdaz(time, conc)
#'
#' time <- c(0.5, 1, 2)
#' conc <- c(12, 8, 5)
#' force_find_lambdaz(time, conc)
#'
#'
#' @export
#'
#'
force_find_lambdaz <- function(time, conc, ...) {

  dots <- list(...)
  slope_args <-
    dots[names(dots) %in% names(formals(find_best_lambdaz))]

  result <- do.call(find_best_lambdaz, c(list(time = time, conc = conc),
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
      fit <- try(lm(log(fallback_points$conc) ~ fallback_points$time), silent = TRUE)

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
          message_text <- "Fallback regression successful with 2 points."
          slopefit <- fit
        }
      }
    } else {
      for (k in seq(length(time) - 1, 2, by = -1)) {
        fallback_points <- tail(data.frame(time = time, conc = conc), k)
        fit <- try(lm(log(fallback_points$conc) ~ fallback_points$time), silent = TRUE)
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

