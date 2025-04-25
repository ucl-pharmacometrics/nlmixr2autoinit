#' Estimate half-life from pooled pharmacokinetic data
#'
#' Estimates the terminal half-life of a drug using pooled and binned
#' pharmacokinetic data. It supports analysis based on first-dose data, repeated-dose data,
#' or a combined profile that includes both. The estimation is performed by applying
#' linear regression on the terminal phase of log-transformed concentration-time data.
#'
#' @param dat A data frame containing the pharmacokinetic data. The required columns
#'   depend on the specified \code{data_type} and are validated within \code{\link{get_pooled_data}}.
#' @param data_type Character string specifying the type of data to analyze. One of:
#'   \itemize{
#'     \item{\code{"first_dose"}}: Analysis based on first-dose data only.
#'     \item{\code{"repeated_doses"}}: Analysis based on repeated-dose data only.
#'     \item{\code{"combined_doses"}}: Combined analysis using both first and repeated doses.
#'   }
#' @param pooled (Optional) A pooled data object as returned by \code{\link{get_pooled_data}}.
#'   If not supplied, the function will internally generate pooled data using \code{dat}
#'   and any applicable arguments to \code{bin.time}.
#' @param ... Additional arguments passed to either \code{\link{bin.time}} (for pooling)
#'   or \code{\link{find_best_lambdaz}} (for elimination slope calculation). The function
#'   automatically splits these arguments based on their intended function.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item{\code{half_life_median}}: The median of all positive half-life estimates from the subsets.
#'     \item{\code{half_life_fd}}: The half-life estimate based on first-dose data, if available.
#'     \item{\code{half_life_md}}: The half-life estimate based on repeated-dose data, if available.
#'     \item{\code{half_life_all}}: The half-life estimate based on the full dataset, if available.
#'   }
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Generating (or using) pooled data using \code{\link{get_pooled_data}}, which applies
#'         time binning to normalize the time-course concentration data.
#'   \item Estimating the terminal elimination slope (\code{lamdaz}) using \code{\link{find_best_lambdaz}}.
#'   \item Computing the terminal half-life as \code{log(2) / lamdaz} for each data subset.
#'   \item Returning the median of all valid (positive, non-missing) half-life estimates.
#' }
#'
#' @examples
#' \dontrun{
#' # Example: half-life estimation from a combined profile (first + repeated doses)
#' dat <- Bolus_1CPT
#' dat <- processData(dat)$dat
#'
#' # For half-life estimation, the number of terminal points for regression can start from 2.
#' # (NCA methods use 3 or more points.)
#' get_hf(dat, data_type = "combined_doses", nlastpoints = 3)
#'
#' # Using externally pooled data
#' pooled <- get_pooled_data(dat, data_type = "combined_doses")
#' get_hf(dat, data_type = "combined_doses", pooled = pooled, nlastpoints = 2)
#' }
#'
#' @seealso \code{\link{get_pooled_data}}, \code{\link{bin.time}}, \code{\link{find_best_lambdaz}}
#' @importFrom stats lm
#' @export
#'
get_hf <- function(dat,
                   data_type = "first_dose",
                   pooled = NULL,
                   ...) {

  message(crayon::black(paste0("Estimating half-life", strrep(".", 20))))

  dots <- list(...)
  bin_args <- dots[names(dots) %in% names(formals(bin.time))]
  slope_args <- dots[names(dots) %in% names(formals(find_best_lambdaz))]

  if (is.null(pooled)) {
    pooled <- do.call(get_pooled_data,
                      c(list(
                        dat = dat, data_type = data_type
                      ), bin_args))
  }

  # Initialize results
  half_life_fd <- NA
  half_life_md <- NA
  half_life_all <- NA

  # First dose
  if (!is.null(pooled$datpooled_fd) &&
      "binned.df" %in% names(pooled$datpooled_fd)) {
    slope_results <- do.call(find_best_lambdaz, c(
      list(
        time = pooled$datpooled_fd$binned.df$Time,
        conc = pooled$datpooled_fd$binned.df$Conc
      ),
      slope_args
    ))
    ke <- slope_results$lamdaz
    half_life_fd <- ifelse(ke > 0, log(2) / ke, NA)
  }

  # Repeated doses
  if (!is.null(pooled$datpooled_efd) &&
      "binned.df" %in% names(pooled$datpooled_efd)) {
    slope_results <- do.call(find_best_lambdaz, c(
      list(
        time = pooled$datpooled_efd$binned.df$Time,
        conc = pooled$datpooled_efd$binned.df$Conc
      ),
      slope_args
    ))
    ke <- slope_results$lamdaz
    half_life_md <- ifelse(ke > 0, log(2) / ke, NA)
  }

  # Combined
  if (!is.null(pooled$datpooled_all) &&
      "binned.df" %in% names(pooled$datpooled_all)) {
    slope_results <-  do.call(find_best_lambdaz, c(
      list(
        time = pooled$datpooled_all$binned.df$Time,
        conc = pooled$datpooled_all$binned.df$Conc
      ),
      slope_args
    ))
    ke <- slope_results$lamdaz
    half_life_all <- ifelse(ke > 0, log(2) / ke, NA)
  }

  # Summarize
  half_life_values <- c(half_life_fd, half_life_md, half_life_all)
  positive_values <-
    half_life_values[half_life_values > 0 & !is.na(half_life_values)]
  half_life_median <-
    if (length(positive_values) > 0)
      round(median(positive_values), 2)
  else
    NA

  return(
    list(
      half_life_median = half_life_median,
      half_life_fd = half_life_fd,
      half_life_md = half_life_md,
      half_life_all = half_life_all
    )
  )
}
