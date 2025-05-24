#' Control settings for pooled analysis
#'
#' Controls time binning and preprocessing used in pooled pharmacokinetic analysis.
#'
#' @param nbins Integer or "auto". Number of time bins. Default is 10.
#' @param bin_method Character. Binning method to use. One of:
#'   "quantile", "jenks", "kmeans", "pretty", "sd", "equal", "density".
#'
#' @return A named list of pooled control parameters.
#' @export
#'
#' @examples
#' pooled_control(nbins = 8, bin_method = "jenks")
pooled_control <- function(nbins = 10,
                           bin_method = c("quantile",
                                          "jenks",
                                          "kmeans",
                                          "pretty",
                                          "sd",
                                          "equal",
                                          "density"),
                           tad_rounding = TRUE) {
  bin_method <- tryCatch(
    match.arg(
      bin_method,
      choices = c(
        "quantile",
        "jenks",
        "kmeans",
        "pretty",
        "sd",
        "equal",
        "density"
      )
    ),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `bin_method`: '%s'. Must be one of: %s.",
          as.character(bin_method),
          paste(shQuote(
            c(
              "quantile",
              "jenks",
              "kmeans",
              "pretty",
              "sd",
              "equal",
              "density"
            )
          ), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # Defensive check for nbins
  if (!(is.numeric(nbins) &&
        nbins > 0 && nbins == as.integer(nbins)) &&
      !(is.character(nbins) && nbins == "auto")) {
    stop(
      "`nbins` must be a positive integer or the string 'auto'. ",
      "You provided: ",
      deparse(substitute(nbins)),
      " = ",
      nbins
    )
  }

  # Safe logical check
  if (!is.logical(tad_rounding) ||
      length(tad_rounding) != 1 || is.na(tad_rounding)) {
    stop("`tad_rounding` must be a single logical value: TRUE or FALSE.",
         call. = FALSE)
  }

  list(nbins = nbins,
       bin_method = bin_method,
       tad_rounding = tad_rounding)
}


#' Generate pooled data for pharmacokinetic analysis
#'
#' Processes pharmacokinetic (PK) data and produces pooled datasets based on
#' the specified dosing scenario, with flexible control over binning and time-alignment behavior.
#'
#' @param dat A data frame containing pharmacokinetic data. The required columns depend on `dose_type`:
#'   - `"first_dose"`: `dose_number`, `iiobs`, `ID`, `TIME`, `DV`
#'   - `"repeated_doses"`: `EVID`, `AMT`, `tad`, `dose_number`, `iiobs`, `ID`, `TIME`, `DV`, `resetflag`
#'   - `"combined_doses"`: `EVID`, `tad`, `dose_number`, `ID`, `TIME`, `DV`, `resetflag`
#'
#' @param dose_type Character string indicating the type of data to pool.
#'   One of `"first_dose"`, `"repeated_doses"`, or `"combined_doses"`.
#'   - `"first_dose"` processes only the first dose of each subject.
#'   - `"repeated_doses"` processes doses beyond the first and their post-dose observations.
#'   - `"combined_doses"` merges both first and repeated dose data into a single dataset for pooled analysis.
#'
#' @param pooled_ctrl A list of control options created by `pooled_control()`. Includes:
#'   \describe{
#'     \item{`nbins`}{Number of bins to use when binning data.}
#'     \item{`bin_method`}{Method used for time binning (e.g., `"equal"`, `"quantile"`).}
#'     \item{`tad_rounding`}{Logical. If TRUE (default), both `tad` and the most common dosing interval
#'                           are rounded to the nearest whole unit before comparing. This allows for
#'                           small deviations (e.g., a `tad` of 24.3 is treated as within a 24-unit interval).}
#'   }
#'
#' @param ... Additional arguments passed to `bin.time()`, such as `nbins` or `bin.method`.
#'
#' @details
#' For `"repeated_doses"` and `"combined_doses"` data types, the function determines the most
#' common dosing interval (`most_common_ii`) based on dosing records (`EVID == 1`). This interval is
#' used to assess whether post-dose observations fall within the relevant dosing window.
#' If `tad_rounding` is TRUE, both the observed time after dose (`tad`) and the interval are rounded
#' before comparison, making the method more robust to minor time variations.
#'
#' The function then filters, bins, and returns the appropriate datasets for analysis.
#'
#' @return A list containing one or more of the following elements, depending on the `dose_type`:
#'   \describe{
#'     \item{`datpooled_fd`}{Binned data for the first dose (single-dose profiles).}
#'     \item{`datpooled_efd`}{Binned data excluding the first dose (repeated-dose profiles).}
#'     \item{`datpooled_all`}{Combined first and repeated dose data, for overall pooled analysis.}
#'   }
#'
#' @examples
#' \dontrun{
#' processed_data <- processData(Bolus_1CPT)$dat
#' pooled <- get_pooled_data(processed_data, dose_type = "combined_doses")
#' pooled$datpooled_all$binned.df
#' }
#'
#' @export

get_pooled_data <- function(dat,
                            dose_type = c("first_dose", "repeated_doses", "combined_doses"),
                            pooled_ctrl = pooled_control()) {
  # Safe matching
  dose_type <- tryCatch(
    match.arg(
      dose_type,
      choices = c("first_dose", "repeated_doses", "combined_doses")
    ),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "dose_type",
          as.character(dose_type),
          paste(shQuote(
            c("first_dose", "repeated_doses", "combined_doses")
          ), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # Extract binning parameters
  bin_args <- list(nbins = pooled_ctrl$nbins,
                   bin.method = pooled_ctrl$bin_method)

  pooled_results <- list(
    datpooled_fd = NA,
    datpooled_efd = NA,
    datpooled_all = NA
  )

  if (dose_type == "first_dose" || dose_type == "combined_doses") {
    fd_data <- dat[dat$dose_number == 1 & dat$iiobs == 0, ]
    fd_data <- fd_data[order(fd_data$ID, fd_data$TIME), ]
    pooled_results$datpooled_fd <-
      do.call(bin.time, c(list(dat = fd_data), bin_args))
  }

  if (dose_type == "repeated_doses" ||
      dose_type == "combined_doses") {
    most_common_ii <- dat %>%
      dplyr::filter(EVID == 1, AMT > 0) %>%
      dplyr::arrange(ID, TIME) %>%
      dplyr::group_by(ID) %>%
      dplyr::mutate(interval = c(NA, diff(TIME))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(interval)) %>%
      dplyr::mutate(interval = round(interval)) %>%
      dplyr::count(interval, sort = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::pull(interval) %>%
      as.numeric()

    # Compute tad_check based on logical control
    dat <- dat %>%
      dplyr::mutate(tad_check = if (pooled_ctrl$tad_rounding) {
        round(tad) <= round(most_common_ii)
      } else {
        tad <= most_common_ii
      })

    # Filter based on tad_check
    efd_data <- dat %>%
      dplyr::filter(
        (dose_number != 1 & EVID == 0 & tad_check) |
          (dose_number == 1 & iiobs > 0 & tad_check)
      ) %>%
      dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

    pooled_results$datpooled_efd <-
      do.call(bin.time, c(list(dat = efd_data), bin_args))

    if (dose_type == "combined_doses") {
      all_data <- dat %>%
        dplyr::filter((EVID == 0 & tad_check) | EVID == 1) %>%
        dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

      pooled_results$datpooled_all <-
        do.call(bin.time, c(list(dat = all_data), bin_args))
    }
  }

  return(pooled_results)
}




