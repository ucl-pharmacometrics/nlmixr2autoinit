#' Generate pooled data for pharmacokinetic analysis
#'
#' Processes pharmacokinetic (PK) data and produces pooled datasets for based on
#' the specified dosing scenario.
#'
#' @param dat A data frame containing pharmacokinetic data. The required columns depend on `data_type`:
#'   - `"first_dose"`: `dose_number`, `iiobs`, `ID`, `TIME`, `DV`
#'   - `"repeated_doses"`: `EVID`, `AMT`, `tad`, `dose_number`, `iiobs`, `ID`, `TIME`, `DV`, `resetflag`
#'   - `"combined_doses"`: `EVID`, `tad`, `dose_number`, `ID`, `TIME`, `DV`, `resetflag`
#'
#' @param data_type Character string indicating the type of data to pool.
#'   One of `"first_dose"`, `"repeated_doses"`, or `"combined_doses"`.
#'   - `"first_dose"` processes only the first dose of each subject.
#'   - `"repeated_doses"` processes doses beyond the first and their post-dose observations.
#'   - `"combined_doses"` merges both first and repeated dose data into a single dataset for pooled analysis.
#'
#' @param ... Additional arguments passed to `bin.time()`, such as `nbins` or `bin.method`.
#'
#' @details
#' For `"repeated_doses"` and `"combined_doses"` analyses, the function calculates the most
#' common dosing interval (`most_common_ii`) from dosing events (`EVID == 1`) across subjects.
#' This is done by computing the time difference between successive doses within each subject
#' and selecting the mode (most frequent value) of the rounded intervals.
#'
#' This interval is used to define a time window for pooling observed concentrations
#' (`tad <= most_common_ii * 1.2`) to ensure data consistency within a typical dosing period.
#' This approach helps filter out excessively long or irregular post-dose observations
#' that may bias the pooled profile, especially under real-world or non-ideal dosing conditions.
#'
#' @return A list containing one or more elements, depending on the specified `data_type`:
#'   \describe{
#'     \item{`datpooled_fd`}{Binned data corresponding to the first dose only, typically used for single-dose NCA.}
#'     \item{`datpooled_efd`}{Binned data excluding the first dose, representing repeated-dose profiles.}
#'     \item{`datpooled_all`}{Binned data that combines both first and repeated doses into a single dataset,
#'                            used when `data_type = "combined_doses"` for pooled overall profile analysis.}
#'   }
#'
#' @examples
#' \dontrun{
#' processed_data <- processData(Bolus_1CPT)$dat
#' pooled <- get_pooled_data(processed_data, data_type = "combined_doses")
#' pooled$datpooled_all$binned.df
#' }
#'
#' @export
#'
get_pooled_data <- function(dat, data_type = "first_dose", ...) {
  dots <- list(...)
  bin_args <- dots[names(dots) %in% names(formals(bin.time))]

  valid_types <- c("first_dose", "repeated_doses", "combined_doses")
  if (!data_type %in% valid_types) {
    stop(
      "Invalid data_type: '",
      data_type,
      "'. Valid options are: ",
      paste(valid_types, collapse = ", ")
    )
  }

  pooled_results <- list(
    datpooled_fd = NA,
    datpooled_efd = NA,
    datpooled_all = NA
  )

  if (data_type == "first_dose" || data_type == "combined_doses") {
    fd_data <- dat[dat$dose_number == 1 & dat$iiobs == 0,]
    fd_data <- fd_data[order(fd_data$ID, fd_data$TIME),]
    pooled_results$datpooled_fd <-
      do.call(bin.time, c(list(dat = fd_data), bin_args))
  }

  if (data_type == "repeated_doses" ||
      data_type == "combined_doses") {
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

    efd_data <- dat %>%
      dplyr::filter(
        (dose_number != 1 & EVID == 0 & tad <= most_common_ii * 1.2) |
          (dose_number != 1 & EVID == 1) |
          (dose_number == 1 & iiobs > 0)
      ) %>%
      dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

    pooled_results$datpooled_efd <-
      do.call(bin.time, c(list(dat = efd_data), bin_args))

    if (data_type == "combined_doses") {
      all_data <- dat %>%
        dplyr::filter((EVID == 0 &
                         tad <= most_common_ii * 1.2) |
                        EVID == 1) %>%
        dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

      pooled_results$datpooled_all <-
        do.call(bin.time, c(list(dat = all_data), bin_args))
    }
  }

  return(pooled_results)
}
