#' Control settings for pooled data analysis
#'
#' Defines control parameters for time binning and preprocessing in pooled
#' data analysis. These parameters are typically passed to 'get_pooled_data'.
#'
#' @param nbins Integer or the character string auto. Number of time bins used
#'   to group observations. Default is 10.
#' @param bin_method Character string specifying the binning method. Must be one
#'   of "quantile", "jenks", "kmeans", "pretty", "sd", "equal", or "density".
#' @param tad_rounding Logical value indicating whether tad and the most common
#'   dosing interval should be rounded to the nearest whole unit before
#'   comparison. Default is TRUE, allowing small deviations (for example, a tad
#'   of 24.3 is treated as within a 24-unit interval).
#'
#' @return
#' A named list containing control parameters for pooled pharmacokinetic
#' analysis.
#'
#' @seealso \link{get_pooled_data}
#'
#' @examples
#' pooled_control()
#'
#' @export


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
#' Processes pharmacokinetic data and produces pooled datasets according to
#' the dosing context. Data can be grouped based on first dose, repeated dosing,
#' or a combination of both, with control over binning and time alignment.
#'
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @param dose_type Specifies the dosing context of the pharmacokinetic
#'   observations. Classified as:
#'   \itemize{
#'     \item first_dose: data include only observations following the initial
#'           administration
#'     \item repeated_doses: data include only observations during repeated or
#'           steady-state dosing
#'     \item combined_doses: data include observations from both first-dose and
#'           repeated-dose intervals
#'   }
#'
#' @param pooled_ctrl A list of control parameters created by
#'   'pooled_control', including settings for binning and time rounding.
#'
#' @details
#' For repeated-doses and combined-doses classifications, the most common
#' interdose interval is identified from dosing records and used to determine
#' whether observations fall within the relevant interval. If tad_rounding is
#' TRUE, both time after dose and dosing interval are rounded before comparison.
#'
#' @return
#' A list containing pooled pharmacokinetic datasets depending on the specified
#' dose type:
#'   - datpooled_fd: pooled data for first-dose observations
#'   - datpooled_efd: pooled data for repeated dosing
#'   - datpooled_all: pooled data combining first-dose and repeated-dose
#'     observations
#'
#' @author Zhonghui Huang
#'
#' @examples
#' dat <- processData(Bolus_1CPT)$dat
#' get_pooled_data(dat, dose_type = "combined_doses")
#'
#' @seealso \link{pooled_control}, \link{trimmed_geom_mean}
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




