#' Perform non-compartmental analysis (NCA) for pooled pharmacokinetic Data
#'
#' Performs non-compartmental analysis on pooled pharmacokinetic data
#' or different dosing scenarios, including first dose, repeated doses, or
#' combined analysis. It supports both bolus and oral administration routes.
#'
#' @param dat A data frame containing pharmacokinetic data. Required columns depend on
#'            the analysis type (see Details).
#' @param data_type Type of analysis to perform. Must be one of:
#'   \itemize{
#'     \item{\code{"first_dose"}}{First dose analysis (default)}
#'     \item{\code{"repeated_doses"}}{Steady-state/repeated dose analysis}
#'     \item{\code{"combined_doses"}}{Combined first and repeated dose analysis}
#'   }
#' @param nlastpoints Number of points for terminal slope calculation (default: 3)
#' @param trapezoidal.rule Integration method:
#'   \itemize{
#'     \item{1}{Linear trapezoidal rule (default)}
#'     \item{2}{Log-linear trapezoidal rule}
#'   }
#' @param nbins Number of time bins for data pooling (default: 8)
#' @param route Administration route:
#'   \itemize{
#'     \item{\code{"bolus"}}{Bolus administration (default)}
#'     \item{\code{"oral"}}{Oral administration}
#'   }
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item{nca.fd.results}{NCA results for first dose analysis (data frame)}
#'     \item{nca.efd.results}{NCA results for repeated doses analysis (data frame)}
#'     \item{nca.all.results}{NCA results for combined analysis (data frame)}
#'     \item{datpooled_fd}{Binned data for first dose analysis}
#'     \item{datpooled_efd}{Binned data for repeated doses analysis}
#'     \item{datpooled_all}{Binned data for combined analysis}
#'     \item{time.spent}{Calculation time (difftime object)}
#'   }
#'   Unused analysis slots will contain \code{NA}.
#' @details
#' Required columns in \code{dat} depend on the \code{data_type}:
#'   \itemize{
#'     \item{\code{"first_dose"}:}{\code{dose_number}, \code{iiobs}, \code{ID}, \code{TIME}, \code{DV}}
#'     \item{\code{"repeated_doses"}:}{\code{EVID}, \code{AMT}, \code{tad}, \code{dose_number},
#'                                    \code{iiobs}, \code{ID}, \code{TIME}, \code{DV}, \code{resetflag}}
#'     \item{\code{"combined_doses"}:}{\code{EVID}, \code{tad}, \code{dose_number}, \code{ID},
#'                                    \code{TIME}, \code{DV}, \code{resetflag}}
#'   }
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # First dose analysis for bolus data
#' processed_data <- processData(Bolus_1CPT)$dat
#' run_pooled_nca(dat = processed_data, data_type = "combined_doses")
#'
#' # Combined analysis for oral data
#' data(Oral_1CPT)
#' processed_oral <- processData(Oral_1CPT)$dat
#' results_oral <- run_pooled_nca(processed_oral,
#'                           data_type = "combined_doses",
#'                           route = "oral")
#' }
#'
#' @seealso \code{\link{bin.time}}, \code{\link{getnca}}
#' @export

run_pooled_nca <- function(dat,
                           data_type = "first_dose",
                           ...) {
  start.time <- Sys.time()

  dots <- list(...)

  # Split user-supplied arguments into those meant for bin.time and getnca
  bin_args <- dots[names(dots) %in% names(formals(bin.time))]
  nca_args <- dots[names(dots) %in% names(formals(getnca))]

  # Initialize results structure with NA placeholders
  results <- list(
    nca.fd.results = NA,
    nca.efd.results = NA,
    nca.all.results = NA,
    datpooled_fd = NA,
    datpooled_efd = NA,
    datpooled_all = NA,
    time.spent  = NA
  )

  # Validate analysis type parameter
  valid_types <- c("first_dose", "repeated_doses", "combined_doses")
  if (!data_type %in% valid_types) {
    stop(
      "Invalid data_type: '",
      data_type,
      "'. Valid options are: ",
      paste(valid_types, collapse = ", ")
    )
  }

  # Define required columns based on analysis type
  required_cols <- switch(
    data_type,
    "first_dose" = c("dose_number", "iiobs", "ID", "TIME", "DV"),
    "repeated_doses" = c(
      "EVID",
      "AMT",
      "tad",
      "dose_number",
      "iiobs",
      "ID",
      "TIME",
      "DV",
      "resetflag"
    ),
    "combined" = c("EVID", "tad", "dose_number", "ID", "TIME", "DV", "resetflag")
  )

  # Check for missing columns
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns for ",
      data_type,
      " analysis: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Main analysis
  if (data_type == "first_dose" || data_type == "combined_doses") {
    # Filter and sort first dose observations
    fd_data <- dat[dat$dose_number == 1 & dat$iiobs == 0, ]
    fd_data <- fd_data[order(fd_data$ID, fd_data$TIME), ]

    pooled_data <-
      do.call(bin.time, c(list(dat = fd_data), bin_args))

    results$nca.fd.results <- do.call(getnca, c(
      list(
        x = pooled_data$binned.df$Time,
        y = pooled_data$binned.df$Conc,
        ss = 0
      ),
      nca_args
    ))

    results$datpooled_fd <- pooled_data
  }

  if (data_type == "repeated_doses" ||
      data_type == "combined_doses") {
    # Calculate most common dosing interval
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

    # Prepare repeated-doses data
    efd_data <- dat %>%
      dplyr::filter(
        (dose_number != 1 & EVID == 0 & tad <= most_common_ii * 1.2) |
          (dose_number != 1 & EVID == 1) |
          (dose_number == 1 & iiobs > 0)
      ) %>%
      dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

    pooled_data <-
      do.call(bin.time, c(list(dat = efd_data), bin_args))

    results$nca.efd.results <- do.call(getnca, c(
      list(
        x = pooled_data$binned.df$Time,
        y = pooled_data$binned.df$Conc,
        ss = 1
      ),
      nca_args
    ))

    results$datpooled_efd <- pooled_data

    if (data_type == "combined_doses") {
      # Prepare combined data
      all_data <- dat %>%
        dplyr::filter((EVID == 0 &
                         tad <= most_common_ii * 1.2) |
                        EVID == 1) %>%
        dplyr::arrange(ID, resetflag, TIME, dplyr::desc(AMT))

      pooled_data <-
        do.call(bin.time, c(list(dat = all_data), bin_args))

      results$nca.all.results <- do.call(getnca, c(
        list(
          x = pooled_data$binned.df$Time,
          y = pooled_data$binned.df$Conc,
          ss = 1
        ),
        nca_args
      ))

    }
  }
  end.time <- Sys.time()
  time.spent <- round(difftime(end.time, start.time), 4)
  results$time.spent <- time.spent

  return(results)
}
