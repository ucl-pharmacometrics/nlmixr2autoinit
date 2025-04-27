#' Perform Non-Compartmental Analysis (NCA) on Pooled Pharmacokinetic Data
#' Performs non-compartmental analysis (NCA) using pooled pharmacokinetic data
#' for various dosing scenarios, including first dose, repeated doses, or a combined
#' analysis of both.
#'
#' @param dat A data frame containing pharmacokinetic data. Required columns depend on
#'            the specified \code{data_type} (see Details).
#' @param data_type Type of analysis to perform. One of:
#'   \itemize{
#'     \item{\code{"first_dose"}}: First-dose analysis (default)
#'     \item{\code{"repeated_doses"}}: Analysis of doses beyond the first (e.g., steady-state)
#'     \item{\code{"combined_doses"}}: Analysis combining first and repeated doses
#'   }
#' @param pooled (Optional) A precomputed pooled data list as returned by \code{\link{get_pooled_data}}.
#'               This must be a list that may contain one or more of the following named elements:
#'   \describe{
#'     \item{\code{datpooled_fd}}{Binned data for first-dose analysis (output of \code{bin.time})}
#'     \item{\code{datpooled_efd}}{Binned data for repeated-dose analysis}
#'     \item{\code{datpooled_all}}{Binned data combining first and repeated doses}
#'   }
#'               If not supplied, the function will compute pooled data automatically.
#' @param ... Additional arguments passed to \code{\link{bin.time}} and \code{\link{getnca}}.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item{\code{nca.fd.results}}: NCA results for first-dose data (if applicable)
#'     \item{\code{nca.efd.results}}: NCA results for repeated doses (if applicable)
#'     \item{\code{nca.all.results}}: NCA results for combined first and repeated doses (if applicable)
#'   }
#'   Any result not relevant to the \code{data_type} will return \code{NA}.
#'
#' @details
#' This function supports flexible NCA workflows where pooled data can either be computed
#' internally using \code{\link{get_pooled_data}}, or supplied externally (e.g., reused across calls).
#' Required columns in \code{dat} are automatically checked within \code{get_pooled_data}.
#'
#' @examples
#' \dontrun{
#' # First-dose NCA analysis (pre-filtered dataset)
#' dat <- Bolus_1CPT[Bolus_1CPT$SD == 1, ]
#' dat <- processData(dat)$dat
#'
#' # Option 1: Pre-pooling externally using get_pooled_data()
#' pooled <- get_pooled_data(dat, data_type = "first_dose", binsize = 1)
#' result <- run_pooled_nca(dat, data_type = "first_dose", pooled = pooled)
#'
#' # Option 2: Automatic pooling within run_pooled_nca()
#' result <- run_pooled_nca(dat, data_type = "first_dose")
#'
#' # Combined analysis: first-dose + repeated-dose profiles
#' dat <- Bolus_1CPT
#' dat <- processData(dat)$dat
#' result <- run_pooled_nca(dat, data_type = "combined_doses")
#' }
#'
#' @seealso \code{\link{get_pooled_data}}, \code{\link{bin.time}}, \code{\link{getnca}}
#' @importFrom magrittr %>%
#' @export

run_pooled_nca <- function(dat,
                           data_type = "first_dose",
                           pooled = NULL,
                           ...) {
  dots <- list(...)

  # Split user-supplied arguments into those meant for bin.time and getnca
  bin_args <- dots[names(dots) %in% names(formals(bin.time))]
  nca_args <- dots[names(dots) %in% names(formals(getnca))]

  # Initialize results structure with NA placeholders
  results <- list(
    nca.fd.results = list(
      clobs = NA,
      vzobs = NA,
      half_life = NA,
      auct = NA,
      auc0_inf = NA,
      C_last = NA,
      lambdaz = NA,
      aumc_0_inf = NA,
      used_points = NA,
      adj.r.squared = NA,
      messages = NA,
      time.spent = NA
    ),
    nca.efd.results = list(
      clobs = NA,
      vzobs = NA,
      half_life = NA,
      auct = NA,
      auc0_inf = NA,
      C_last = NA,
      lambdaz = NA,
      aumc_0_inf = NA,
      used_points = NA,
      adj.r.squared = NA,
      messages = NA,
      time.spent = NA
    ),
    nca.all.results = list(
      clobs = NA,
      vzobs = NA,
      half_life = NA,
      auct = NA,
      auc0_inf = NA,
      C_last = NA,
      lambdaz = NA,
      aumc_0_inf = NA,
      used_points = NA,
      adj.r.squared = NA,
      messages = NA,
      time.spent = NA
    )
  )

  # Generate pooled data if not already supplied
  if (is.null(pooled)) {
    pooled <- do.call(get_pooled_data,
                      c(list(
                        dat = dat, data_type = data_type
                      ),
                      bin_args))
  }

  # Perform NCA analysis where data is available
  if (!is.null(pooled$datpooled_fd) &&
      "binned.df" %in% names(pooled$datpooled_fd)) {
    results$nca.fd.results <- do.call(getnca, c(
      list(
        x = pooled$datpooled_fd$binned.df$Time,
        y = pooled$datpooled_fd$binned.df$Conc,
        ss = 0
      ),
      nca_args
    ))
  }

  if (!is.null(pooled$datpooled_efd) &&
      "binned.df" %in% names(pooled$datpooled_efd)) {
    results$nca.efd.results <- do.call(getnca, c(
      list(
        x = pooled$datpooled_efd$binned.df$Time,
        y = pooled$datpooled_efd$binned.df$Conc,
        ss = 1
      ),
      nca_args
    ))
  }

  if (!is.null(pooled$datpooled_all) &&
      "binned.df" %in% names(pooled$datpooled_all)) {
    results$nca.all.results <- do.call(getnca, c(
      list(
        x = pooled$datpooled_all$binned.df$Time,
        y = pooled$datpooled_all$binned.df$Conc,
        ss = 1
      ),
      nca_args
    ))
  }

  return(results)
}
