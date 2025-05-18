#' Control settings for non-compartmental analysis (NCA)
#'
#' @param trapezoidal_rule Integer. 1 = linear, 2 = linear-up log-down (default = 1).
#' @param half_life Numeric or NA. Optional fixed half-life override.
#' @param nlastpoints Integer. Number of terminal points for half-life regression.
#'
#' @return A list with NCA control parameters.
#' @export
#' @examples
#' nca_control()
nca_control <-
  function(trapezoidal.rule = c("linear_up_log_down", "linear"),
           duration=NULL,
           nlastpoints = 3) {
    list(trapezoidal.rule = trapezoidal.rule ,
         duration = duration,
         nlastpoints = nlastpoints)
  }

#' Perform Non-Compartmental Analysis (NCA) on Pooled Pharmacokinetic Data
#' Performs non-compartmental analysis (NCA) using pooled pharmacokinetic data
#' for various dosing scenarios, including first dose, repeated doses, or a combined
#' analysis of both.
#'
#' @param dat A data frame containing pharmacokinetic data. Required columns depend on
#'            the specified \code{dose_type} (see Details).
#' @param dose_type Type of analysis to perform. One of:
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
#'   Any result not relevant to the \code{dose_type} will return \code{NA}.
#'
#' @details
#' This function supports flexible NCA workflows where pooled data can either be computed
#' internally using \code{\link{get_pooled_data}}, or supplied externally (e.g., reused across calls).
#' Required columns in \code{dat} are automatically checked within \code{get_pooled_data}.
#'
#' For infusion analyses (\code{route = "infusion"}), the representative infusion duration is determined
#' based on observation records (\code{EVID == 0}). Each sampling time point is associated with a corresponding
#' infusion duration, recorded in the \code{durationobs} column. This value reflects the infusion responsible
#' for the observed concentration at that time. When estimating a pooled duration (e.g., for lambda-z estimation),
#' the most commonly used value in \code{durationobs} is selected to ensure consistency with the analyzed profile.
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: Pre-pooling externally using get_pooled_data()
#' dat <- Bolus_1CPT[Bolus_1CPT$SD == 1, ]
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' pooled <- get_pooled_data(fdat, dose_type = "first_dose")
#' run_pooled_nca(dat, dose_type = "first_dose", route=froute, pooled = pooled)
#'
#' # Example 2: Automatic pooling within run_pooled_nca()
#' dat <- Bolus_1CPT
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' pooled <- get_pooled_data(dat, dose_type = "first_dose")
#' run_pooled_nca(fdat, dose_type = "combined_doses", route=froute)
#'
#'
#' # Example 4: Infusion case
#' dat <- Infusion_1CPT
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' fdosetype<-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Type"]
#' pooled <- get_pooled_data(fdat, dose_type = fdosetype)
#' run_pooled_nca(fdat, dose_type = fdosetype, route=froute)
#'
#' }
#'
#' @seealso \code{\link{get_pooled_data}}, \code{\link{bin.time}}, \code{\link{getnca}}
#' @importFrom magrittr %>%
#' @export

run_pooled_nca <- function(dat,
                           route = c("bolus", "oral", "infusion"),
                           dose_type = c("first_dose", "repeated_doses", "combined_doses"),
                           pooled = NULL,
                           pooled_ctrl=pooled_control(),
                           nca_ctrl=nca_control()) {

  `%>%` <- magrittr::`%>%`

  # --- Check for infusion route and extract most common II from EVID == 1 ---
  if (!is.null(route) && route == "infusion") {
    if (!"durationobs" %in% names(dat) || !"EVID" %in% names(dat)) {
      stop(
        "To determine infusion duration, `dat` must contain 'durationobs' and 'EVID' columns.",
        call. = FALSE
      )
    }

    # Only use EVID == 0 (observation rows)
    obs_rows <- dat[dat$EVID == 0 & !is.na(dat$durationobs),]
    if (nrow(obs_rows) == 0) {
      stop("No valid EVID == 0 rows with non-missing 'durationobs' found.",
           call. = FALSE)
    }

    # Get the mode (most frequent durationobs value)
    most_common_duration <-
      as.numeric(names(which.max(table(
        obs_rows$durationobs
      ))))

    message(crayon::black(
      sprintf(
        "Duration for NCA was set to %g, based on the most commonly observed duration among sampling records.",
        most_common_duration
      )
    ))

    # Assign only if user hasn't manually specified it
    if (is.null(nca_ctrl$duration)) {
      nca_ctrl$duration <- most_common_duration
    }
  }

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
      aumc_0_t = NA,
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
      aumc_0_t = NA,
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
      aumc_0_t = NA,
      aumc_0_inf = NA,
      used_points = NA,
      adj.r.squared = NA,
      messages = NA,
      time.spent = NA
    )
  )

  # Generate pooled data if not already supplied
  if (is.null(pooled)) {
    pooled <- get_pooled_data(dat = dat,
                              dose_type = dose_type,
                              pooled_ctrl=pooled_ctrl)
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
      nca_ctrl
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
      nca_ctrl
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
      nca_ctrl
    ))
  }

  return(results)
}
