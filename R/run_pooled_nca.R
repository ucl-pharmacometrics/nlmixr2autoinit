#' Control options for non-compartmental analysis
#'
#' Control options for non-compartmental analysis (NCA)
#'
#' @param trapezoidal.rule Character. Method for trapezoidal AUC integration:
#' \itemize{
#'   \item \code{"linear"} - Linear trapezoidal rule (default)
#'   \item \code{"linear_up_log_down"} - Linear-up / log-down rule
#' }
#'
#' @param duration Numeric. Optional. Duration of the observation window (same
#'   units as time). Used to restrict the integration or define the evaluation
#'   range.
#'
#' @param nlastpoints Integer. Number of terminal points for half-life
#'   regression (default = 3).
#'
#' @param slope.method Character. Method for estimating the terminal slope
#'   (lambdaz):
#' \itemize{
#'   \item \code{"bestfitforce"} - Force estimation using decreasing number of
#'     terminal points if best-fit fails (default)
#'   \item \code{"bestfit"} - Use automated best-fit selection based on adjusted
#'     R-squared
#' }
#'
#' @return A list with NCA control parameters.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' nca_control()
#' @export

nca_control <-
  function(trapezoidal.rule = c("linear_up_log_down", "linear"),
           duration=NULL,
           nlastpoints = 3,
           slope.method= "bestfitforce") {
    list(trapezoidal.rule = trapezoidal.rule ,
         duration = duration,
         nlastpoints = nlastpoints,
         slope.method=slope.method)
  }


#' Performs non-compartmental analysis on pooled data
#'
#' Implements pooled concentration–time profiling followed by non-compartmental
#' analysis (NCA) to derive pharmacokinetic parameters across single-dose, multiple-dose,
#' or combined dosing scenarios under bolus, oral, or infusion routes.
#'
##' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#' @param dose_type Classified as first_dose, repeated_doses, or combined_doses
#'   based on whether observed concentrations occur following the first
#'   administration, during repeated dosing, or across both contexts.
#' @param route Administration route: "bolus", "oral", or "infusion".
#' @param pooled Optional pre-pooled data returned by `get_pooled_data`.
#' @param pooled_ctrl Optional list of control parameters used by `get_pooled_data()`
#'   for pooling observations. Defaults to output from `pooled_control()`.
#' @param nca_ctrl List of options created by `nca_control` for NCA settings.
#'
#' @return A list containing NCA results according to the selected dose_type.
#'
#' @details
#' The function first pools individual subject data into representative
#' concentration–time profiles using `get_pooled_data` based on the settings in
#' `pooled_ctrl`. The pooled profiles are then passed to `getnca`, which computes
#' non-compartmental parameters using rules specified in `nca_ctrl`.
#'
#' @examples
#' out   <- processData(Bolus_1CPT)
#' dat   <- out$dat
#' route <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#'
#' run_pooled_nca(
#'   dat       = dat,
#'   dose_type = "first_dose",
#'   route     = route
#' )$nca.fd.results
#'
#' @seealso \code{\link{get_pooled_data}}, \code{\link{bin.time}}, \code{\link{getnca}}
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
