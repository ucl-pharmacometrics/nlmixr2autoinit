#' Run graphical analysis of PK Parameters
#'
#' Performs graphical estimation of pharmacokinetic parameters based on the input data and route of administration,
#' returning detailed results in a list format.
#'
#' @param dat A data frame containing pharmacokinetic data.
#' @param route Administration route, must be one of \code{"bolus"}, \code{"infusion"}, or \code{"oral"}.
#' @param dose_type Type of analysis to perform. One of:
#'   \itemize{
#'     \item{\code{"first_dose"}}: First-dose analysis (default)
#'     \item{\code{"repeated_doses"}}: Analysis of doses beyond the first (e.g., steady-state)
#'     \item{\code{"combined_doses"}}: Analysis combining first and repeated doses
#'   }
#' @param pooled Optional externally provided pooled data. If \code{NULL}, the data will be pooled internally.
#' @param ... Additional arguments passed to \code{bin.time} or to the graphical calculation functions.
#'
#' @return A list containing graphical estimation results, including clearance, volume of distribution, terminal slope,
#' extrapolated concentration, and for oral data, absorption rate constant.
#'
#' @seealso \code{\link{graphcal_iv}}, \code{\link{graphcal_oral}}, \code{\link{get_pooled_data}}
#'
#' @examples
#' # Example 1 (iv case)
#' dat <- Bolus_1CPT
#' dat <- processData(dat)$dat
#' run_graphcal(dat, route="bolus")
#'
#' # Example 2 (oral case)
#' dat <- Oral_1CPT
#' dat <- processData(dat)$dat
#' run_graphcal(dat, route="oral")
#'
#' # Example 3 (infusion case).
#' Approximate calculation. only use when the infusion duration is very short
#'
#' dat <- Infusion_1CPT
#' dat <- processData(dat)$dat
#' run_graphcal(dat, route="infusion")
#'
#' @export

run_graphcal <- function(dat,
                         route,
                         dose_type = c("first_dose", "repeated_doses", "combined_doses"),
                         pooled = NULL,
                         pooled_ctrl=pooled_control(),
                         ...) {

  dots <- list(...)

  graph_args <- dots[names(dots) %in% names(formals(graphcal_iv))]

  if (route %in% c("bolus","infusion")){
  graph.fd.output <-list(
    cl = NA,
    vd = NA,
    slope = NA,
    C0 = NA,
    method = NA,
    slopefit = NA,
    time.spent = NA
  )
  }

  if (route %in% c("oral")){
    graph.fd.output <-list(
      ka= NA,
      cl = NA,
      vd = NA,
      slope = NA,
      C0 = NA,
      method = NA,
      slopefit = NA,
      time.spent = NA
    )
  }
    # Generate pooled data if not already supplied
  if (is.null(pooled)) {
    pooled <- get_pooled_data(dat = dat,
                              dose_type = dose_type,
                              pooled_ctrl=pooled_ctrl)
  }

  # Initialize default output
  if (route == "bolus" || route == "infusion") {

    if (!is.null(pooled$datpooled_fd) &&
        "binned.df" %in% names(pooled$datpooled_fd)) {

      graph.fd.output <- do.call(graphcal_iv, c(
        list(dat = pooled$datpooled_fd$binned.df, dose = 1),
        graph_args
      ))
    }
  }

  if (route == "oral") {

    if (!is.null(pooled$datpooled_fd) &&
        "binned.df" %in% names(pooled$datpooled_fd)) {

      graph.fd.output <- do.call(graphcal_oral, c(
        list(dat = pooled$datpooled_fd$binned.df, dose = 1),
        graph_args
      ))
    }
  }

  return(graph.fd.output)
}




#' Graphical Calculation of Clearance and Volume of Distribution (IV Route)
#'
#' Performs graphical calculation of pharmacokinetic parameters including clearance (CL),
#' volume of distribution (Vd), terminal slope (lambda_z), and extrapolated concentration at time zero (C0)
#' based on intravenous (IV) pharmacokinetic data.
#'
#' @param dat A data frame containing at least two columns: \code{TIME} (time after dosing) and \code{DV} (measured drug concentration).
#' @param dose Administered dose amount (default is \code{1}).
#' @param ... Additional arguments passed to \code{\link{force_find_lambdaz}}, such as \code{nlastpoints}, \code{adj_r_squared_threshold}, and \code{tolerance}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{cl}}{Estimated clearance (CL)}
#'   \item{\code{vd}}{Estimated volume of distribution (Vd)}
#'   \item{\code{slope}}{Estimated negative terminal phase slope (lambda_z)}
#'   \item{\code{C0}}{Extrapolated concentration at time zero}
#'   \item{\code{method}}{Slope estimation method used ("find_best_lambdaz" or "fallback_regression")}
#'   \item{\code{time.spent}}{Elapsed computation time (in seconds)}
#' }
#' If sufficient valid data points are not available, returned values will be \code{NA}.
#'
#' @details
#' Estimation of the terminal elimination slope (lambda_z) is performed using the \code{\link{force_find_lambdaz}} function,
#' which applies a two-step strategy combining optimal phase selection and fallback regression if necessary.
#'
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
#'                   DV = c(12, 8, 5, 3, 2, 1.5, 1))
#' graphcal_iv(dat, dose = 100)
#'
#' @seealso \code{\link{force_find_lambdaz}}
#' @export
#'
graphcal_iv <- function(dat,
                        dose = 1,
                        ...) {
  start.time <- Sys.time()

  dots <- list(...)
  slope_args <-
    dots[names(dots) %in% names(formals(find_best_lambdaz))]

  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  cl <- NA
  vd <- NA
  slope <- NA
  C0 <- NA
  method <- NA
  slopefit <- NULL

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <- dat[max_index:nrow(dat),]

  if (nrow(temp1) >= 2) {
    result <- force_find_lambdaz(time = temp1$TIME,
                                 conc = temp1$DV, ...)
    kel <- result$lambdaz
    method <- result$method
    slopefit <- result$slopefit

    if (!is.na(kel)) {
      slope <- -kel
      if (!is.na(result$intercept)) {
        C0 <- exp(result$intercept)
      }
      if (!is.na(C0) && C0 > 0) {
        vd <- dose / C0
        cl <- kel * vd
      }
    }
  }

  end.time <- Sys.time()
  time.spent <-
    round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

  return(
    list(
      cl = cl,
      vd = vd,
      slope = slope,
      C0 = C0,
      method = method,
      slopefit = slopefit,
      time.spent = time.spent
    )
  )
}


#' Graphical calculation of PK parameters for oral administration
#'
#' Calculates key pharmacokinetic parameters from oral pharmacokinetic data using graphical methods,
#' including absorption rate constant (ka), elimination rate constant (kel), terminal slope,
#' extrapolated concentration (C0), apparent volume of distribution (Vd/F), and clearance (Cl/F).
#'
#' @param dat A data frame containing at least two columns: \code{TIME} (time after dosing) and \code{DV} (measured drug concentration).
#' @param dose Administered dose amount (default is \code{1}).
#' @param ... Additional arguments passed to \code{\link{find_best_lambdaz}}, such as \code{nlastpoints}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{ka}}{Estimated absorption rate constant (1/h)}
#'   \item{\code{kel}}{Estimated elimination rate constant (1/h)}
#'   \item{\code{slope}}{Negative terminal phase slope (lambda_z)}
#'   \item{\code{C0}}{Extrapolated concentration at the start of elimination phase}
#'   \item{\code{cl}}{Estimated clearance normalized by bioavailability (Cl/F)}
#'   \item{\code{vd}}{Estimated apparent volume of distribution normalized by bioavailability (Vd/F)}
#'   \item{\code{method}}{Method used for terminal phase slope estimation ("find_best_lambdaz" or "fallback_regression")}
#'   \item{\code{time.spent}}{Elapsed computation time in seconds}
#' }
#'
#' @details
#' Estimation of the terminal elimination slope (lambda_z) is performed using the \code{\link{force_find_lambdaz}} function,
#' which applies a two-step strategy combining optimal phase selection and fallback regression if necessary.
#'
#' \deqn{
#' Vd/F = \frac{Dose \times ka}{C_0 \times (ka - kel)}
#' }
#'
#' \deqn{
#' Cl/F = kel \times Vd/F
#' }
#'
#' where \code{ka} is estimated separately using residual analysis of the absorption phase.
#'
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
#'                 DV = c(1, 2, 5, 3, 2, 1.5, 1))
#' graphcal_oral(dat, dose = 100,route="oral")
#'
#' @seealso \code{\link{find_best_lambdaz}}
#' @export


graphcal_oral <- function(dat,
                          dose = 1,
                          ...) {
  start.time <- Sys.time()

  dots <- list(...)
  slope_args <-
    dots[names(dots) %in% names(formals(find_best_lambdaz))]

  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"

  ka <- NA
  kel <- NA
  slope <- NA
  C0 <- NA
  vd <- NA
  cl <- NA
  method <- NA
  slopefit <- NULL

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <-
    dat[(max_index + 1):nrow(dat), ]  # After Tmax (exclude Tmax)

  if (nrow(temp1) >= 2) {
    result <-
      force_find_lambdaz(time = temp1$TIME, conc = temp1$DV, ...)

    kel <- result$lambdaz
    method <- result$method
    slopefit <- result$slopefit

    if (!is.na(kel)) {
      slope <- -kel
      if (!is.na(result$intercept)) {
        C0exp <- exp(result$intercept)
      }
    }

    if (!is.na(kel) && !is.na(C0exp)) {
      Cmax_point <- which(dat$DV == max(dat$DV), arr.ind = TRUE)
      absorb_phase <- dat[1:Cmax_point, ]
      absorb_phase <- absorb_phase[absorb_phase$TIME > 0, ]

      if (nrow(absorb_phase) > 0) {
        absorb_phase$IVconc <- C0exp * exp(-kel * absorb_phase$TIME)
        absorb_phase$residuals <-
          absorb_phase$IVconc - absorb_phase$DV

        absorb_phase <- rbind(data.frame(
          TIME = 0,
          DV = 0,
          IVconc = C0exp,
          residuals = C0exp
        ),
        absorb_phase)

        absorb_phase <- absorb_phase[absorb_phase$residuals > 0, ]

        if (nrow(absorb_phase) >= 2) {
          abslinear <- lm(log(absorb_phase$residuals) ~ absorb_phase$TIME)
          slope_ka <- summary(abslinear)[[4]][[2]]
          if (!is.na(slope_ka)) {
            ka <- -slope_ka
          }
        } else {
          ka <- NA
        }
      }
    }

    if (!is.na(ka) && !is.na(kel) && !is.na(C0exp)) {
      if (abs(ka - kel) > .Machine$double.eps) {
        vd <- (dose * ka) / (C0exp * (ka - kel))
        cl <- kel * vd
      }
    }
  }

  end.time <- Sys.time()
  time.spent <-
    round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

  return(
    list(
      ka = ka,
      kel = kel,
      slope = slope,
      C0exp = C0exp,
      cl = cl,
      vd = vd,
      method = method,
      slopefit = slopefit,
      time.spent = time.spent
    )
  )
}
