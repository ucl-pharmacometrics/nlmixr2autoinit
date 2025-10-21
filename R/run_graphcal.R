#' Run graphical analysis of pharmacokinetic parameters
#'
#' Performs graphical estimation of pharmacokinetic parameters based on pooled
#' concentration–time data and the specified route of administration.
#'
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#' @param route Route of administration. Must be one of bolus, oral, or infusion.
#' @param dose_type Specifies the dosing context of the pharmacokinetic
#'   observations. Classified as first_dose, repeated_doses, or combined_doses
#'   based on whether observed concentrations occur following the first
#'   administration, during repeated dosing, or across both contexts.
#' @param pooled Optional pooled dataset. If NULL, pooling is performed internally.
#' @param pooled_ctrl Control settings created by `pooled_control()` for time binning and pooling.
#' @param ... Additional arguments passed to graphical calculation functions.
#'
#' @details
#' The function pools individual profiles using `get_pooled_data()` when needed,
#' and then applies route-specific graphical methods (`graphcal_iv` or `graphcal_oral`)
#' to estimate parameters such as clearance, volume of distribution, terminal slope,
#' and absorption rate constant (for oral data).
#'
#' @return A list containing graphical estimates of key pharmacokinetic parameters.
#'
#' @seealso \code{\link{graphcal_iv}}, \code{\link{graphcal_oral}}, \code{\link{get_pooled_data}}

#' @author Zhonghui Huang
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
#' # Approximate calculation. only use when the infusion duration is very short
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
    C0exp = NA,
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
      C0exp = NA,
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

#' Graphical calculation of clearance and volume of distribution (IV route)
#'
#' Estimates clearance (CL), volume of distribution (Vd), terminal slope
#' (lambdaz), and extrapolated concentration at time zero (C0exp) from
#' intravenous pharmacokinetic data using graphical methods.
#'
#' @param dat A data frame containing TIME (time after dosing) and DV
#'   (observed concentration).
#' @param dose Administered dose amount. Defaults to 1.
#' @param ... Additional arguments passed to `force_find_lambdaz()`.
#'
#' @details
#' Terminal slope (lambdaz) is estimated using `force_find_lambdaz()`, which
#' applies an automated phase selection strategy with fallback regression when
#' required.
#'
#' @return A list containing graphical estimates of CL, Vd, lambda_z, and C0exp.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
#'                   DV = c(12, 8, 5, 3, 2, 1.5, 1))
#' graphcal_iv(dat, dose = 100)
#'
#' @seealso \code{\link{force_find_lambdaz}}
#' @export

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
  C0exp <- NA
  method <- NA
  slopefit <- NULL

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <- dat[max_index:nrow(dat), ]

  if (nrow(temp1) >= 2) {
    result <- force_find_lambdaz(time = temp1$TIME,
                                 conc = temp1$DV, ...)
    kel <- result$lambdaz
    method <- result$method
    slopefit <- result$slopefit

    if (!is.na(kel)) {
      slope <- -kel
      if (!is.na(result$intercept)) {
        C0exp <- exp(result$intercept)
      }
      if (!is.na(C0exp) && C0exp > 0) {
        vd <- dose / C0exp
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
      C0exp = C0exp,
      method = method,
      slopefit = slopefit,
      time.spent = time.spent
    )
  )
}

#' Graphical calculation of pharmacokinetic parameters for oral administration
#'
#' Estimates key pharmacokinetic parameters from oral concentration–time data using
#' graphical methods, including absorption rate constant (ka), elimination rate
#' constant (kel), terminal slope, extrapolated concentration (C0exp), apparent
#' volume of distribution (Vd/F), and clearance (Cl/F).
#'
#' @param dat A data frame containing TIME (time after dosing) and DV (observed
#'   concentration).
#' @param dose Administered dose amount. Defaults to 1.
#' @param ... Additional arguments passed to `find_best_lambdaz()`.
#'
#' @details
#' The terminal slope (lambdaz) is estimated using `force_find_lambdaz()`. The
#' apparent volume of distribution and clearance are computed using the
#' following relationships:
#' \deqn{Vd/F = \frac{Dose \times ka}{C_0 \times (ka - kel)}}
#' \deqn{Cl/F = kel \times Vd/F}
#' where \code{ka} is estimated from the absorption phase.
#'
#' @return A list containing graphical estimates of ka, kel, lambda_z, C0exp,
#' Vd/F, and Cl/F.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
#'                   DV = c(1, 2, 5, 3, 2, 1.5, 1))
#' graphcal_oral(dat, dose = 100, route = "oral")
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
  C0exp <- NA
  vd <- NA
  cl <- NA
  method <- NA
  slopefit <- NULL

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <-
    dat[(max_index + 1):nrow(dat),]  # After Tmax (exclude Tmax)

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
      absorb_phase <- dat[1:Cmax_point,]
      absorb_phase <- absorb_phase[absorb_phase$TIME > 0,]

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

        absorb_phase <- absorb_phase[absorb_phase$residuals > 0,]

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
