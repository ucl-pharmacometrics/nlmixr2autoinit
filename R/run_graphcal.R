#' Run graphical analysis of PK Parameters
#'
#' Performs graphical estimation of pharmacokinetic parameters based on the input data and route of administration,
#' returning detailed results in a list format.
#'
#' @param dat A data frame containing pharmacokinetic data.
#' @param route Administration route, must be one of \code{"bolus"}, \code{"infusion"}, or \code{"oral"}.
#' @param data_type Type of dataset to use for binning (default: \code{"first_dose"}).
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
                         data_type = "first_dose",
                         pooled = NULL,
                         ...) {

  dots <- list(...)
  bin_args <- dots[names(dots) %in% names(formals(bin.time))]
  graph_args <- dots[names(dots) %in% names(formals(graphcal_iv))]

  if (is.null(pooled)) {
    pooled <- do.call(get_pooled_data,
                      c(list(dat = dat, data_type = data_type),
                        bin_args))
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
#' volume of distribution (Vd), slope of terminal phase, and extrapolated concentration at time zero (C0)
#' based on intravenous (IV) pharmacokinetic data.
#'
#' @param dat A data frame containing at least two columns: \code{TIME} (time after dosing) and \code{DV} (measured drug concentration).
#' @param ... Additional arguments passed to \code{\link{find_best_lambdaz}}, such as
#'   \code{nlastpoints}, \code{adj_r_squared_threshold}, and \code{tolerance}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{cl}}{Estimated clearance (CL)}
#'   \item{\code{vd}}{Estimated volume of distribution (Vd)}
#'   \item{\code{slope}}{Estimated negative terminal phase slope (lambda_z)}
#'   \item{\code{C0}}{Extrapolated concentration at time zero}
#' }
#' If sufficient valid data points are not available, returned values will be \code{NA}.
#'
#' @details
#' The function applies a two-step fallback strategy to estimate terminal phase parameters:
#' \enumerate{
#'   \item It first attempts to identify the optimal elimination phase using \code{\link{find_best_lambdaz}}.
#'         If a valid slope (lambda_z) is found, the intercept from the regression is used to calculate C0, and Vd and CL are derived accordingly.
#'   \item If no valid lambda_z is found (i.e., \code{NA}), and if there are at least 2 data points available after \code{Tmax},
#'         the function falls back to simple linear regressions, starting from \code{n-1} points down to 2 points.
#'         It sequentially attempts to fit log-linear regressions; the first successful fit with a negative slope is accepted
#'         for recalculating lambda_z, C0, Vd, and CL.
#'   \item If fewer than 2 points are available, or no valid regression is found, the function returns \code{NA} for all parameters.
#' }
#'
#' @examples
#' dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10),
#'                DV = c(12, 8, 5, 3, 2, 1.5, 1))

#' graphcal_iv(dat,dose=100)
#'
#' @seealso \code{\link{find_best_lambdaz}}
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

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <- dat[max_index:nrow(dat),]

  if (nrow(temp1) >= 2) {
    # Step 1: Try find_best_lambdaz
    result <- do.call(find_best_lambdaz, c(
      list(
        time = temp1$TIME,
        conc = temp1$DV,
        route = "bolus"
      ),
      slope_args
    ))
    kel <- result$lamdaz

    if (!is.na(kel)) {
      # Step 2: kel found successfully
      method <- "find_best_lambdaz"
      slope <- -kel
      if (!is.null(result$slopefit)) {
        Intercept <- summary(result$slopefit)[[4]][[1]]
        C0 <- exp(Intercept)
      } else {
        C0 <- NA
      }
      if (!is.na(C0) && C0 > 0) {
        vd <- dose / C0
        cl <- kel * vd
      }

    } else {
      # Step 3: Fallback - try from (nrow(temp1)-1) down to 2 points
      fallback_success <- FALSE

      for (k in seq(nrow(temp1) - 1, 2, by = -1)) {
        fallback_points <- tail(temp1, k)
        fit <-
          try(lm(log(fallback_points$DV) ~ fallback_points$TIME),
              silent = TRUE)

        if (inherits(fit, "try-error"))
          next

        coefs <- summary(fit)[[4]]
        slope_val <- coefs[2]
        Intercept <- coefs[1]

        if (!is.na(slope_val) && slope_val < 0) {
          method <- "fallback_regression"
          slope <- slope_val
          kel <- -slope_val
          C0 <- exp(Intercept)

          if (!is.na(C0) && C0 > 0) {
            vd <- dose / C0
            cl <- kel * vd

            fallback_success <- TRUE
            break
          }
        }
      }
    }
  }

  end.time <- Sys.time()
  time.spent <-
    round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

  return(list(
    cl = cl,
    vd = vd,
    slope = slope,
    C0 = C0,
    method = method,
    time.spent = time.spent
  ))
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
#' The function first attempts to identify the elimination phase using \code{\link{find_best_lambdaz}},
#' based on adjusted R-squared criteria. If a valid elimination phase cannot be identified,
#' it falls back to simple linear regressions using decreasing numbers of points (from \code{n-1} down to 2).
#'
#' The extrapolated intercept (C0) from the log-linear regression represents the concentration at the start
#' of the elimination phase. The apparent volume of distribution (Vd/F) and clearance (Cl/F) are estimated
#' under the assumption of complete absorption (F=1) and based on the relationship:
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
#'                  DV = c(1, 2, 5, 3, 2, 1.5, 1))
#' graphcal_oral(dat, dose = 100)
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
  method <- "NA"

  # Identify Tmax
  max_index <- which.max(dat$DV)
  temp1 <-
    dat[(max_index + 1):nrow(dat), ]  # After Tmax (exclude Tmax)

  if (nrow(temp1) >= 2) {
    # Step 1: Try find_best_lambdaz
    result <- do.call(find_best_lambdaz, c(
      list(
        time = temp1$TIME,
        conc = temp1$DV,
        route = "oral"
      ),
      slope_args
    ))
    kel <- result$lamdaz

    if (!is.na(kel)) {
      method <- "find_best_lambdaz"
      slope <- -kel
      if (!is.null(result$slopefit)) {
        Intercept <- summary(result$slopefit)[[4]][[1]]
        C0exp <- exp(Intercept)
      } else {
        C0exp <- NA
      }
    } else {
      # Step 2: Fallback manual regression
      for (k in seq(nrow(temp1) - 1, 2, by = -1)) {
        fallback_points <- tail(temp1, k)
        fit <-
          try(lm(log(fallback_points$DV) ~ fallback_points$TIME), silent = TRUE)
        if (inherits(fit, "try-error"))
          next

        coefs <- summary(fit)[[4]]
        slope_val <- coefs[2]
        Intercept <- coefs[1]

        if (!is.na(slope_val) && slope_val < 0) {
          method <- "fallback_regression"
          slope <- slope_val
          kel <- -slope_val
          C0exp <- exp(Intercept)
          break
        }
      }
    }

    # Step 3: Estimate ka
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

        # Remove negative or zero residuals
        absorb_phase <- absorb_phase[absorb_phase$residuals > 0,]

        # Only proceed if enough points remain
        if (nrow(absorb_phase) >= 2) {
          abslinear <- lm(log(absorb_phase$residuals) ~ absorb_phase$TIME)
          slope_ka <- summary(abslinear)[[4]][[2]]
          if (!is.na(slope_ka)) {
            ka <- -slope_ka
          }
        } else {
          ka <- NA  # Not enough points left
        }
      }
    }

    # Step 4: Estimate Vd/F and Cl/F if ka estimated
    if (!is.na(ka) && !is.na(kel) && !is.na(C0exp)) {
      if (abs(ka - kel) > .Machine$double.eps) {
        vd <- (dose * ka) / (C0exp * (ka - kel))  # Vd/F
        cl <- kel * vd                            # Cl/F
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
      time.spent = time.spent
    )
  )
}
