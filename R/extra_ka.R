#' Calculate the absorption rate constant using the Wagner-Nelson method
#'
#' Calculates absorption rate constant using the Wagner–Nelson method for
#' single-dose extravascular pharmacokinetics.
#'
#' @param dat A data frame containing two columns: 'TIME' for sampling time
#'   points and 'DV' for observed plasma drug concentrations.
#'
#' @param nca.out Optional object containing results from a previous
#'   noncompartmental analysis. It must include 'auc0_inf' for the area under
#'   the concentration-time curve extrapolated to infinity and 'lambdaz' for the
#'   terminal elimination rate constant. If not provided, the function calls
#'   'getnca' internally using the input data.
#'
#' @details
#' The Wagner-Nelson method estimates the fraction of drug absorbed over time
#' based on the principle of mass balance, where the unabsorbed fraction is
#' quantified as the proportion of the administered dose that has not yet
#' entered systemic circulation. A linear regression is applied to the natural
#' logarithm of the unabsorbed fraction versus time, and the negative slope of
#' this regression corresponds to the first-order absorption rate constant 'ka'.
#'
#' Key assumptions:
#' * Single-dose oral or extravascular administration
#' * First-order absorption and first-order elimination
#' * Linear pharmacokinetics with 'ka' greater than 'ke'
#'
#' Computational steps:
#' - AUC is calculated using trapezoidal integration.
#' - The fraction absorbed is calculated from AUC and the terminal elimination
#'   rate constant.
#' - The remaining fraction is transformed using the natural logarithm.
#' - Linear regression of log(remaining fraction) against time yields 'ka'.
#'
#' @return A list containing:
#'   - ka: Estimated absorption rate constant
#'   - dat_out_wanger_nelson: Input data frame augmented with calculated
#'     pharmacokinetic variables including cumulative AUC, fraction absorbed,
#'     and fraction remaining
#'
#' @references
#' Wagner JG and Nelson E (1963). Percent absorbed time plots derived from
#' blood level and/or urinary excretion data. Journal of Pharmaceutical
#' Sciences, 52(6), 610-611.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Simulated one-compartment oral absorption data
#' Dose <- 100
#' Fbio <- 1
#' Vd   <- 70
#' CL   <- 4
#' ka   <- 1.2
#' ke   <- CL / Vd
#' t  <- seq(0.5, 8, by = 0.5)
#' Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
#'
#' dat <- data.frame(TIME = t, DV = Ct)
#'
#' ka_wanger_nelson(dat)
#' @export

ka_wanger_nelson <- function(dat,
                             nca.out = NULL) {
  colnames(dat)[1] <- "TIME"
  colnames(dat)[2] <- "DV"
  x <- dat$TIME
  y <- dat$DV

  if (missing(nca.out)) {
    nca.out <- getnca(
      x = dat$TIME,
      y = dat$DV,
      ss = 0,
      route = "oral"
    )
  }

  auc_ <- nca.out$auc0_inf
  ke <- nca.out$lambdaz
  auc_intervals <- diff(x) * (y[-1] + y[-length(y)]) / 2

  # use a triangle to calculate 0 - first sample
  dat$auc_intervals <- c(dat[1, ]$DV * dat[1, ]$TIME / 2, auc_intervals)
  dat$auc_accumulate <- NA

  for (k in 1:nrow(dat)) {
    if (k == 1) {
      dat[k, ]$auc_accumulate <- dat[1, ]$auc_intervals
    }
    else{
      dat[k, ]$auc_accumulate <-
        sum(dat[dat$TIME <= dat[k, ]$TIME,]$auc_intervals)
    }
  }
  dat$frac_abs <-  (dat$DV + ke * dat$auc_accumulate) / (ke * auc_)

  # max_abs_time <- which(dat$frac_abs == max(dat$frac_abs ), arr.ind = TRUE)
  # dat_absorb_phase <- dat[1:max_abs_time, ]
  # dat_absorb_phase$frac_remained<-1-  dat_absorb_phase$frac_abs
  dat$frac_remained <- 1 -  dat$frac_abs
  dat_absorb_phase <- dat[dat$frac_remained > 0.1, ]
  # dat_absorb_phase<- dat_absorb_phase[dat_absorb_phase$frac_remained>0.05,]
  if (nrow(dat_absorb_phase) < 2) {
    ka = NA
  }

  if (nrow(dat_absorb_phase) > 1) {
    abc2 <-
      lm(log(dat_absorb_phase$frac_remained) ~  dat_absorb_phase$TIME)
    slope_ka <- summary(abc2)[[4]][[2]]
    ka = -slope_ka
  }
  return(list(ka = ka, dat_out_wanger_nelson = dat))
}



#' Estimate absorption rate constant in a one-compartment oral model
#'
#' This estimates the absorption rate constant in a single-dose oral model using
#' first-order pharmacokinetics.
#'
#' @param cl Numeric. Clearance of the drug.
#' @param ke Numeric. Elimination rate constant.
#' @param t Numeric. Time after administration.
#' @param Ct Numeric. Observed plasma concentration at time t.
#' @param Fbio Numeric. Absolute bioavailability fraction. Default is 1.
#' @param Dose Numeric. Administered oral dose.
#'
#' @details
#' The model assumes a one-compartment structure with first-order absorption and
#' first-order elimination.
#'
#' The concentration-time relationship is:
#' \deqn{Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left( e^{-ke \cdot t} - e^{-ka \cdot t} \right)}
#' where the volume of distribution is defined as:
#' \deqn{Vd = \frac{cl}{ke}}
#'
#' ka is estimated using `uniroot()`, which solves for the root of the residual
#' function (predicted Ct - observed Ct) within a bounded interval (ka > ke and ka <= 1000)
#'
#' @return A list containing:
#' \item{ka}{Estimated absorption rate constant.}
#' \item{full_solution}{The full result object returned by the root-finding process.}
#' \item{message}{A character string indicating the status of the estimation or any warnings.}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example from Oral_1CPT dataset (ID = 1, 1st dose, t = 0.5 h)
#' ka_calculation_sd(cl = 3.62, ke = 0.0556, t = 0.5, Ct = 310.6, Dose = 60000)
#'
#' @export

ka_calculation_sd <- function(cl,
                              ke,
                              t,
                              Ct,
                              Fbio = 1,
                              Dose) {
  Vd <- cl / ke
  ka.equation_sd <- function(ka) {
    predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
      (exp(-ke * t) - exp(-ka * t))
    return(predicted_Ct - Ct)
  }
  # Use the uniroot function to numerically solve for ka
  # The reasonable range for ka is provided as [lower, upper] for root finding
  # To avoid flip-flop, lower boundry of ka was set as higher than ke
  ka_lower <- ke + 1e-4
  ka_upper <- 1000

  # Theoretical upper limit of Ct
  max_pred_Ct <- (Fbio * Dose / Vd) * exp(-ke * t)

  # Return early if Ct is too high
  if (!is.finite(max_pred_Ct) || Ct > max_pred_Ct) {
    msg <-
      "Observed Ct exceeds model-predicted maximum concentration; ka estimation skipped."
    return(list(
      ka = NA,
      full_solution = NULL,
      message = msg
    ))
  }

  # Check boundary function values
  f_lower <-
    tryCatch(
      ka.equation_sd(ka_lower),
      error = function(e)
        NA
    )
  f_upper <-
    tryCatch(
      ka.equation_sd(ka_upper),
      error = function(e)
        NA
    )

  if (!is.finite(f_lower) ||
      !is.finite(f_upper) || f_lower * f_upper > 0) {
    msg <-
      "Function values at uniroot boundaries are invalid or not of opposite sign."
    return(list(
      ka = NA,
      full_solution = NULL,
      message = msg
    ))
  }

  # Try solving via uniroot safely
  solution <- tryCatch({
    stats::uniroot(ka.equation_sd, lower = ka_lower, upper = ka_upper)
  }, error = function(e) {
    return(structure(NULL, class = "try-error", message = e$message))
  })

  if (inherits(solution, "try-error")) {
    return(list(
      ka = NA,
      full_solution = NULL,
      message = paste("uniroot error:", attr(solution, "message"))
    ))
  }

  return(list(
    ka = solution$root,
    full_solution = solution,
    message = "complete"
  ))
}



#' Calculate absorption rate constant (ka) in a multiple-dose one-compartment model
#'
#' This estimates the absorption rate constant in a multiple-dose oral model using
#' first-order pharmacokinetics.
#'
#' @param cl Numeric. Clearance of the drug (in L/hr).
#' @param ke Numeric. Elimination rate constant (in 1/hr).
#' @param t Numeric. Time after the last dose (in hours) at which the concentration is measured.
#' @param Ct Numeric. Observed concentration of the drug at time \code{t} (in mg/L).
#' @param Fbio Numeric. Bioavailability fraction (default = 1, meaning 100% bioavailability).
#' @param Dose Numeric. Administered dose of the drug (in mg).
#' @param tau Numeric. Dosing interval (in hours) between successive doses.
#'
#' @return A list containing the following components:
#' \item{ka}{The calculated absorption rate constant.}
#' \item{full_solution}{The full solution object returned by the \code{uniroot()}
#' function, which includes additional details about the root-finding process.}
#'
#' @details
#' The value of ka is obtained numerically using the uniroot unction by solving
#' the following equation:
#'
#' \deqn{Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left(
#' \frac{e^{-ke \cdot t}}{1 - e^{-ke \cdot \tau}} - \frac{e^{-ka \cdot t}}{1 -
#' e^{-ka \cdot \tau}} \right)}
#'
#' ka is estimated using `uniroot()`, which solves for the root of the residual
#' function (predicted Ct - observed Ct) within a bounded interval (ka > ke and ka <= 1000)
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example from Oral_1CPT dataset (ID = 1, 5th dose, t = 2 h)
#' ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 852, Dose = 60000, tau = 24)
#'
#' @export
ka_calculation_md <- function(cl,
                              ke,
                              t,
                              Ct,
                              Fbio = 1,
                              Dose,
                              tau) {

  Vd <- cl / ke
  ka.equation_md <- function(ka) {
    predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
      ((exp(-ke * t) / (1 - exp(-ke * tau))) - (exp(-ka * t) / (1 - exp(-ka * tau))))

    return(predicted_Ct - Ct)
  }

  # Use the uniroot function to numerically solve for ka
  # The reasonable range for ka is provided as [lower, upper] for root finding
  # To avoid flip-flop, lower boundry of ka was set as higher than ke
  ka_lower <- ke + 1e-4
  ka_upper <- 1000

  # Estimate theoretical max Ct (ka to infinity)
  max_pred_Ct <-
    (Fbio * Dose) / Vd * (exp(-ke * t) / (1 - exp(-ke * tau)))

  if (!is.finite(max_pred_Ct) || Ct > max_pred_Ct) {
    msg <-
      "Observed Ct exceeds model-predicted maximum concentration (multiple-dose); ka estimation skipped."
    return(list(
      ka = NA,
      full_solution = NULL,
      message = msg
    ))
  }

  # Check signs at boundaries
  f_lower <-
    tryCatch(
      ka.equation_md(ka_lower),
      error = function(e)
        NA
    )
  f_upper <-
    tryCatch(
      ka.equation_md(ka_upper),
      error = function(e)
        NA
    )

  if (!is.finite(f_lower) ||
      !is.finite(f_upper) || f_lower * f_upper > 0) {
    msg <-
      "Function values at uniroot boundaries are invalid or not of opposite sign (multiple-dose)."
    return(list(
      ka = NA,
      full_solution = NULL,
      message = msg
    ))
  }

  # Try solving with uniroot
  solution <- tryCatch({
    stats::uniroot(ka.equation_md, lower = ka_lower, upper = ka_upper)
  }, error = function(e) {
    return(structure(NULL, class = "try-error", message = e$message))
  })

  if (inherits(solution, "try-error")) {
    return(list(
      ka = NA,
      full_solution = NULL,
      message = paste(
        "uniroot error (multiple-dose):",
        attr(solution, "message")
      )
    ))
  }

  # Success
  return(list(
    ka = solution$root,
    full_solution = solution,
    message = "complete"
  ))
}

#' Estimate the absorption rate constant using pointwise methods
#'
#' It implements pointwise estimation of absorption rate constants for single-dose
#' and multiple-dose pharmacokinetic models.
#'
#' @param df A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @param cl A numeric value for drug clearance. It is assumed constant across subjects unless pre-specified per subject.
#' @param ke A numeric value for the elimination rate constant. This is assumed known or estimated from the terminal phase.
#' @param Fbio A numeric value for bioavailability (F). Default is 1.
#'
#' @details
#' For each subject, the time of maximum observed concentration (Tmax) is identified as the time corresponding to the highest DV.
#' Only records with TIME less than or equal to Tmax are retained, representing the absorption phase.
#'
#' Two scenario-specific calculations are implemented: single-dose and multiple-dose at steady state.
#' \itemize{
#'   \item For single-dose data (dose_number == 1 and SteadyState == FALSE), the function uses `ka_calculation_sd()`,
#'   which applies a one-compartment oral absorption model under first-order absorption and elimination.
#'   \item For steady-state multiple-dose data (dose_number > 1 and SteadyState == TRUE), the function uses `ka_calculation_md()`,
#'   which accounts for accumulation using the dosing interval (dose_interval).
#' }
#'
#' This function does not perform model fitting. The median is recommended for use in pharmacokinetic modeling.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{ka_calc_median}{Median ka value across all valid observations.}
#'   \item{ka_calc_dat_sd}{Data frame containing absorption-phase single-dose data, with estimated ka and diagnostic messages.}
#'   \item{ka_calc_dat_md}{Data frame containing absorption-phase steady-state multiple-dose data, with ka estimates.}
#'   \item{ka_calc_dat}{Combined data frame (single-dose and multiple-dose) containing all ka estimates.}
#' }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Single-dose
#' df <- Oral_1CPT[Oral_1CPT$SD == 1, ]
#' df <- processData(df)$dat
#' df <- is_ss(df)
#' run_ka_solution(df = df, cl = 4, ke = 4/70, Fbio = 1)$ka_calc_median
#'
#' # Mixed doses
#' dat <- Oral_1CPT
#' df_ss <- processData(dat)$dat
#' df_ss <- is_ss(df_ss)
#' run_ka_solution(df = df_ss, cl = 4, ke = 4/70, Fbio = 1)$ka_calc_median
#'
#' @export

run_ka_solution <- function(df,
                            cl,
                            ke,
                            Fbio = 1) {
  # Step 1: Find Tmax for each individual
  tmax_df <- df %>%
    dplyr::group_by(ID) %>%
    dplyr::summarize(Tmax = TIME[which.max(DV)])

  # Step 2: Join Tmax with the original data and filter rows where TIME <= Tmax
  data_before_tmax <- df %>%
    dplyr::left_join(tmax_df, by = "ID") %>%
    dplyr::filter(TIME <= Tmax) %>%
    dplyr::mutate(
      cl = cl,
      ke = ke,
      Fbio = Fbio,
      ka_calc = NA_real_
    )

  data_before_tmax_sd <- NULL
  data_before_tmax_md <- NULL

  if (any(
    data_before_tmax$dose_number == 1 &
    data_before_tmax$EVID == 0 &
    data_before_tmax$SSflag == 0
  )) {
    data_before_tmax_sd <- data_before_tmax %>%
      dplyr::filter(dose_number == 1, EVID == 0, SteadyState == FALSE)

    # Use mapply to compute ka and diagnostic message for each row in single-dose data
    ka_sd_results <- mapply(
      function(cl, ke, t, Ct, Fbio, Dose) {
        res <- try(ka_calculation_sd(
          cl = cl,
          ke = ke,
          t = t,
          Ct = Ct,
          Fbio = Fbio,
          Dose = Dose
        ),
        silent = TRUE)

        # Return NA and message if error occurred
        if (inherits(res, "try-error")) {
          return(list(ka = NA, message = "try-error"))
        }

        # Return estimated ka and status message
        list(ka = res$ka, message = res$message)
      },
      cl    = data_before_tmax_sd$cl,
      ke    = data_before_tmax_sd$ke,
      t     = data_before_tmax_sd$TIME,
      # time after single dose
      Ct    = data_before_tmax_sd$DV,
      # observed concentration
      Fbio  = data_before_tmax_sd$Fbio,
      Dose  = data_before_tmax_sd$dose,
      SIMPLIFY = FALSE
    )

    # Extract ka and message columns
    data_before_tmax_sd$ka_calc    <-
      sapply(ka_sd_results, `[[`, "ka")
    data_before_tmax_sd$ka_message <-
      sapply(ka_sd_results, `[[`, "message")

    # Optional: convert ka column to numeric explicitly, suppressing warnings
    data_before_tmax_sd$ka_calcv <-
      suppressWarnings(suppressMessages(as.numeric(data_before_tmax_sd$ka_calc)))
  }

  if (any(
    data_before_tmax$dose_number > 1 &
    data_before_tmax$EVID == 0 &
    data_before_tmax$SteadyState
  )) {
    data_before_tmax_md <- data_before_tmax %>%
      dplyr::filter(dose_number > 1, EVID == 0, SteadyState)

    ka_md_results <- mapply(
      function(cl, ke, t, Ct, Fbio, Dose, tau) {
        # Try to estimate ka using ka_calculation_md()
        res <- try(ka_calculation_md(
          cl = cl,
          ke = ke,
          t = t,
          Ct = Ct,
          Fbio = Fbio,
          Dose = Dose,
          tau = tau
        ),
        silent = TRUE)

        # If error occurs, return NA and a message
        if (inherits(res, "try-error")) {
          return(list(ka = NA, message = "try-error"))
        }

        # Return both ka and message from the result
        list(ka = res$ka, message = res$message)
      },
      cl    = data_before_tmax_md$cl,
      ke    = data_before_tmax_md$ke,
      t     = data_before_tmax_md$tad,
      # time after last dose
      Ct    = data_before_tmax_md$DV,
      # observed concentration
      Fbio  = data_before_tmax_md$Fbio,
      Dose  = data_before_tmax_md$dose,
      tau   = data_before_tmax_md$dose_interval,
      # dosing interval
      SIMPLIFY = FALSE
    )

    # Extract ka and message columns from result list
    data_before_tmax_md$ka_calc    <-
      sapply(ka_md_results, `[[`, "ka")
    data_before_tmax_md$ka_message <-
      sapply(ka_md_results, `[[`, "message")

    # Optional: convert ka column to numeric explicitly, suppressing warnings
    data_before_tmax_md$ka_calcv <-
      suppressWarnings(suppressMessages(as.numeric(data_before_tmax_md$ka_calc)))
  }

  data_before_tmax <-
    rbind(data_before_tmax_sd, data_before_tmax_md)

  ka_calc_median <-
    suppressWarnings(suppressMessages(median(data_before_tmax$ka_calcv, na.rm = TRUE)))

  if (is.null(ka_calc_median)) {
    ka_calc_median <- NA
  }
  # trimmed_mean_ka <-
  #   tryCatch(
  #     trimmed_geom_mean(data_before_tmax$ka_calc, trim = 0.05, na.rm = TRUE),
  #     error = function(e) {
  #       NA
  #     }
  #   )

  return(
    list(
      ka_calc_median = ka_calc_median,
      ka_calc_dat_sd = data_before_tmax_sd,
      ka_calc_dat_md = data_before_tmax_md,
      ka_calc_dat = data_before_tmax
    )
  )
}
