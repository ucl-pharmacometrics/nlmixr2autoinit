#' Calculate the Absorption Rate Constant (ka) using the Wanger-Nelson Method
#'
#' Computes the absorption rate constant (ka) using the Wanger-Nelson method.
#' The function makes use of non-compartmental analysis (NCA) function part to calculate the area under the curve (AUC) and the elimination rate constant (ke).
#' Linear regression on the logarithm of the remaining fraction of drug is used to estimate the absorption rate constant (ka).
#'
#' @param dat A data frame with at least two columns: `DV` (drug concentration)
#' and `TIME` (time points of measurement).
#'
#' @details
#' The function calculates the absorption rate constant (ka) as follows:
#' \itemize{
#'   \item The function uses the `nca.iv.normalised` method to calculate the area under the curve (AUC) and elimination rate constant (ke).
#'   \item AUC is calculated using trapezoidal integration for the observed data points.
#'   \item The fraction of drug absorbed is calculated based on ke and the cumulative AUC.
#'   \item The logarithm of the remaining fraction of the drug is then modeled using linear regression to determine the absorption rate constant (ka).
#' }
#'
#' This method assumes that the drug absorption follows first-order kinetics.
#'
#' @return A numeric value representing the estimated absorption rate constant (ka).
#'
#' @examples
#'
#' # data frame from nlmixr2data
#' dat<-Oral_1CPT[Oral_1CPT$ID==1 &Oral_1CPT$SD==1&Oral_1CPT$EVID==0,]
#' dat<-data.frame(TIME=dat$TIME,DV=dat$DV)
#' ka_wanger_nelson(dat=dat)$ka
#' ka_wanger_nelson(dat=dat)$dat_out_wanger_nelson
#'
#' @export
#'
ka_wanger_nelson<-function(dat,
                           nca.out=NULL){

  colnames(dat)[1]<-"TIME"
  colnames(dat)[2]<-"DV"
  x<-dat$TIME
  y<-dat$DV

  if (missing(nca.out)) {
    nca.out <- getnca(x = dat$TIME, y = dat$DV,ss = 0)
  }

  auc_ <- nca.out$auc0_inf
  ke<- nca.out$lambdaz
  auc_intervals<- diff(x) * (y[-1] + y[-length(y)])/2

  # use a triangle to calculate 0 - first sample
  dat$auc_intervals<-c( dat[1,]$DV * dat[1,]$TIME/2, auc_intervals)
  dat$auc_accumulate<-NA

  for (k in 1:nrow(dat)){
    if (k==1){
      dat[k,]$auc_accumulate<-dat[1,]$auc_intervals
    }
    else{
      dat[k,]$auc_accumulate<-sum(dat[dat$TIME<=dat[k,]$TIME, ]$auc_intervals)
    }
  }
  dat$frac_abs<-  (dat$DV+ ke*dat$auc_accumulate)/ (ke*auc_)

  # max_abs_time <- which(dat$frac_abs == max(dat$frac_abs ), arr.ind = T)
  # dat_absorb_phase <- dat[1:max_abs_time, ]
  # dat_absorb_phase$frac_remained<-1-  dat_absorb_phase$frac_abs
  dat$frac_remained<-1-  dat$frac_abs
  dat_absorb_phase<- dat[dat$frac_remained>0.1,]
  # dat_absorb_phase<- dat_absorb_phase[dat_absorb_phase$frac_remained>0.05,]
  if (nrow(dat_absorb_phase)<2){
    ka=NA
  }

  if (nrow(dat_absorb_phase)>1){
  abc2 <- lm(log(dat_absorb_phase$frac_remained) ~  dat_absorb_phase$TIME)
  slope_ka<-summary(abc2)[[4]][[2]]
  ka= -slope_ka
  }
  return(list(ka=ka,dat_out_wanger_nelson=dat))
}



#' Calculate absorption rate constant (ka) in a single-dose one-compartment model
#'
#' Calculates the absorption rate constant (\code{ka}) for a drug administered orally
#' using a one-compartment model with first-order absorption and elimination kinetics.
#' The calculation is based on known values of clearance (\code{cl}), elimination rate constant
#' (\code{ke}), time (\code{t}), observed concentration (\code{Ct}), bioavailability (\code{Fbio}),
#' and dose (\code{Dose}).
#'
#' @param cl Numeric. Clearance of the drug (in L/hr).
#' @param ke Numeric. Elimination rate constant (in 1/hr).
#' @param t Numeric. Time after drug administration (in hours) at which the concentration is measured.
#' @param Ct Numeric. Observed concentration of the drug at time \code{t} (in mg/L).
#' @param Fbio Numeric. Bioavailability fraction (default = 1, meaning 100% bioavailability).
#' @param Dose Numeric. Administered dose of the drug (in mg).
#'
#' @return A list containing the following components:
#' \item{ka}{The calculated absorption rate constant.}
#' \item{full_solution}{The full solution object returned by the \code{uniroot()} function, which includes additional details about the root-finding process.}
#' \item{message}{A character string providing status information about the estimation, including warnings if the observed concentration exceeds model-predicted limits or if numerical errors occurred.}
#'
#' @details
#' This function uses a one-compartment model with first-order absorption and elimination to
#' estimate the absorption rate constant (\code{ka}). The function solves the following equation
#' for \code{ka} numerically using \code{uniroot()}:
#'
#' \deqn{Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left( e^{-ke \cdot t} - e^{-ka \cdot t} \right)}
#'
#' The \code{uniroot()} function is used to find the value of \code{ka} that makes the difference
#' between the predicted and observed concentrations equal to zero. The reasonable range for
#' \code{ka} is set as [0.01, 1000] for root-finding.
#'
#' Additional safeguards include ensuring that \code{ka > ke} to avoid flip-flop kinetics, and rejecting
#' values of \code{Ct} that exceed the model-predicted maximum concentration. In such cases, the function
#' returns \code{NA} for \code{ka} and a descriptive \code{message}.
#'
#' @examples
#' # Example usage:
#' # Data from first point in simulated dataset Oral_1CPT
#' ka_result <- ka_calculation_sd(cl = 3.62, ke = 0.0556, t = 0.5, Ct = 310, Dose = 60000)
#' print(ka_result$ka)
#'
#' @export
ka_calculation_sd <- function(cl,     # Clearance of the drug (L/hr)
                              ke,     # Elimination rate constant (1/hr)
                              t,      # Time (hr) after drug administration at which concentration is measured
                              Ct,     # Observed concentration of the drug at time 't' (mg/L)
                              Fbio = 1, # Bioavailability fraction, default is 1 (100% bioavailability)
                              Dose) {
  Vd <- cl / ke
  # Define the equation to solve for the absorption rate constant (ka)
  ka.equation_sd <- function(ka) {
    # Predicted concentration using the one-compartment model with first-order absorption and elimination
    predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
      (exp(-ke * t) - exp(-ka * t))

    # Return the difference between predicted concentration and observed concentration
    return(predicted_Ct - Ct)
  }

  # Use the uniroot function to numerically solve for ka (absorption rate constant)
  # The reasonable range for ka is provided as [lower, upper] for root finding
  # To avoid flip-flop, lower boundry of ka was set as higher than ke
  ka_lower <- ke + 1e-4
  ka_upper <- 1000

  # Theoretical upper limit of Ct as ka → ∞
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
    uniroot(ka.equation_sd, lower = ka_lower, upper = ka_upper)
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

  # Successful result
  return(list(
    ka = solution$root,
    full_solution = solution,
    message = "complete"
  ))
}



#' Calculate absorption rate constant (ka) in a multiple-dose one-compartment model
#'
#' Estimates the absorption rate constant (\code{ka}) for a drug administered orally
#' in a multiple-dose regimen using a one-compartment model with first-order absorption
#' and elimination kinetics. The calculation is based on known values of clearance (\code{cl}),
#' elimination rate constant (\code{ke}), time since last dose (\code{t}), observed concentration
#' (\code{Ct}), bioavailability (\code{Fbio}), dose (\code{Dose}), dosing interval (\code{tau}),
#' and the number of doses administered (\code{n}).
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
#' \item{full_solution}{The full solution object returned by the \code{uniroot()} function, which includes additional details about the root-finding process.}
#'
#' @details
#' This function uses a multiple-dose one-compartment model with first-order absorption and
#' elimination to estimate the absorption rate constant (\code{ka}). The function solves the
#' following equation for \code{ka} numerically using \code{uniroot()}:
#'
#' \deqn{Ct = \frac{Fbio \cdot Dose \cdot ka}{Vd \cdot (ka - ke)} \left( \frac{e^{-ke \cdot t}}{1 - e^{-ke \cdot \tau}} - \frac{e^{-ka \cdot t}}{1 - e^{-ka \cdot \tau}} \right)}
#'
#' The \code{uniroot()} function is used to find the value of \code{ka} that makes the difference
#' between the predicted and observed concentrations equal to zero. The reasonable range for
#' \code{ka} is set as [0.01, 1000] for root-finding.
#'
#' @examples
#' # Example usage:
#' ka_result <- ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 852, Dose = 60000, tau = 24)
#' ka_result <- ka_calculation_md(cl = 4, ke = 0.057, t = 2, Ct = 189.6, Dose = 10000, tau = 24)
#' print(ka_result$ka)
#'
#' @export
ka_calculation_md <- function(cl,      # Clearance of the drug (L/hr)
                              ke,      # Elimination rate constant (1/hr)
                              t,       # Time (hr) after last dose at which concentration is measured
                              Ct,      # Observed concentration of the drug at time 't' (mg/L)
                              Fbio = 1, # Bioavailability fraction, default is 1 (100% bioavailability)
                              Dose,    # Administered dose of the drug (mg)
                              tau) {
  # Dosing interval (hr) between doses

  Vd <- cl / ke  # Calculate the volume of distribution

  # Define the equation to solve for the absorption rate constant (ka)
  ka.equation_md <- function(ka) {
    # Predicted concentration using the multiple-dose formula
    predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
      ((exp(-ke * t) / (1 - exp(-ke * tau))) - (exp(-ka * t) / (1 - exp(-ka * tau))))

    # predicted_Ct <- (Fbio * Dose * ka) / (Vd * (ka - ke)) *
    #   ((1 / (1 - exp(-ke * tau))) - (1 / (1 - exp(-ka * tau))))

    # Return the difference between predicted concentration and observed concentration
    return(predicted_Ct - Ct)
  }

  # Use the uniroot function to numerically solve for ka (absorption rate constant)
  # The reasonable range for ka is provided as [lower, upper] for root finding
  # To avoid flip-flop, lower boundry of ka was set as higher than ke
  ka_lower <- ke + 1e-4
  ka_upper <- 1000


  # Estimate theoretical max Ct (ka → ∞), i.e., first term dominates
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
    uniroot(ka.equation_md, lower = ka_lower, upper = ka_upper)
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

#' Calculate absorption rate constant (ka) for oral administration Data
#'
#' Calculate the absorption rate constant (ka) using sampling points from the absorption phase (t<Tmax). The single-point method is applied, with each plasma concentration point used to calculate a corresponding ka value
#'
#' @param df A data frame containing pharmacokinetic data with columns: `ID` (subject identifier), `TIME` (time point), `DV` (drug concentration), and `DOSE` (dose amount). The data frame is typically filtered to include data only from the first dosing interval (e.g., single dose and evid = 0).
#' @param cl A numeric value representing the clearance (CL) of the drug.
#' @param ke A numeric value representing the elimination rate constant (ke).
#' @param Fbio A numeric value representing the fraction of the drug absorbed (bioavailability, F).
#'
#' @details
#' The function first calculates the time to maximum concentration (Tmax) for each individual in the dataset. It then filters the data to include only the time points up to Tmax for each individual. The absorption rate constant (ka) is calculated for each time point using the `ka_calculation` function.
#'
#' @return A list containing:
#' \item{ka_calc_median}{The median ka value across all individuals.}
#' \item{data_before_tmax}{The original data frame filtered for times before Tmax, with an added `ka_calc` column for individual ka values.}
#'
#' @examples
#' # Example usage:
#'  df<-Oral_1CPT[Oral_1CPT$SD==1,]
#'  df<-processData(df)$dat
#'  result <- run_ka_solution(df = df, cl = 4, ke = 4/70, Fbio = 1)
#'  ka_median <- result[[1]]
#'  data_with_ka <- result[[2]]
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

  if (any(data_before_tmax$dose_number == 1 &
          data_before_tmax$EVID == 0 &
          data_before_tmax$SSflag == 0)) {

    data_before_tmax_sd <- data_before_tmax %>%
      dplyr::filter(dose_number == 1, EVID == 0, SteadyState==F)

    # Use mapply to compute ka and diagnostic message for each row in single-dose data
    ka_sd_results <- mapply(
      function(cl, ke, t, Ct, Fbio, Dose) {
        # Attempt to estimate ka using ka_calculation_sd()
        res <- try(ka_calculation_sd(
          cl = cl,
          ke = ke,
          t = t,
          Ct = Ct,
          Fbio = Fbio,
          Dose = Dose
        ), silent = TRUE)

        # Return NA and message if error occurred
        if (inherits(res, "try-error")) {
          return(list(ka = NA, message = "try-error"))
        }

        # Return estimated ka and status message
        list(ka = res$ka, message = res$message)
      },
      cl    = data_before_tmax_sd$cl,
      ke    = data_before_tmax_sd$ke,
      t     = data_before_tmax_sd$TIME,   # time after single dose
      Ct    = data_before_tmax_sd$DV,     # observed concentration
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
      suppressWarnings(suppressMessages(as.numeric(data_before_tmax_sd$ka_calc)
))
  }

  if (any(data_before_tmax$dose_number > 1 &
          data_before_tmax$EVID == 0 &
          data_before_tmax$SteadyState)) {

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
        ), silent = TRUE)

        # If error occurs, return NA and a message
        if (inherits(res, "try-error")) {
          return(list(ka = NA, message = "try-error"))
        }

        # Return both ka and message from the result
        list(ka = res$ka, message = res$message)
      },
      cl    = data_before_tmax_md$cl,
      ke    = data_before_tmax_md$ke,
      t     = data_before_tmax_md$tad,        # time after last dose
      Ct    = data_before_tmax_md$DV,         # observed concentration
      Fbio  = data_before_tmax_md$Fbio,
      Dose  = data_before_tmax_md$dose,
      tau   = data_before_tmax_md$dose_interval,  # dosing interval
      SIMPLIFY = FALSE
    )

    # Extract ka and message columns from result list
    data_before_tmax_md$ka_calc    <- sapply(ka_md_results, `[[`, "ka")
    data_before_tmax_md$ka_message <- sapply(ka_md_results, `[[`, "message")

    # Optional: convert ka column to numeric explicitly, suppressing warnings
    data_before_tmax_md$ka_calcv <-
      suppressWarnings(suppressMessages(as.numeric(data_before_tmax_md$ka_calc)))
  }

  data_before_tmax <- rbind(data_before_tmax_sd, data_before_tmax_md)

  ka_calc_median <-
    suppressWarnings(suppressMessages(median(data_before_tmax$ka_calcv, na.rm = T)))

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
