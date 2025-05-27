#' Parameter sweeping for a one-compartment nonlinear (Michaelis-Menten) model
#'
#' Performs sensitivity analysis via parameter sweeping in a one-compartment pharmacokinetic model
#' with Michaelis-Menten elimination. By systematically varying user-defined or automatically
#' generated values for key pharmacokinetic parameters (e.g., Vmax, Km, Vd, and Ka), the function
#' evaluates model sensitivity and fit across a defined grid of parameter combinations.
#'
#' @param dat A data frame containing the pharmacokinetic dataset.
#' @param sim_vmax A list for Vmax simulation. Default: \code{list(mode = "auto", values = NULL)}.
#'   \describe{
#'     \item{mode}{`"manual"` or `"auto"`.}
#'     \item{values}{Vector of candidate Vmax values if `"manual"`.}
#'     \item{est.cl}{Estimated clearance, required if `"auto"`.}
#'   }
#' @param sim_km A list for Km simulation. Default: \code{list(mode = "auto", values = NULL)}.
#'   \describe{
#'     \item{mode}{`"manual"` or `"auto"`.}
#'     \item{values}{Vector of candidate Km values if `"manual"`.}
#'   }
#' @param sim_vd A list for volume of distribution simulation. Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of Vd values. Required.}
#'   }
#' @param sim_ka A list for absorption rate constant (Ka). Used only if \code{route = "oral"}. Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of Ka values. Required for oral.}
#'   }
#' @param route A character string indicating administration route. One of `"iv"` or `"oral"`. Default is `"iv"`.
#'
#' @return A data frame with simulated parameter combinations and evaluation metrics:
#'   \describe{
#'     \item{Vmax, Km, Vd, Ka}{The parameter values used.}
#'     \item{APE, MAE, MAPE, RMSE, rRMSE}{Model fit evaluation metrics.}
#'     \item{Cumulative.Time.Sec}{Total elapsed time (in seconds) for each simulation.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example 1: IV route
#' out <- sim_sens_1cmpt_mm(
#'   dat = Bolus_1CPTMM,
#'   sim_vmax = list(mode = "auto",est.cl=4),
#'   sim_km   = list(mode = "auto"),
#'   sim_vd   = list(mode = "manual", values = c(70)),
#'   sim_ka   = list(mode = "manual", values = NA),  # NA for IV route
#'   route = "iv"
#' )
#' head(out)
#'
#' # Example 2: Oral route
#' #out <- sim_sens_1cmpt_mm(
#' #  dat = Oral_1CPTMM,
#' #  sim_vmax = list(mode = "auto",est.cl=1),
#' #  sim_km   = list(mode = "auto"),
#' #  sim_vd   = list(mode = "manual", values = c(48)),
#' #  sim_ka   = list(mode = "manual", values = c(1)),  # 1 for oral route
#' #  route = "oral"
#' #)
#' #head(out)
#' }
#'
#' @export
sim_sens_1cmpt_mm <- function(dat,
                              sim_vmax = list(mode = "auto",
                                              values = NULL,
                                              est.cl = NULL),
                              sim_km   = list(mode = "auto",
                                              values = NULL),
                              sim_vd   = list(mode = "manual", values = NULL),
                              sim_ka   = list(mode = "manual", values = NULL),
                              route = c("iv", "oral")) {
  `%>%` <- magrittr::`%>%`
  route <- tryCatch(
    match.arg(route, choices = c("iv", "oral")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `route`: '%s'. Must be 'iv' or 'oral'.",
          as.character(route)
        ),
        call. = FALSE
      )
    }
  )

  sim_vmax$mode <- tryCatch(
    match.arg(sim_vmax$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_vmax$mode`. Must be 'manual' or 'auto'.",
           call. = FALSE)
    }
  )
  sim_km$mode <- tryCatch(
    match.arg(sim_km$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_km$mode`. Must be 'manual' or 'auto'.", call. = FALSE)
    }
  )
  sim_vd$mode <- tryCatch(
    match.arg(sim_vd$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_vd$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  vd_values <- sim_vd$values
  if (is.null(vd_values) || all(is.na(vd_values))) {
    stop("No candidate Vd values available for parameter sweeping.",
         call. = FALSE)
  }

  # Deduplicate similar values if values are close
  vd_values <- na.omit(vd_values) %>%
    sort() %>%
    tibble::tibble(value = .) %>%
    dplyr::mutate(prev = dplyr::lag(value),
                  rel_diff = abs(value - prev) / prev) %>%
    dplyr::filter(is.na(rel_diff) | rel_diff > 0.2) %>%
    dplyr::pull(value)

  if (route == "oral") {
    sim_ka$mode <- tryCatch(
      match.arg(sim_ka$mode, choices = c("manual")),
      error = function(e) {
        stop("Invalid `sim_ka$mode`. Must be 'manual'.", call. = FALSE)
      }
    )
    if (is.null(sim_ka$values) || all(is.na(sim_ka$values))) {
      stop("No candidate Ka values available for parameter sweeping (oral route).",
           call. = FALSE)
    }
  }

  # Handle Vmax/Km grid
  if (sim_vmax$mode == "auto" || sim_km$mode == "auto") {
    if (is.null(sim_vmax$est.cl)) {
      stop("Auto mode for Vmax/Km requires `sim_vmax$est.cl`.",
           call. = FALSE)
    }
    cl_values<- sim_vmax$est.cl
    # Deduplicate similar values if valus are close
    cl_values <- na.omit(cl_values)%>%
      sort() %>%
      tibble::tibble(value = .) %>%
      dplyr::mutate(prev = dplyr::lag(value),
                    rel_diff = abs(value - prev) / prev) %>%
      dplyr::filter(is.na(rel_diff) | rel_diff > 0.2) %>%
      dplyr::pull(value)

    # obtain the observed cmax
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs$DV,
                          list(dat.obs$ID),
                          FUN = max,
                          na.rm = T)
    cmax <- mean(pop.cmax$x, trim = 0.05, na.rm = T)
    km_range <- c(4, 2, 1, 0.5, 0.25, 0.125, 0.1, 0.05) * cmax
    conc_range <- c(0.05, 0.1, 0.25, 0.5, 0.75) * cmax
    combs <- expand.grid(Km = km_range, Conc = conc_range, cl=cl_values )
    combs$Vmax <- (combs$Km + combs$Conc) * combs$cl
    param_grid <- combs %>% dplyr::select(Vmax, Km)
   #  Filter combinations similar to each other
    keep <- rep(TRUE, nrow(param_grid))
    for (i in 1:(nrow(param_grid) - 1)) {
      if (!keep[i]) next
      for (j in (i + 1):nrow(param_grid)) {
        if (!keep[j]) next
        vmax_diff <- abs((param_grid$Vmax[i] - param_grid$Vmax[j]) / param_grid$Vmax[j])
        km_diff <- abs((param_grid$Km[i] - param_grid$Km[j]) / param_grid$Km[j])
        if (vmax_diff <= 0.2 && km_diff <= 0.2) {
          keep[j] <- FALSE
        }
      }
    }
    param_grid <- param_grid[keep, ]

  } else {
    if (is.null(sim_vmax$values) || all(is.na(sim_vmax$values))) {
      stop("No candidate Vmax values available for parameter sweeping.",
           call. = FALSE)
    }
    if (is.null(sim_km$values) || all(is.na(sim_km$values))) {
      stop("No candidate Km values available for parameter sweeping.",
           call. = FALSE)
    }
    param_grid <-
      tidyr::crossing(Vmax = sim_vmax$values, Km = sim_km$values)
  }

  if (route == "oral") {
    param_grid <-
      tidyr::crossing(param_grid, Vd = vd_values, Ka = sim_ka$values)
  } else {
    param_grid <- tidyr::crossing(param_grid, Vd = vd_values) %>%
      dplyr::mutate(Ka = NA_real_)
  }

  start_time <- Sys.time()
  progressr::handlers(
    progressr::handler_progress(format = ":message [:bar] :percent (:current/:total)", width = 80)
  )

  sim.1cmpt.mm.results.all <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(param_grid))
    param_grid %>%
      dplyr::mutate(row = dplyr::row_number()) %>%
      purrr::pmap_dfr(function(Vmax, Km, Vd, Ka, row) {
        p(sprintf("Running simulation: Vmax=%.2f, Km=%.2f", Vmax, Km))
        sim_out <- if (route == "iv") {
          Fit_1cmpt_mm_iv(
            data = dat[dat$EVID != 2, ],
            est.method = "rxSolve",
            input.vmax = Vmax,
            input.km = Km,
            input.vd = Vd,
            input.add = 0
          )
        } else {
          Fit_1cmpt_mm_oral(
            data = dat[dat$EVID != 2, ],
            est.method = "rxSolve",
            input.ka = Ka,
            input.vmax = Vmax,
            input.km = Km,
            input.vd = Vd,
            input.add = 0
          )
        }

        met <-
          metrics.(pred.x = sim_out$cp, obs.y = dat[dat$EVID == 0, ]$DV)
        elapsed <-
          round(difftime(Sys.time(), start_time, units = "secs"), 2)
        rm(sim_out)
        gc()

        tibble::tibble(
          Vmax = Vmax,
          Km = Km,
          Vd = Vd,
          Ka = if (route == "oral")
            Ka
          else
            NA_real_,
          APE = round(met[1], 2),
          MAE = round(met[2], 2),
          MAPE = round(met[3], 2),
          RMSE = round(met[4], 2),
          rRMSE1 = round(met[5], 2),
          rRMSE2 = round(met[6], 2),
          Cumulative.Time.Sec = as.numeric(elapsed)
        )
      })
  })

  return(sim.1cmpt.mm.results.all)
}

#' Parameter sweeping for a two-compartment pharmacokinetic model
#'
#' Performs sensitivity analysis via parameter sweeping in a two-compartment pharmacokinetic model.
#' By systematically varying user-defined or automatically generated values for key pharmacokinetic
#' parameters (e.g., CL, Vc, Vp, Q, and Ka), the function evaluates model sensitivity and fit
#' across a defined grid of parameter combinations.

#' @param dat A data frame containing the pharmacokinetic dataset.
#' @param sim_ka A list with Ka simulation settings (used only if `route = "oral"`). Default:
#'   \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`; only used when `route = "oral"`.}
#'     \item{values}{Vector of Ka values; required for oral.}
#'   }
#' @param sim_cl A list for clearance (CL) simulation. Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of candidate CL values. Required.}
#'   }
#' @param sim_vc A list for central compartment volume (Vc) simulation. Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of candidate Vc values. Required.}
#'   }
#' @param sim_vp A list for peripheral compartment volume (Vp) simulation. Default: \code{list(mode = c("auto", "manual"), values = NULL)}.
#'   \describe{
#'     \item{mode}{`"auto"` or `"manual"`.}
#'     \item{values}{If `"manual"`, a vector of Vp values. If `"auto"`, generated based on Vc/ratio.}
#'   }
#' @param sim_q A list for inter-compartmental clearance (Q). Default:
#'   \code{list(mode = c("auto", "manual"), values = NULL, auto.strategy = c("scaled", "fixed"))}.
#'   \describe{
#'     \item{mode}{`"auto"` or `"manual"`.}
#'     \item{values}{Required if `"manual"`; a vector of Q values.}
#'     \item{auto.strategy}{Only used if `"auto"`; either `"scaled"` (relative to CL) or `"fixed"` (fixed values).}
#'   }
#' @param route A character string indicating administration route. One of `"iv"` (intravenous) or `"oral"` (extravascular). Default is `"iv"`.
#'
#' @return A data frame with simulated parameter combinations and evaluation metrics:
#'   \describe{
#'     \item{Vc, Vp, Q, CL, Ka}{The parameter values used.}
#'     \item{APE, MAE, MAPE, RMSE, rRMSE}{Model fit evaluation metrics.}
#'     \item{Cumulative.Time.Sec}{Total elapsed time (in seconds) for each simulation.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example 1: IV route
#' out <- sim_sens_2cmpt(
#'   dat = Bolus_1CPT,
#'   sim_cl = list(mode = "manual", values = c(4)),
#'   sim_vc = list(mode = "manual", values = c(50, 70)),
#'   sim_vp = list(mode = "auto"),
#'   sim_q  = list(mode = "auto", auto.strategy = "scaled"),
#'   sim_ka = list(mode = "manual", values = NA),  # NA for IV route
#'   route = "iv"
#' )
#' head(out)
#'
#' # Example 2: Oral route
#' out <- sim_sens_2cmpt(
#'   dat = Oral_2CPT,
#'   sim_cl = list(mode = "manual", values = c(4)),
#'   sim_vc = list(mode = "manual", values = c(70)),
#'   sim_vp = list(mode = "manual", values = c(20, 30)),
#'   sim_q  = list(mode = "manual", values = c(10)),
#'   sim_ka = list(mode = "manual", values = c(0.5, 1, 2)),
#'   route = "oral"
#' )
#' head(out)
#' }
#'
#' @export

sim_sens_2cmpt <- function(dat,
                           sim_ka = list(mode = "manual", values = NULL),
                           sim_cl = list(mode = "manual", values = NULL),
                           sim_vc = list(mode = "manual", values = NULL),
                           sim_vp = list(mode = c("auto", "manual"), values = NULL),
                           sim_q  = list(
                             mode = c("auto", "manual"),
                             values = NULL,
                             auto.strategy = c("scaled", "fixed")
                           ),
                           route = c("iv", "oral")) {
  `%>%` <- magrittr::`%>%`

  # Safe check
  # Defensive validation for `route`
  route <- tryCatch(
    match.arg(route, choices = c("iv", "oral")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `route`: '%s'. Must be one of: %s.",
          as.character(route),
          paste(shQuote(c("iv", "oral")), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # Safe check and set sim ka
  sim_ka$mode <- tryCatch(
    match.arg(sim_ka$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_ka$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  if (route == "oral" &&
      (is.null(sim_ka$values) ||
       length(sim_ka$values) == 0 || all(is.na(sim_ka$values)))) {
    stop("No candidate Ka values available for parameter sweeping.",
         call. = FALSE)
  }

  # Defensive checks for sim_cl
  sim_cl$mode <- tryCatch(
    match.arg(sim_cl$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_cl$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  if (is.null(sim_cl$values) ||
      length(sim_cl$values) == 0 || all(is.na(sim_cl$values))) {
    stop("No candidate CL values available for parameter sweeping.",
         call. = FALSE)
  }

  # Safe check and set sim vc
  sim_vc$mode <- tryCatch(
    match.arg(sim_vc$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_vc$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  vc_values <- sim_vc$values
  if (is.null(vc_values) ||
      length(vc_values) == 0 || all(is.na(vc_values))) {
    stop("No candidate Vc values available for parameter sweeping.",
         call. = FALSE)
  }

  # Deduplicate similar values if valus are close
  vc_values <- na.omit(vc_values)%>%
    sort() %>%
    tibble::tibble(value = .) %>%
    dplyr::mutate(prev = dplyr::lag(value),
                  rel_diff = abs(value - prev) / prev) %>%
    dplyr::filter(is.na(rel_diff) | rel_diff > 0.2) %>%
    dplyr::pull(value)

  # Safe check and set sim vp
  sim_vp$mode <- tryCatch(
    match.arg(sim_vp$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_vp$mode`. Must be 'manual' or 'auto'.", call. = FALSE)
    }
  )

  if (sim_vp$mode == "auto") {
    vp_ratios <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
    vp_combs <- expand.grid(vc = vc_values, ratio = vp_ratios)
    vp_combs$vp <- vp_combs$vc / vp_combs$ratio
  } else {
    if (is.null(sim_vp$values) ||
        length(sim_vp$values) == 0 || all(is.na(sim_vp$values))) {
      stop("No candidate Vp values available for parameter sweeping.",
           call. = FALSE)
    }
    vp_combs <- expand.grid(vc = vc_values, vp = sim_vp$values)
  }

  #  Safe check and set sim q
  sim_q$mode <- tryCatch(
    match.arg(sim_q$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_q$mode`. Must be 'manual' or 'auto'.", call. = FALSE)
    }
  )

  if (sim_q$mode == "auto") {
    strategy <- tryCatch(
      match.arg(sim_q$auto.strategy, choices = c("scaled", "fixed")),
      error = function(e) {
        stop("Invalid `sim_q$auto.strategy`. Must be 'scaled' or 'fixed'.",
             call. = FALSE)
      }
    )
    cl_value <- sim_cl$values[1]
    if (strategy == "scaled") {
      q_values <- c(cl_value / 2, cl_value, cl_value * 2)
    } else {
      q_values <- c(1, 10, 100)
    }
  } else {
    if (is.null(sim_q$values) ||
        length(sim_q$values) == 0 || all(is.na(sim_q$values))) {
      stop("No candidate Q values available for parameter sweeping.",
           call. = FALSE)
    }
    q_values <- sim_q$values
  }

  vp_combs <- vp_combs %>%
    dplyr::rename(Vc = vc, Vp = vp)

  if (route == "oral") {
    param_grid <- tidyr::crossing(vp_combs,
                                  Q = q_values,
                                  CL = sim_cl$values,
                                  Ka = sim_ka$values)
  } else {
    param_grid <- tidyr::crossing(vp_combs,
                                  Q = q_values,
                                  CL = sim_cl$values) %>%
      dplyr::mutate(Ka = NA_real_)
  }
  param_grid <- param_grid %>%
    dplyr::select(Vc, Vp, Q, CL, Ka)

  # === Begin Simulation ===
  start_time <- Sys.time()

  # Use a custom handler for pretty output in R console
  progressr::handlers(
    progressr::handler_progress(format = ":message [:bar] :percent (:current/:total)",
                                width = 80)
  )

  sim.2cmpt.results.all <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(param_grid))
    param_grid %>%
      dplyr::mutate(row = dplyr::row_number()) %>%
      purrr::pmap_dfr(function(Vc, Vp, Q, CL, Ka, row) {
        p(sprintf("Running simulation: Vc=%.2f, Vp=%.2f", Vc, Vp))
        sim_out <- if (route == "iv") {
          Fit_2cmpt_iv(
            data = dat[dat$EVID != 2,],
            est.method = "rxSolve",
            input.cl = CL,
            input.vc2cmpt = Vc,
            input.vp2cmpt = Vp,
            input.q2cmpt = Q,
            input.add = 0
          )

        } else {
          Fit_2cmpt_oral(
            data = dat[dat$EVID != 2,],
            est.method = "rxSolve",
            input.ka = Ka,
            input.cl = CL,
            input.vc2cmpt = Vc,
            input.vp2cmpt = Vp,
            input.q2cmpt = Q,
            input.add = 0
          )
        }

        met <-
          metrics.(pred.x = sim_out$cp, obs.y = dat[dat$EVID == 0,]$DV)
        elapsed <-
          round(difftime(Sys.time(), start_time, units = "secs"), 2)

        rm(sim_out)
        gc()

        tibble::tibble(
          Vc = Vc,
          Vp = Vp,
          Q = Q,
          CL = CL,
          Ka = if (route == "oral")
            Ka
          else
            NA,
          APE = round(met[1], 2),
          MAE = round(met[2], 2),
          MAPE = round(met[3], 2),
          RMSE = round(met[4], 2),
          rRMSE1 = round(met[5], 2),
          rRMSE2 = round(met[6], 2),
          Cumulative.Time.Sec = as.numeric(elapsed)
        )
      })
  })

  return(sim.2cmpt.results.all)
}


#' Parameter sweeping for a three-compartment pharmacokinetic model
#'
#' Performs sensitivity analysis via parameter sweeping in a three-compartment pharmacokinetic model.
#' This function systematically varies key pharmacokinetic parameters including central volume (Vc),
#' peripheral volumes (Vp1, Vp2), intercompartmental clearances (Q1, Q2), systemic clearance (CL),
#' and absorption rate constant (Ka, for oral route).
#' The function supports both `"iv"` (intravenous) and `"oral"` administration routes.
#'
#' @param dat A data frame containing the pharmacokinetic dataset.
#' @param sim_vc A list for central compartment volume (Vc) simulation. Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of candidate Vc values. Required.}
#'   }
#' @param sim_vp A list for first peripheral compartment volume (Vp1) simulation. Default: \code{list(mode = c("auto", "manual"), values = NULL)}.
#'   \describe{
#'     \item{mode}{Either `"auto"` or `"manual"`.}
#'     \item{values}{If `"manual"`, a vector of Vp1 values. If `"auto"`, values are computed from Vc and Vc/Vp1 ratios.}
#'   }
#' @param sim_vp2 A list for second peripheral compartment volume (Vp2) simulation. Same structure as \code{sim_vp}.
#' @param sim_q A list for intercompartmental clearance between Vc and Vp1 (Q1). Default: \code{list(mode = c("auto", "manual"), values = NULL, auto.strategy = c("scaled", "fixed"))}.
#'   \describe{
#'     \item{mode}{Either `"auto"` or `"manual"`.}
#'     \item{values}{If `"manual"`, a vector of Q1 values.}
#'     \item{auto.strategy}{If `"auto"`, `"scaled"` (relative to CL) or `"fixed"` (fixed candidates).}
#'   }
#' @param sim_q2 A list for intercompartmental clearance between Vc and Vp2 (Q2). Same structure as \code{sim_q}.
#' @param sim_cl A list for systemic clearance (CL). Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of candidate CL values. Required.}
#'   }
#' @param sim_ka A list with Ka simulation settings (used only when \code{route = "oral"}). Default: \code{list(mode = "manual", values = NULL)}.
#'   \describe{
#'     \item{mode}{Must be `"manual"`.}
#'     \item{values}{Vector of Ka values; required for oral route.}
#'   }
#' @param route A character string indicating administration route. One of `"iv"` or `"oral"`. Default is `"iv"`.
#'
#' @return A data frame with all simulated parameter combinations and associated model performance metrics:
#'   \describe{
#'     \item{Vc, Vp1, Vp2, Q1, Q2, CL, Ka}{Parameter values used in simulation.}
#'     \item{APE, MAE, MAPE, RMSE, rRMSE}{Model fit metrics between predicted and observed data.}
#'     \item{Cumulative.Time.Sec}{Time taken (in seconds) to complete each simulation.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Example 1: IV route (no Ka required)
#' dat <- Bolus_2CPT
#' out <- sim_sens_3cmpt(
#'   dat = dat,
#'   sim_cl = list(mode = "manual", values = c(4)),
#'   sim_vc = list(mode = "manual", values = c(50, 70)),
#'   sim_vp = list(mode = "auto"),
#'   sim_vp2 = list(mode = "auto"),
#'   sim_q = list(mode = "auto", auto.strategy = "scaled"),
#'   sim_q2 = list(mode = "auto", auto.strategy = "scaled"),
#'   route = "iv"
#' )
#' head(out)
#'
#' # Example 2: Oral route
#' dat <- Oral_2CPT
#' out <- sim_sens_3cmpt(
#'   dat = dat,
#'   sim_ka = list(mode = "manual", values = c(1)),
#'   sim_cl = list(mode = "manual", values = c(4)),
#'   sim_vc = list(mode = "manual", values = c(50, 70)),
#'   sim_vp = list(mode = "auto"),
#'   sim_vp2 = list(mode = "auto"),
#'   sim_q = list(mode = "auto", auto.strategy = "scaled"),
#'   sim_q2 = list(mode = "auto", auto.strategy = "scaled"),
#'   route = "oral"
#' )
#' head(out)
#' }
#'
#' @export


sim_sens_3cmpt <- function(dat,
                           sim_vc = list(mode = "manual", values = NULL),
                           sim_vp = list(mode = c("auto", "manual"), values = NULL),
                           sim_vp2 = list(mode = c("auto", "manual"), values = NULL),
                           sim_q = list(
                             mode = c("auto", "manual"),
                             values = NULL,
                             auto.strategy = c("scaled", "fixed")
                           ),
                           sim_q2 = list(
                             mode = c("auto", "manual"),
                             values = NULL,
                             auto.strategy = c("scaled", "fixed")
                           ),
                           sim_cl = list(mode = "manual", values = NULL),
                           sim_ka = list(mode = "manual", values = NULL),
                           route = c("iv", "oral")) {
  `%>%` <- magrittr::`%>%`

  # --- Route check ---
  route <- tryCatch(
    match.arg(route, choices = c("iv", "oral")),
    error = function(e) {
      stop(sprintf(
        "Invalid `route`: '%s'. Must be one of: %s.",
        as.character(route),
        paste(shQuote(c("iv", "oral")), collapse = ", ")
      ),
      call. = FALSE)
    }
  )

  # --- sim_ka ---
  sim_ka$mode <- tryCatch(
    match.arg(sim_ka$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_ka$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  if (route == "oral" &&
      (is.null(sim_ka$values) ||
       length(sim_ka$values) == 0 || all(is.na(sim_ka$values)))) {
    stop("No candidate Ka values available for parameter sweeping.",
         call. = FALSE)
  }

  # --- sim_cl ---
  sim_cl$mode <- tryCatch(
    match.arg(sim_cl$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_cl$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  if (is.null(sim_cl$values) ||
      length(sim_cl$values) == 0 || all(is.na(sim_cl$values))) {
    stop("No candidate CL values available for parameter sweeping.",
         call. = FALSE)
  }

  # --- sim_vc ---
  sim_vc$mode <- tryCatch(
    match.arg(sim_vc$mode, choices = c("manual")),
    error = function(e) {
      stop("Invalid `sim_vc$mode`. Must be 'manual'.", call. = FALSE)
    }
  )
  vc_values <- sim_vc$values
  if (is.null(vc_values) ||
      length(vc_values) == 0 || all(is.na(vc_values))) {
    stop("No candidate Vc values available for parameter sweeping.",
         call. = FALSE)
  }

  vc_values <- na.omit(vc_values) %>%
    sort() %>%
    tibble::tibble(value = .) %>%
    dplyr::mutate(prev = dplyr::lag(value),
                  rel_diff = abs(value - prev) / prev) %>%
    dplyr::filter(is.na(rel_diff) | rel_diff > 0.2) %>%
    dplyr::pull(value)

  # --- sim_vp & sim_vp2 ---
  sim_vp$mode <- tryCatch(
    match.arg(sim_vp$mode, choices = c("auto", "manual")),
    error = function(e) {
      stop("Invalid `sim_vp$mode`. Must be 'auto' or 'manual'.", call. = FALSE)
    }
  )
  sim_vp2$mode <- tryCatch(
    match.arg(sim_vp2$mode, choices = c("auto", "manual")),
    error = function(e) {
      stop("Invalid `sim_vp2$mode`. Must be 'auto' or 'manual'.", call. = FALSE)
    }
  )

  if (sim_vp$mode == "manual" &&
      (is.null(sim_vp$values) || all(is.na(sim_vp$values)))) {
    stop("No candidate Vp1 values in manual mode.", call. = FALSE)
  }
  if (sim_vp2$mode == "manual" &&
      (is.null(sim_vp2$values) || all(is.na(sim_vp2$values)))) {
    stop("No candidate Vp2 values in manual mode.", call. = FALSE)
  }

  if (sim_vp$mode == "auto" && sim_vp2$mode == "auto") {
    vc_vp_ratio_range <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
    combs <-
      expand.grid(estvc = vc_values, vc_vp_ratio = vc_vp_ratio_range)
    combs$estvp <- signif(combs$estvc / combs$vc_vp_ratio)
    combs2 <- combs
    combs$estvp2 <- combs$estvc
    combs2$estvp2 <- combs2$estvp
    diagonal_combs <-
      combs %>% dplyr::transmute(estvc = estvc,
                                 estvp = estvp,
                                 estvp2 = estvp)
    combs_df <- dplyr::bind_rows(combs, combs2, diagonal_combs) %>%
      dplyr::distinct(estvc, estvp, estvp2, .keep_all = TRUE) %>%
      dplyr::transmute(Vc = estvc, Vp1 = estvp, Vp2 = estvp2)

  } else if (sim_vp$mode == "manual" && sim_vp2$mode == "manual") {
    combs_df <-
      expand.grid(Vc = vc_values,
                  Vp1 = sim_vp$values,
                  Vp2 = sim_vp2$values)

  } else if (sim_vp$mode == "manual" && sim_vp2$mode == "auto") {
    vp2_df <-
      expand.grid(vc = vc_values, ratio = c(10, 5, 2, 1, 0.5, 0.2, 0.1)) %>%
      dplyr::mutate(vp2 = vc / ratio)
    combs_df <- expand.grid(Vc = vc_values, Vp1 = sim_vp$values) %>%
      dplyr::left_join(vp2_df %>% dplyr::rename(Vc = vc), by = "Vc") %>%
      dplyr::rename(Vp2 = vp2)

  } else if (sim_vp$mode == "auto" && sim_vp2$mode == "manual") {
    vp1_df <-
      expand.grid(vc = vc_values, ratio = c(10, 5, 2, 1, 0.5, 0.2, 0.1)) %>%
      dplyr::mutate(vp1 = vc / ratio)
    combs_df <-
      expand.grid(Vc = vc_values, Vp2 = sim_vp2$values) %>%
      dplyr::left_join(vp1_df %>% dplyr::rename(Vc = vc), by = "Vc") %>%
      dplyr::rename(Vp1 = vp1)
  }

  # --- sim_q ---
  sim_q$mode <- tryCatch(
    match.arg(sim_q$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_q$mode`. Must be 'manual' or 'auto'.", call. = FALSE)
    }
  )
  sim_q$auto.strategy <- tryCatch(
    match.arg(sim_q$auto.strategy, choices = c("scaled", "fixed")),
    error = function(e) {
      stop("Invalid `sim_q$auto.strategy`. Must be 'scaled' or 'fixed'.",
           call. = FALSE)
    }
  )
  q1_values <- if (sim_q$mode == "auto") {
    cl_val <- sim_cl$values[1]
    if (sim_q$auto.strategy == "scaled")
      c(cl_val / 2, cl_val, cl_val * 2)
    else
      c(1, 10, 100)
  } else {
    if (is.null(sim_q$values) || all(is.na(sim_q$values))) {
      stop("No candidate Q values available for parameter sweeping.",
           call. = FALSE)
    }
    sim_q$values
  }

  # --- sim_q2 ---
  sim_q2$mode <- tryCatch(
    match.arg(sim_q2$mode, choices = c("manual", "auto")),
    error = function(e) {
      stop("Invalid `sim_q2$mode`. Must be 'manual' or 'auto'.", call. = FALSE)
    }
  )
  sim_q2$auto.strategy <- tryCatch(
    match.arg(sim_q2$auto.strategy, choices = c("scaled", "fixed")),
    error = function(e) {
      stop("Invalid `sim_q2$auto.strategy`. Must be 'scaled' or 'fixed'.",
           call. = FALSE)
    }
  )

  q2_strategy_scaled <- FALSE
  q2_values <- if (sim_q2$mode == "auto") {
    if (sim_q2$auto.strategy == "scaled") {
      q2_strategy_scaled <- TRUE
      NULL
    } else {
      c(1, 10, 100)
    }
  } else {
    if (is.null(sim_q2$values) || all(is.na(sim_q2$values))) {
      stop("No candidate Q2 values available for parameter sweeping.",
           call. = FALSE)
    }
    sim_q2$values
  }

  # --- Param grid ---
  if (q2_strategy_scaled) {
    param_grid <- tidyr::crossing(
      combs_df,
      CL = sim_cl$values,
      Q1 = q1_values,
      Ka = if (route == "oral")
        sim_ka$values
      else
        NA_real_
    )
    param_grid$Q2 <- param_grid$Q1
    param_grid <-
      param_grid %>% dplyr::select(Vc, Vp1, Vp2, Q1, Q2, CL, Ka)
  } else {
    param_grid <- tidyr::crossing(
      combs_df,
      CL = sim_cl$values,
      Q1 = q1_values,
      Q2 = q2_values,
      Ka = if (route == "oral")
        sim_ka$values
      else
        NA_real_
    ) %>% dplyr::select(Vc, Vp1, Vp2, Q1, Q2, CL, Ka)
  }

  # --- Simulations ---
  start_time <- Sys.time()
  progressr::handlers(
    progressr::handler_progress(format = ":message [:bar] :percent (:current/:total)", width = 80)
  )
  results <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(param_grid))
    param_grid %>%
      dplyr::mutate(row = dplyr::row_number()) %>%
      purrr::pmap_dfr(function(Vc, Vp1, Vp2, Q1, Q2, CL, Ka, row) {
        p(sprintf(
          "Running simulation: Vc=%.0f Vp1=%.0f Vp2=%.0f",
          Vc,
          Vp1,
          Vp2
        ))
        sim_out <- if (route == "iv") {
          Fit_3cmpt_iv(dat[dat$EVID != 2, ], "rxSolve", CL, Vc, Vp1, Vp2, Q1, Q2, input.add = 0)
        } else {
          Fit_3cmpt_oral(dat[dat$EVID != 2, ], "rxSolve", Ka, CL, Vc, Vp1, Vp2, Q1, Q2, input.add = 0)
        }
        met <-
          metrics.(pred.x = sim_out$cp, obs.y = dat[dat$EVID == 0,]$DV)
        elapsed <-
          round(difftime(Sys.time(), start_time, units = "secs"), 2)

        tibble::tibble(
          Vc = Vc,
          Vp1 = Vp1,
          Vp2 = Vp2,
          Q1 = Q1,
          Q2 = Q2,
          CL = CL,
          Ka = if (route == "oral")
            Ka
          else
            NA_real_,
          APE = round(met[1], 2),
          MAE = round(met[2], 2),
          MAPE = round(met[3], 2),
          RMSE = round(met[4], 2),
          rRMSE1 = round(met[5], 2),
          rRMSE2 = round(met[6], 2),
          Cumulative.Time.Sec = as.numeric(elapsed)
        )
      })
  })

  return(results)
}
