#' Perform extended single-point pharmacokinetic calculations
#'
#' Extends the base single-point pharmacokinetic calculations by incorporating
#' additional logic to derive clearance (CL), volume of distribution (Vd), and absorption rate
#' constant (ka) based on availability of pharmacokinetic data.
#
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @param half_life Optional numeric value representing the elimination half-life of the drug.
#'   If not provided, half-life is estimated within run_single_point_base() using get_hf() applied
#'   to the pharmacokinetic observations.
#'
#' @param single_point_base.lst A list object returned from the base single-point calculation
#'   that includes input data, preprocessing information, and initial estimates of CL and Vd.
#'   If not supplied, the function will internally call the base routine using dat.
#'
#' @param route Character string specifying the route of administration. Must be one of
#'   "bolus", "oral", or "infusion".
#'
#' @param dose_type Classified as "first_dose", "repeated_doses", or "combined_doses"
#'   based on whether observed concentrations occur following the first administration,
#'   during repeated dosing, or across both contexts. This parameter is passed directly to
#'   run_single_point_base().
#'
#' @param pooled_ctrl A list of pooled-analysis control options created using pooled_control().
#'   These control time binning and time-after-dose rounding when pooled summaries are required.
#'
#' @param ssctrl A list of steady-state control options created using ss_control().
#'   These govern assumptions and thresholds used when interpreting steady-state behavior
#'   in single-point logic.
#'
#' @return A list containing:
#' \itemize{
#'   \item singlepoint.results: A data frame with estimated ka, CL, Vd, computation start time,
#'   processing time, and a descriptive message of the applied logic.
#'   \item dat: The dataset used for processing.
#'   \item single_point.ka.out: Output used for estimating the absorption rate constant
#'   (for oral administration), if applicable.
#'   \item single_point_cl_df: Data used for clearance estimation.
#'   \item single_point_vd_df: Data used for volume of distribution estimation.
#'   \item approx.vc.out: Data used for estimating the volume of distribution from maximum
#'   concentration (Cmax) and dose when needed.
#' }
#'
#' @details
#' The function derives pharmacokinetic parameters using the following logic:
#'
#' - When both clearance (CL) and volume of distribution (Vd) are available from the base
#'   calculation, these values are directly used without modification.
#'
#' - If Vd is missing but CL and elimination half-life are provided, Vd is calculated using:
#'   \deqn{V_d = \frac{CL \cdot t_{1/2}}{\ln(2)}}
#'
#' - If CL is missing but Vd and half-life are available, CL is calculated using:
#'   \deqn{CL = \frac{V_d \cdot \ln(2)}{t_{1/2}}}
#'
#' - If both CL and Vd are unavailable but half-life is provided, Vd is estimated using dose and Cmax:
#'   \deqn{V_d = \frac{\mathrm{Dose}}{C_{\mathrm{max}}}}
#'   and CL is subsequently derived:
#'   \deqn{CL = \frac{V_d \cdot \ln(2)}{t_{1/2}}}
#'
#' - For oral administration, the absorption rate constant (ka) is estimated using concentration-time
#'   data collected prior to Tmax and restricted to the absorption phase, defined as time after dose
#'   (tad) less than 0.2 multiplied by the half-life. The value of ka is obtained using a solution-based
#'   estimation method.
#'
#' The function supports bolus, oral, and infusion administration routes. For oral dosing,
#' only data within the absorption phase are used to estimate the absorption rate constant (ka),
#' specifically using concentration-time observations prior to the maximum concentration (Tmax).
#'
#' @seealso
#'   \itemize{
#'     \item run_single_point_base(): Performs the initial single-point estimation of CL and Vd.
#'     \item run_single_point(): Wrapper function that sequentially runs both base and extended steps.
#'   }
#'
#' @examples
#' dat <- Bolus_1CPT
#' out <- processData(dat)
#' fdat <- out$dat
#' froute <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' run_single_point_extra(
#'   dat = fdat,
#'   half_life = half_life,
#'   single_point_base.lst = run_single_point_base(
#'     dat = fdat,
#'     half_life = half_life,
#'     route = froute
#'   ),
#'   route = froute
#' )
#'
#' @export

run_single_point_extra <- function(dat = NULL,
                                   half_life = NULL,
                                   single_point_base.lst = NULL,
                                   route = c("bolus", "oral", "infusion"),
                                   dose_type = NULL,
                                   pooled_ctrl = pooled_control(),
                                   ssctrl = ss_control()) {
  `%>%` <- magrittr::`%>%`

  start.time <-
    single_point_base.lst$summary$start.time

  # If base result not supplied, compute it
  if (is.null(single_point_base.lst)) {
    if (is.null(dat))
      stop("You must supply either`dat`.")
    single_point_base.lst <- run_single_point_base(
      dat = dat,
      half_life = half_life,
      route = route,
      dose_type = dose_type,
      pooled_ctrl = pooled_ctrl,
      ssctrl = ssctrl
    )
  }

  median_ka <- NA
  single_point.ka.out <- NA

  # obtain run_single_point_base information
  dat <- single_point_base.lst$dat
  trimmed_mean_cl <-
    single_point_base.lst$summary$cl
  trimmed_mean_vd <-
    single_point_base.lst$summary$vd

  dat.ss.obs <- single_point_base.lst$cl_df
  dat.fd.obs <- single_point_base.lst$vd_df

  # obtain volume of distribution through cmax
  approx.vc.out <- approx.vc(
    dat = dat,
    route = route,
    single_point_base.lst = single_point_base.lst,
    half_life = half_life
  )

  cl_vd_final <- tibble::tibble(
    trimmed_mean_cl = trimmed_mean_cl,
    trimmed_mean_vd = trimmed_mean_vd,
    half_life = half_life,
    approx_vd = approx.vc.out$approx.vc.value
  ) %>%
    dplyr::mutate(
      trimmed_mean_vd = dplyr::case_when(
        !is.na(trimmed_mean_vd) ~ trimmed_mean_vd,
        !is.na(trimmed_mean_cl) &
          !is.na(half_life) ~
          signif(trimmed_mean_cl * half_life / log(2), 3),
        is.na(trimmed_mean_cl) &
          is.na(trimmed_mean_vd) & !is.na(half_life) ~
          approx_vd,
        TRUE ~ NA_real_
      ),
      trimmed_mean_cl = dplyr::case_when(
        !is.na(trimmed_mean_cl) ~ trimmed_mean_cl,
        is.na(trimmed_mean_cl) &
          !is.na(trimmed_mean_vd) & !is.na(half_life) ~
          signif(trimmed_mean_vd * log(2) / half_life, 3),
        is.na(trimmed_mean_cl) &
          is.na(trimmed_mean_vd) & !is.na(half_life) ~
          signif(approx_vd * log(2) / half_life, 3),
        TRUE ~ NA_real_
      ),
      single_point.message = dplyr::case_when(
        !is.na(trimmed_mean_cl) & !is.na(trimmed_mean_vd) ~
          "CL and Vd were calculated directly from steady-state and single-dose data.",
        !is.na(trimmed_mean_cl) &
          is.na(trimmed_mean_vd) & !is.na(half_life) ~
          "CL from steady-state; Vd derived from CL and estimated half-life.",
        is.na(trimmed_mean_cl) &
          !is.na(trimmed_mean_vd) & !is.na(half_life) ~
          "Vd from early data; CL derived from Vd and estimated half-life.",
        is.na(trimmed_mean_cl) &
          is.na(trimmed_mean_vd) & !is.na(half_life) ~
          "Vd from Cmax; CL derived from Vd and estimated half-life.",
        TRUE ~ "Unable to compute CL or Vd"
      )
    )

  trimmed_mean_cl <- cl_vd_final$trimmed_mean_cl
  trimmed_mean_vd <- cl_vd_final$trimmed_mean_vd
  single_point.message <- cl_vd_final$single_point.message

  # extra ka for oral case
  if (route == "oral") {
    if (is.na(trimmed_mean_cl) == F & is.na(trimmed_mean_vd) == F) {
      datobs <- dat[dat$EVID == 0, ]

      cmax_by_group <- datobs %>%
        dplyr::group_by(ID, dose_number) %>%
        dplyr::mutate(Tmax = TIME[which.max(DV)]) %>%
        dplyr::filter(TIME < Tmax |
                        (
                          dplyr::n() == 1 &
                            tad < 0.2 * log(2) * trimmed_mean_cl / trimmed_mean_vd
                        )) %>%
        dplyr::slice_max(order_by = DV, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Tmax)

      if (nrow(cmax_by_group[cmax_by_group$EVID == 0  &
                             cmax_by_group$iiobs == 0, ]) > 0) {
        single_point.ka.out <- run_ka_solution(df = cmax_by_group,
                                               cl = trimmed_mean_cl,
                                               ke = trimmed_mean_cl / trimmed_mean_vd)
        median_ka <- single_point.ka.out$ka_calc_median
      }
    }
  }

  end.time <- Sys.time()

  time.spent <-
    round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

  # Only selected the key columns
  singlepoint.results <- data.frame(
    ka = signif(median_ka, 3),
    cl = signif(trimmed_mean_cl, 3),
    vd = signif(trimmed_mean_vd, 3),
    starttime = start.time,
    time.spent = time.spent,
    single_point.message = single_point.message
  )

  single.point.lst <- list(
    singlepoint.results = singlepoint.results,
    dat = dat,
    single_point.ka.out =   single_point.ka.out,
    single_point_cl_df =  dat.ss.obs,
    single_point_vd_df =  dat.fd.obs,
    approx.vc.out = approx.vc.out
  )

  class(single.point.lst) <- "single.point.lst"

  return(single.point.lst)
}



#' Run full adaptive single-point PK analysis (CL, Vd, ka)
#'
#' Runs the base and extended single-point pharmacokinetic calculations to estimate
#' clearance (CL), volume of distribution (Vd), and, when applicable, the absorption
#' rate constant (ka).
#'
#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#'
#' @param route Character string specifying the route of administration. Must be one of
#'   "bolus", "oral", or "infusion".
#'
#' @param half_life Optional numeric value representing the elimination half-life of the drug.
#'   If not provided, half-life is estimated within run_single_point_base() using get_hf() applied
#'   to the pharmacokinetic observations.
#'
#' @param dose_type Classified as "first_dose", "repeated_doses", or "combined_doses"
#'   based on whether observed concentrations occur following the first administration,
#'   during repeated dosing, or across both contexts. This parameter is passed directly to
#'   run_single_point_base().
#'
#' @param pooled_ctrl A list of pooled-analysis control options created using pooled_control().
#'   These control time binning and time-after-dose rounding when pooled summaries are required.
#'
#' @param ssctrl A list of steady-state control options created using ss_control().
#'   These govern assumptions and thresholds used when interpreting steady-state behavior
#'   in single-point logic.
#'
#' @return An object of class "single.point.lst" with results from the base and extended steps.
#'
#' @examples
#' # Example: Adaptive single-point PK analysis for bolus data
#' # Step 1: Preprocess the data
#' dat <- Oral_1CPT
#' out <- processData(dat)
#' # Step 2: Extract route and dose type info
#' froute <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' fdose_type <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Type"]
#' # Step 3: Estimate half life
#' half_life <- get_hf(dat = out$dat)$half_life_median
#' # Step 4: Run single-point analysis (CL, Vd, Ka if oral)
#' result <- run_single_point(
#'   dat = out$dat,
#'   route = froute,
#'   dose_type = fdose_type,
#'   half_life = half_life
#' )
#' # Step 5: View results
#' print(result$singlepoint.results)
#'
#' @seealso
#'   \itemize{
#'     \item run_single_point_base(): Performs the initial single-point estimation of CL and Vd.
#'     \item run_single_point_extra(): Extends the base estimation to derive additional parameters such as ka.
#'   }
#'
#' @export
run_single_point <- function(dat,
                             route = c("bolus", "oral", "infusion"),
                             half_life = NULL,
                             dose_type = NULL,
                             pooled_ctrl = pooled_control(),
                             ssctrl = ss_control()) {
  # Run base first
  base_out <- run_single_point_base(
    dat = dat,
    route = route,
    half_life = half_life,
    dose_type = dose_type,
    pooled_ctrl = pooled_ctrl,
    ssctrl = ssctrl
  )

  # Then run extra step (uses base output)
  extra_out <- run_single_point_extra(
    dat = dat,
    half_life = half_life,
    route = route,
    dose_type = dose_type,
    pooled_ctrl = pooled_ctrl,
    ssctrl = ssctrl,
    single_point_base.lst = base_out
  )

  return(extra_out)
}




