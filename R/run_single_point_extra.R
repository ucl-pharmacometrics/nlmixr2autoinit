#' Perform extended single-point pharmacokinetic calculations
#'
#' Extends the single-point pharmacokinetic calculations by incorporating additional logic to derive clearance (\eqn{CL}), volume of distribution (\eqn{V_d}), and absorption rate constant (\eqn{k_a}) based on the availability of data. The function evaluates which parameters were not successfully calculated in the \code{single_point_base} step and uses appropriate methods to estimate missing parameters.
#'
#' @param single_point_base.lst A list object returned by \code{single_point_base}, containing preprocessed data and initial calculations for \eqn{CL} and \eqn{V_d}.
#' @param dat A data frame of pharmacokinetic observations.
#' @param half_life Numeric. half-life (in the same time units as
#' the data). Used to derive missing \eqn{CL} or \eqn{V_d} when needed and to
#' restrict oral \eqn{k_a} calculations to the absorption phase.
#' @param route Character. Administration route; one of \code{"bolus"},
#' \code{"oral"}, or \code{"infusion"}. Determines which single-point rules and
#' fallbacks are applied.
#' @param dose_type Optional character string indicating the dosing scenario for
#' alignment with pooled processing, e.g., \code{"first_dose"},
#' \code{"repeated_doses"}, or \code{"combined_doses"}.
#' @param pooled_ctrl A list of pooled-analysis control options created by
#' \code{\link{pooled_control}}. Controls time binning and TAD rounding used if
#' additional pooled summaries are needed.
#' @param ssctrl A list of steady-state control options created by
#' \code{\link{ss_control}}. Governs assumptions and thresholds used when
#' interpreting or approximating steady-state behavior in single-point logic.

#' @return A list containing:
#' \itemize{
#'   \item \code{singlepoint.results}: A data frame with estimated \eqn{k_a}, \eqn{CL}, \eqn{V_d}, and processing time.
#'   \item \code{dat}: The input dataset used in the calculations.
#'   \item \code{single_point_ka_df}: Data used for \eqn{k_a} calculations (if applicable).
#'   \item \code{single_point_cl_df}: Data used for \eqn{CL} calculations.
#'   \item \code{single_point_vd_df}: Data used for \eqn{V_d} calculations.
#'   \item \code{single_point_vd_cmax_df}: \eqn{C_{\mathrm{max}}} data used for \eqn{V_d} estimation.
#' }
#'
#' @details
#' The function uses a series of conditional logic to determine which pharmacokinetic parameters have not been calculated in the \code{single_point_base} part and derives them using the parameters that have already been calculated:
#'
#' - **Successfully Calculated \eqn{CL} and \eqn{V_d}:**
#'   - If both \eqn{CL} and \eqn{V_d} are successfully calculated in the \code{single_point_base} part, they are directly used in the results without further derivations.
#'
#' - **Uncalculated \eqn{V_d}:**
#'   - If \eqn{V_d} cannot be calculated in the \code{single_point_base} part but \eqn{CL} and \eqn{t_{1/2}} are available, \eqn{V_d} is derived using:
#'     \deqn{V_d = \frac{CL \cdot t_{1/2}}{\ln(2)}}
#'
#' - **Uncalculated \eqn{CL}:**
#'   - If \eqn{CL} cannot be calculated in the \code{single_point_base} part but \eqn{V_d} and \eqn{t_{1/2}} are available, \eqn{CL} is derived using:
#'     \deqn{CL = \frac{V_d \cdot \ln(2)}{t_{1/2}}}
#'
#' - **Uncalculated Both \eqn{CL} and \eqn{V_d}:**
#'   - If neither \eqn{CL} nor \eqn{V_d} can be calculated in the \code{single_point_base} part, \eqn{V_d} is estimated from \eqn{C_{\mathrm{max}}} and dose:
#'     \deqn{V_d = \frac{\mathrm{Dose}}{C_{\mathrm{max}}}}
#'     \eqn{CL} is then derived from \eqn{V_d} and \eqn{t_{1/2}}.
#'
#' - **Oral Absorption Rate (\eqn{k_a}):**
#'   - For oral dosing, the calculation focuses on data collected before \eqn{t_{\mathrm{max}}}, ensuring that \eqn{\mathrm{tad} < 0.2 \cdot \mathrm{half-life}} to reflect the absorption phase.
#'   - The absorption rate constant (\eqn{k_a}) is computed using the \code{run_ka_solution} function, which applies a solution-based approach tailored to the selected data points.
#'
#' @examples
#'
#' # Example usage with preprocessed theo_md data
#' dat <- Bolus_1CPT
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' run_single_point_extra( route=froute,half_life=half_life, single_point_base.lst=run_single_point_base(dat = fdat, half_life = half_life,route=froute))
#'
#' @export
#'
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
        !is.na(trimmed_mean_vd) ~ trimmed_mean_vd,!is.na(trimmed_mean_cl) &
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
          "CL and Vd were calculated directly from steady-state and single-dose data.",!is.na(trimmed_mean_cl) &
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
      datobs <- dat[dat$EVID == 0,]

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
                              cmax_by_group$iiobs == 0,]) > 0) {
        single_point.ka.out <- run_ka_solution(df = cmax_by_group,
                                               cl = trimmed_mean_cl,
                                               ke = trimmed_mean_cl / trimmed_mean_vd)
        median_ka<-single_point.ka.out$ka_calc_median
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



#' Run full adaptive single-point PK analysis (CL, Vd, Ka)
#'
#' This function runs both base and extended single-point pharmacokinetic calculations
#' to estimate clearance (CL), volume of distribution (Vd), and absorption rate constant (Ka) when possible.
#'
#' @param dat A pharmacokinetic dataset.
#' @param route Character string: one of `"bolus"`, `"infusion"`, `"oral"`.
#' @param half_life Optional numeric value for the elimination half-life. Estimated if NULL.
#' @param dose_type Optional character string (e.g. `"bolus"` or `"infusion"`). Required if `half_life` is NULL.
#' @param pooled_ctrl Control list from \code{\link{pooled_control}()}.
#' @param ssctrl Control list from \code{\link{ss_control}()}.
#'
#' @return An object of class `"single.point.lst"` containing all results from base and extended estimation.
#'
#' @examples
#' \dontrun{
#' # Example: Adaptive single-point PK analysis for bolus data
#' dat <- Bolus_1CPT
#' # Step 1: Preprocess the data
#' out <- processData(dat)
#' fdat <- out$dat

#' # Step 2: Extract route and dose type info
#' froute <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' fdose_type <- out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Type"]

#' # Step 3: Estimate half life
#' half_life <- get_hf(dat = fdat)$half_life_median

#' # Step 3: Run single-point analysis (CL, Vd, Ka if oral)
#' result <- run_single_point(
#'   dat = fdat,
#'  route = froute,
#'  dose_type = fdose_type,
#'   half_life = half_life
#' )
#' # Step 4: View results
#' print(result$singlepoint.results)
#' head(result$single_point_cl_df)  # Data used for CL estimation
#' head(result$single_point_vd_df)  # Data used for Vd estimation
#' #' head(result$single_point_ka_df)  # Ka (only if oral)
#' }
#'
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




