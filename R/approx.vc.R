#' Approximate volume of distribution from observed Cmax
#'
#' Estimates the volume of distribution (\eqn{V_d}) from observed peak
#' concentrations (\eqn{C_{\mathrm{max}}}) in single-dose, multiple-dose, or mixed datasets.
#'
#' @param dat A data frame containing pharmacokinetic data, including observed
#' concentrations (DV), time after dose (tad), dose, and route information.
#' @param half_life The elimination half-life (\eqn{t_{1/2}}) of the compound,
#' used to identify early-phase \eqn{C_{\mathrm{max}}} values.
#' @param single_point_base.lst Optional list object returned by
#' \link{run_single_point_base}(). If not supplied, the function will generate it internally.
#' @param route Route of administration. One of "bolus", "oral",
#' or "infusion" (default = "bolus").
#' @param dose_type Optional string specifying the dosing type, passed to
#' \link{run_single_point_base}().
#' @param pooled_ctrl Control object created by \link{pooled_control}(),
#' defining data pooling options.
#' @param ssctrl Control object created by \link{ss_control}(),
#' defining steady-state control options.
#'
#' @return A list containing individual and population Vd estimates and related dose-level data.
#'
#' @details
#' Estimates individual apparent volumes of distribution (\eqn{V_d})
#' from observed peak concentrations (\eqn{C_{\mathrm{max}}}).
#' Individual estimates are then summarized to obtain a population-level value.
#'
#' For single-dose data, \eqn{V_d} is calculated according to the route of administration:
#' \itemize{
#'   \item Bolus: \eqn{V_d = \mathrm{Dose} / C_{\mathrm{max}}}
#'   \item Infusion: \eqn{V_d = (\mathrm{Rate} \times t_{\mathrm{inf}}) / C_{\mathrm{max}}}
#'   \item Oral: \eqn{V_d = (\mathrm{Dose} \times F) / C_{\mathrm{max}}}, where \eqn{F = 1 - e^{-k_a t}}
#' }
#'
#' For multiple-dose data, observed \eqn{C_{\mathrm{max}}} values are adjusted to
#' single-dose equivalents using the accumulation ratio:
#' \deqn{R_{\mathrm{ac}} = \frac{1}{1 - e^{-k_e \tau}}, \quad k_e = \ln(2)/t_{1/2}}
#'
#' Adjusted values are used to estimate \eqn{V_d} using the same route-specific equations.
#'
#' The function returns separate subsets for single- and multiple-dose data
#' (`approx.vc.dat.sd` and `approx.vc.dat.md`), their combined dataset (`approx.vc.dat`),
#' and the final summary statistic (`approx.vc.value`).
#'
#' @seealso
#' \code{\link{run_single_point_base}}, \code{\link{trimmed_geom_mean}},
#' \code{\link{pooled_control}}, \code{\link{ss_control}}
#'
#' @examples
#' out <- processData(Bolus_1CPT)
#' approx.vc(
#'   dat = out$dat,
#'   half_life = get_hf(dat = out$dat)$half_life_median,
#'   route = out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' )$approx.vc.value
#' @export

approx.vc <- function(dat = NULL,
                      half_life = NULL,
                      single_point_base.lst = NULL,
                      route = c("bolus", "oral", "infusion"),
                      dose_type = NULL,
                      pooled_ctrl = pooled_control(),
                      ssctrl = ss_control()) {
  # Ensure half-life is provided
  if (is.null(half_life)) {
    stop("No half-life provided for approximate volume of distribution calculation.")
  }

  # Defensive check for route
  route <- tryCatch(
    match.arg(route, choices = c("bolus", "oral", "infusion")),
    error = function(e) {
      stop(sprintf(
        "Invalid `route`: '%s'. Must be one of: %s.",
        as.character(route),
        paste(shQuote(c(
          "bolus", "oral", "infusion"
        )), collapse = ", ")
      ),
      call. = FALSE)
    }
  )

  # Run base analysis if not supplied
  if (is.null(single_point_base.lst)) {
    if (is.null(dat)) {
      stop("You must supply either `single_point_base.lst` or `dat`.")
    }

    single_point_base.lst <- run_single_point_base(
      dat = dat,
      half_life = half_life,
      route = route,
      dose_type = dose_type,
      pooled_ctrl = pooled_ctrl,
      ssctrl = ssctrl
    )
  }

  cmax_by_group_sd <- NULL
  cmax_by_group_md <- NULL
  approx.vc.value <- NA

  if (any(dat$EVID == 0 & dat$dose_number == 1, na.rm = TRUE)) {
    datobs_fd <- dat[dat$EVID == 0 & dat$dose_number == 1,]
    cmax_by_group_sd <- datobs_fd %>%
      dplyr::group_by(ID, dose_number) %>%
      dplyr::slice_max(order_by = DV, with_ties = FALSE) %>%
      dplyr::ungroup()

    cmax_by_group_sd <-
      cmax_by_group_sd[cmax_by_group_sd$tad < half_life * 0.2,]

    if (route == "infusion") {
      cmax_by_group_sd$vd <-
        signif(
          pmin(cmax_by_group_sd$TIME, cmax_by_group_sd$durationobs) *
            cmax_by_group_sd$rateobs / cmax_by_group_sd$DV,
          3
        )
    } else if (route == "bolus") {
      cmax_by_group_sd$vd <-
        signif(cmax_by_group_sd$dose / cmax_by_group_sd$DV, 3)
    } else if (route == "oral") {
      # Compute absorbed fraction: F = 1 - exp(-Ka * t)
      ka_val <- 1.0 # assumed
      absorbed_frac <- 1 - exp(-ka_val * cmax_by_group_sd$TIME)
      cmax_by_group_sd$vd <-
        signif((cmax_by_group_sd$dose * absorbed_frac) / cmax_by_group_sd$DV,
               3)
    }
  }

  # Accumulation ratio was borrowed to approximately estimate the Cmax after single dose
  if (length(single_point_base.lst$cl_df) > 1) {
    dat.ss.obs <- single_point_base.lst$cl_df
    dat.ss.obs$Rac <-
      1 / (1 - exp(-(log(2) / half_life) * dat.ss.obs$dose_interval))

    cmax_by_group_md <-  dat.ss.obs %>%
      dplyr::group_by(ID, dose_number) %>%
      dplyr::slice_max(order_by = DV, with_ties = FALSE) %>%
      dplyr::ungroup()

    cmax_by_group_md <-
      cmax_by_group_md[cmax_by_group_md$tad < half_life * 0.2, ]
    if (route == "infusion") {
      cmax_by_group_md$vd <- signif(
        pmin(cmax_by_group_md$tad, cmax_by_group_md$durationobs) *
          cmax_by_group_md$rateobs /
          (cmax_by_group_md$DV / cmax_by_group_md$Rac),
        3
      )
    } else if (route == "bolus") {
      cmax_by_group_md$vd <-
        signif(cmax_by_group_md$dose /
                 (cmax_by_group_md$DV / cmax_by_group_md$Rac),
               3)
    } else if (route == "oral") {
      ka_val <- 1.0 # assumed
      absorbed_frac <- 1 - exp(-ka_val * cmax_by_group_md$tad)
      cmax_by_group_md$vd <- signif(
        (cmax_by_group_md$dose * absorbed_frac) /
          (cmax_by_group_md$DV / cmax_by_group_md$Rac),
        3
      )
    }
  }

  cmax_by_group <- rbind(
    data.frame(
      ID = cmax_by_group_sd$ID,
      dose_number = cmax_by_group_sd$dose_number,
      vd = cmax_by_group_sd$vd
    ),
    data.frame(
      ID = cmax_by_group_md$ID,
      dose_number = cmax_by_group_md$dose_number,
      vd = cmax_by_group_md$vd
    )
  )

  if (!is.null(cmax_by_group)) {
    # Calculate median Vd for each individual
    individual_mean_vd <-
      tryCatch(
        aggregate(vd ~ ID, data = cmax_by_group, FUN = trimmed_geom_mean),
        error = function(e) {
          NA
        }
      )
    # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
    trimmed_mean_vd <-
      tryCatch(
        trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE),
        error = function(e) {
          NA
        }
      )
    approx.vc.value <- trimmed_mean_vd
  }
  return(
    list = list(
      approx.vc.value = approx.vc.value,
      approx.vc.dat = cmax_by_group,
      approx.vc.dat.sd = cmax_by_group_sd,
      approx.vc.dat.md = cmax_by_group_md
    )
  )

}
