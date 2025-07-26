#' Approximate volume of distribution based on Cmax
#'
#' Approximates the volume of distribution (\eqn{V_d}) using observed \eqn{C_{\mathrm{max}}} values.
#' The calculation assumes rapid absorption (\eqn{k_a \gg k_e}), meaning that \eqn{C_{\mathrm{max}}} occurs shortly after dosing, before significant elimination has taken place.
#' To ensure this assumption, \eqn{C_{\mathrm{max}}} points are filtered to those occurring within \eqn{0.2 \cdot \mathrm{half-life}}, a time window where elimination is minimal.
#'
#' @param dat A data frame containing pharmacokinetic data, including observed concentrations (\code{DV}), time after dose (\code{tad}), dose, and route information.
#' @param half_life The half-life of the drug, used to filter \eqn{C_{\mathrm{max}}} points to the early phase of elimination.
#' @param route The route of administration, either `"bolus"` (default) or `"infusion"`. Determines the formula for \(V_d\) calculation.
#' @return The trimmed geometric mean of the approximated \(V_d\), calculated across individuals.
#' @details
#' - **Rapid Absorption Assumption:**
#'   When \eqn{k_a \gg k_e}, most of the absorption is completed before elimination significantly impacts the drug concentration. Under this assumption, \eqn{C_{\mathrm{max}}} is dominated by absorption rather than elimination.
#'
#' - **Why \eqn{0.2 \cdot \mathrm{half-life}}:**
#'   - At \eqn{0.2 \cdot \mathrm{half-life}}, less than 13% of the drug is eliminated. This is derived from the intravenous (IV) exponential decay model:
#'     \deqn{\mathrm{Elimination\ Fraction} = 1 - e^{-k_e \cdot t}}
#'     Substituting \eqn{t = 0.2 \cdot \mathrm{half-life}} and \eqn{k_e = \ln(2) / \mathrm{half-life}}:
#'     \deqn{\mathrm{Elimination\ Fraction} = 1 - e^{-\ln(2) \cdot 0.2} \approx 0.13}
#'     This ensures the selected \eqn{C_{\mathrm{max}}} points primarily reflect absorption dynamics.
#'
#' - **Data Filtering:**
#'   - Observed data (\code{EVID == 0}) is grouped by \code{ID} and \code{dose_number}.
#'   - The \eqn{C_{\mathrm{max}}} for each group is identified and further filtered to include only points where \eqn{tad < 0.2 \cdot \mathrm{half-life}}.
#'
#' - **Volume of Distribution Calculation (single-dose):**
#'   - For bolus administration:
#'     \deqn{V_d = \frac{\mathrm{Dose}}{C_{\mathrm{max}}}}
#'   - For infusion:
#'     \deqn{V_d = \frac{\mathrm{Rate} \cdot \min(\mathrm{Time}, \mathrm{Duration})}{C_{\mathrm{max}}}}
#'
#' - **Summary of \eqn{V_d} values:**
#'   - Individual \eqn{V_d} values are summarized using a trimmed geometric mean (\eqn{10\%}) to reduce the impact of outliers.
#'
#' - **Accumulation Ratio (\eqn{R_{\mathrm{ac}}}):**
#'   - For multiple-dose scenarios, the accumulation ratio is calculated as:
#'     \deqn{R_{\mathrm{ac}} = \frac{1}{1 - e^{-k_e \cdot \tau}}}
#'     where:
#'       - \eqn{k_e = \frac{\ln(2)}{\mathrm{half-life}}} is the elimination rate constant.
#'       - \eqn{\tau} is the dosing interval (\code{dose_interval}).
#'   - The accumulation ratio represents the steady-state concentration relative to the single-dose concentration.
#'
#' - **Volume of Distribution Calculation (multiple-dose):**
#'   - For infusion routes:
#'     \deqn{V_d = \frac{\mathrm{Rate} \cdot \min(\mathrm{Time}, \mathrm{Duration})}{C_{\mathrm{max}} / R_{\mathrm{ac}}}}
#'       - \eqn{\mathrm{Rate}} is the infusion rate.
#'       - \eqn{\min(\mathrm{Time}, \mathrm{Duration})} represents the shorter of infusion time or observation time.
#'       - \eqn{C_{\mathrm{max}} / R_{\mathrm{ac}}} adjusts the observed concentration to reflect single-dose dynamics.
#'   - For bolus or other routes:
#'     \deqn{V_d = \frac{\mathrm{Dose}}{C_{\mathrm{max}} / R_{\mathrm{ac}}}}
#'       - \eqn{\mathrm{Dose}} is the administered dose.
#'
#'
#' @examples
#' \dontrun{
#' # Example 1: Bolus administration
#' dat<-Bolus_2CPT
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' approx.vc(dat=fdat,half_life = half_life,route=froute)
#'
#' # Example 2: Oral administration
#' dat<-Oral_2CPT
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' approx.vc(dat=fdat,half_life = half_life,route=froute)
#'
#' # Example 3: Theno_sd
#' dat<-theo_sd
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' approx.vc(dat=fdat,half_life = half_life,route=froute)
#'
#' # Example 3: Theno_md
#' dat<-theo_md
#' out <- processData(dat)
#' fdat<- out$dat
#' froute <-out$Datainfo$Value[out$Datainfo$Infometrics == "Dose Route"]
#' half_life <- get_hf(dat = fdat)$half_life_median
#' approx.vc(dat=fdat,half_life = half_life,route=froute)
#' }
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
        paste(shQuote(c("bolus", "oral", "infusion")), collapse = ", ")
      ), call. = FALSE)
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
    # Calculate median vd for each individual
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
