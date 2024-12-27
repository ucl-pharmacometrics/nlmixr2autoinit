#' Run single-point pharmacokinetic analysis
#'
#' This function performs a single-point pharmacokinetic analysis, estimating key parameters such as clearance (\eqn{CL}), volume of distribution (\eqn{V_d}), and absorption rate constant (\eqn{k_a}) based on the provided dataset and optional half-life (\eqn{t_{1/2}}).
#'
#' @param dat A data frame containing pharmacokinetic data, including observed concentrations (\code{DV}), time after dose (\code{tad}), dosing interval (\code{ii}), and dose information.
#' @param half_life The half-life of the drug (\eqn{t_{1/2}}). If not provided, alternative assumptions or data-based rules may be applied for calculations.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{singlepoint.results}: A data frame with estimated \eqn{CL}, \eqn{V_d}, \eqn{k_a}, processing time, and a descriptive message.
#'   \item \code{dat}: The input dataset used in the calculations.
#'   \item \code{single_point_ka_df}: Data used for \eqn{k_a} calculations (if applicable).
#'   \item \code{single_point_cl_df}: Data used for \eqn{CL} calculations.
#'   \item \code{single_point_vd_df}: Data used for \eqn{V_d} calculations.
#'   \item \code{approx.vc.out}: A message or result indicating whether the approximate central volume of distribution (\eqn{V_d}) was successfully estimated.
#' }
#'
#' @details
#' This function integrates two key components:
#' \itemize{
#'   \item \code{single_point_base}: Conducts basic single-point calculations to estimate parameters based on steady-state or dose-interval data.
#'   \item \code{single_point_extra}: Extends the analysis to handle additional cases, such as calculating oral absorption rate (\eqn{k_a}) or estimating parameters when data is incomplete (e.g., when \eqn{CL} or \eqn{V_d} is missing).
#' }
#'
#'
#' @examples
#' # Example 1 (iv case)
#' dat <- Bolus_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' run_single_point(fdat,half_life)
#'
#' # Example 2 (infusion case).
#'
#' dat <- Infusion_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' run_single_point(fdat,half_life)
#'
#' dat <- Oral_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' run_single_point(fdat,half_life)
#'
#' @export


run_single_point <- function(dat,
                              half_life=NA) {

  single.point.base.out<-single_point_base(dat,half_life)
  single.point.out<-single_point_extra(single_point_base.lst = single.point.base.out,half_life = half_life)

  return(single.point.out)
}

#' Perform basic single-point pharmacokinetic calculations
#'
#' Performs pharmacokinetic calculations using a single-point method and calculates clearance (CL) and volume of distribution (Vd) based on the provided pharmacokinetic data.
#'
#' @details
#' The function includes two main calculations:
#'
#' 1. **Clearance (CL):**
#'    - Calculated using steady-state data from multiple dosing regimens.
#'    - Steady-state concentrations are identified and used to estimate clearance based on the dose, observed steady-state concentrations, and dosing intervals.
#'
#' 2. **Volume of Distribution (Vd):**
#'    - Calculated for intravenous single-dose data.
#'    - The first sampling point (C0) after dose administration is used to estimate Vd with the formula `Dose / C0`.
#'
#' The function supports both bolus and short-term infusion routes. It processes input data to identify steady-state and first-dose conditions and calculates trimmed mean values for CL and Vd to reduce the impact of outliers.
#'
#' @param dat A data frame containing pharmacokinetic data, including information on doses, sampling times, and observed concentrations.
#' @param half_life The half-life of the drug, used for determining elimination phase parameters and steady-state conditions.
#' @return A list containing the following:
#'   - `singlepoint.results`: A data frame with the calculated clearance (CL) and volume of distribution (Vd), along with processing start and elapsed time.
#'   - `dat`: The processed dataset with additional columns for steady-state and other calculated flags.
#'   - `dat.ss.obs`: Subset of the data used for steady-state clearance calculation.
#'   - `dat.fd.obs`: Subset of the data used for single-dose volume of distribution calculation.
#' @importFrom dplyr %>% mutate if_else arrange group_by ungroup slice_min filter select
#' @import nlmixr2
#'
#' @examples
#' # Example 1 (iv case)
#' dat <- Bolus_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' run_single_point(fdat,half_life)
#'
#' # Example 2 (infusion case).
#'
#' dat <- Infusion_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' single_point_base(fdat,half_life)
#' run_single_point(fdat,half_life)
#'
#'  # Example 3 (Oral case).
#'
#' dat <- Oral_1CPT
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' run_single_point(fdat,half_life)
#'
#'
#' @export

single_point_base <- function(dat,
                              half_life=NA) {
  start.time <- Sys.time()
  # Default
  trimmed_mean_ka <- NA
  trimmed_mean_cl <- NA
  trimmed_mean_vd <- NA
  ka_single_point.out<-NA
  dat.ss.obs <- NA
  dat.fd.obs <- NA
  cmax_by_group<-NA
  single_point.message<-NA

##############################Single Point Base#############################

 ##################################Calculate of clearance######################
  dat <- is_ss(df = dat,
                 half_life = half_life )


  dat.ss.obs <- dat[dat$SteadyState == TRUE & !is.na(dat$SteadyState), ]

  if (nrow(dat.ss.obs) > 2) {
      # If there are multiple points within the same dose interval,
      # only the minimum and maximum values are selected as the steady-state points for calculation
      dat.ss.obs <- dat.ss.obs %>%
        group_by(ID, dose_number) %>%
        mutate(
          max_value = max(DV),
          min_value = min(DV),
          max_time = tad[which.max(DV)],
          min_time = tad[which.min(DV)],
          avg_value = (max_value + min_value) / 2,
          max_interval = ifelse(DV == max(DV), TRUE, FALSE),
          min_interval = ifelse(DV == min(DV), TRUE, FALSE)
        ) %>%
        ungroup()

      # Extract rows to be marked in dat.ss.obs
      rows_to_mark <- dat.ss.obs %>%
        filter(max_interval == TRUE | min_interval == TRUE) %>%
        select(ID, dose_number, TIME)

      dat <- dat %>%
        rowwise() %>%
        mutate(SteadyState = ifelse(
          any(
            ID == rows_to_mark$ID &
              dose_number == rows_to_mark$dose_number &
              TIME == rows_to_mark$TIME
          ),
          TRUE,
          FALSE
        )) %>%
        ungroup()

      # Same for dat.ss.obs
      dat.ss.obs  <- dat.ss.obs  %>%
        rowwise() %>%
        mutate(SteadyState = ifelse(
          any(
            ID == rows_to_mark$ID &
              dose_number == rows_to_mark$dose_number &
              TIME == rows_to_mark$TIME
          ),
          TRUE,
          FALSE
        )) %>%
        ungroup()

      # Second selection， only select the max, min points
      dat.ss.obs <- dat.ss.obs[dat.ss.obs$SteadyState == T,]

      # avg_value type identification
      # Only applicable for oral case with very fast absorption, and with the majority of absorption occurring before Tmax, and only a short period is required to reach Tmax.

      dat.ss.obs  <-   dat.ss.obs %>%
        mutate(
          # Initialise 'Css_type' to "Css_avg" by default
          Css_type = "Css_avg",

          # Identify Cssmax and Cssmin based on 'tad', 'ii', and additional conditions
          Css_type = case_when(

            # Condition 1: Classify as Css_avg based on exponential condition
            # If e^(-k * tau) >= 0.6667, it implies:
            # 1. Css_max and Css_min differ by less than 1.5 times.
            # 2. Css_avg and Css_min differ by less than 20%.
            # 3. Css_avg and Css_max differ by less than 20%.
            exp(-log(2)/half_life * recent_ii) >= 0.6667 ~ "Css_avg",

            # Condition 2: Classify as Css_avg based on dynamic threshold
            # Use the theoretical decay model to determine if the ratio between
            # max_value and min_value exceeds the expected range based on half-life.
            max_value / min_value > exp(-log(2)/half_life * recent_ii) ~ "Css_avg",

            # Condition 3: Classify as Css_max
            # - Both max_time and min_time must be within the first 20% of the dosing interval (recent_ii).
            # - Neither max_time nor min_time should be 0 (to avoid misclassifying early absorption as Css_max).
            (max_time <= 0.2 * recent_ii & max_time != 0) &
              (min_time <= 0.2 * recent_ii & min_time != 0) ~ "Css_max",

            # Condition 4: Classify as Css_min
            # - Both max_time and min_time must be within the last 80%-100% of the dosing interval (recent_ii).
            # - Either max_time or min_time being 0 (indicating very low concentration) is also classified as Css_min.
            # - Note: The original formula was based on (recent_ii - 0.2 * half_life), but due to concerns about half-life stability,
            #   it has been replaced with recent_ii for better robustness.
            (max_time >= 0.8 * recent_ii | max_time == 0) &
              (min_time >= 0.8 * recent_ii | min_time == 0) ~ "Css_min",

            # Default condition: Classify as Css_avg for all other cases
            TRUE ~ Css_type
          )
        )


      dat.ss.obs$cl <- NA

      # Calculate the clearance by single point and Css_type

      for (i in 1:nrow(dat.ss.obs)) {
        # Only one points in the dose interval
        if (dat.ss.obs[i,]$Css_type=="Css_avg") {
          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / dat.ss.obs[i,]$avg_value /dat.ss.obs[i,]$recent_ii ,
              digits = 3
            )
        }


        if (dat.ss.obs[i,]$Css_type=="Css_max" & dat.ss.obs[i,]$routeobs=="bolus") {

          Css_min_i<-dat.ss.obs[i,]$max_value * exp(-((log(2)/as.numeric(half_life))*dat.ss.obs[i,]$recent_ii))
          Css_avg_i<-(dat.ss.obs[i,]$max_value+   Css_min_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )
        }

        if (dat.ss.obs[i,]$Css_type=="Css_max" & dat.ss.obs[i,]$routeobs=="infusion") {

          Css_min_i<-dat.ss.obs[i,]$max_value * exp(-((log(2)/as.numeric(half_life))*(dat.ss.obs[i,]$recent_ii-dat.ss.obs[i,]$duration_obs)))
          Css_avg_i<-(dat.ss.obs[i,]$max_value+   Css_min_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )
        }

        if (dat.ss.obs[i,]$Css_type=="Css_min" & dat.ss.obs[i,]$routeobs=="bolus") {

          Css_max_i<-dat.ss.obs[i,]$min_value / exp(-((log(2)/as.numeric(half_life))*dat.ss.obs[i,]$recent_ii))
          Css_avg_i<-(dat.ss.obs[i,]$min_value + Css_max_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )

        }


        if (dat.ss.obs[i,]$Css_type=="Css_min" & dat.ss.obs[i,]$routeobs=="infusion") {

          Css_max_i<-dat.ss.obs[i,]$min_value / exp(-((log(2)/as.numeric(half_life))*(dat.ss.obs[i,]$recent_ii-dat.ss.obs[i,]$duration_obs)))
          Css_avg_i<-(dat.ss.obs[i,]$min_value + Css_max_i)/2

          dat.ss.obs[i,]$cl <-
            signif(
              dat.ss.obs[i,]$dose / Css_avg_i/ dat.ss.obs[i,]$recent_ii,
              digits = 3
            )

        }

      } # end loop


      # Calculate geometric mean cl for each individual
      individual_mean_cl <- tryCatch(aggregate(cl ~ ID, data = dat.ss.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_cl <- tryCatch(trimmed_geom_mean(individual_mean_cl$cl, trim=0.05,na.rm = TRUE),error=function(e) {NA})

   }

  ########################### Calculate the Volume of distribution###############

  if (unique(dat[dat$EVID==1,]$route) == "bolus") {

    # remove ssflag=1 case
    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 & dat$tad < half_life * 0.2 & dat$iiobs==0,]) > 0) {

      dat$C_first_flag <- 0
      # dat[dat$EVID == 0 &
      #       dat$dose_number == 1 & dat$tad < half_life * 0.2 & dat$iiobs==0,]$C_first_flag <- 1

      # only retain the first point per ID with C_first_flag = 1
      dat <- dat %>%
        mutate(
          C_first_flag = ifelse(
            EVID == 0 & dose_number == 1 & tad < half_life * 0.2 & iiobs == 0, 1, 0
          )
        ) %>%
        group_by(ID) %>%
        mutate(
          C_first_flag = ifelse(C_first_flag == 1 & TIME == min(TIME[C_first_flag == 1], na.rm = TRUE), 1, 0)
        ) %>%
        ungroup()

    dat.fd.obs <- dat[dat$C_first_flag == 1, ]

    dat.fd.obs$vd <- signif(dat.fd.obs$dose / dat.fd.obs$DV, 3)

    # Calculate geometric mean vd for each individual
    # individual_mean_vd <- aggregate(vd ~ ID, data = dat.fd.obs, FUN = trimmed_geom_mean)

    # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
    # trimmed_mean_vd <-trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE)
    trimmed_mean_vd <-trimmed_geom_mean(dat.fd.obs$vd, trim = 0.05, na.rm = TRUE)
    }
  }

  # Short-term infusion assumed
  if (unique(dat[dat$EVID==1,]$route) == "infusion") {
    # remove ssflag=1 case
    if (nrow(dat[dat$EVID == 0 &
                 dat$dose_number == 1 &
                 dat$tad < half_life * 0.2 &
                 # dat$TIME < dat$duration_obs & # allow sampling outside infusion duration
                 dat$iiobs==0,]) > 0) {

      dat$C_first_flag <- 0

      dat[dat$EVID == 0 &
            dat$dose_number == 1 &
            dat$tad < half_life * 0.2 &
            # dat$TIME < dat$duration_obs  &
            dat$iiobs==0,]$C_first_flag <- 1

      dat.fd.obs <- dat[dat$C_first_flag == 1, ]
      dat.fd.obs$vd <-
        signif(pmin(dat.fd.obs$TIME,dat.fd.obs$duration_obs) * dat.fd.obs$rateobs / dat.fd.obs$DV, 3)

      # Calculate median vd for each individual
      individual_mean_vd <- tryCatch( aggregate(vd ~ ID, data = dat.fd.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

      # Calculate the trimmed mean (e.g., 10% trimmed mean to reduce outlier impact)
      trimmed_mean_vd <-tryCatch(trimmed_geom_mean(individual_mean_vd$vd, trim = 0.05, na.rm = TRUE),error=function(e) {NA})

    }
  }

  # Only selected the key columns
  single_point_base.results <- data.frame(
    cl = signif( trimmed_mean_cl, 3),
    vd = signif( trimmed_mean_vd, 3),
    starttime = start.time
  )

  return(
    list(
      single_point_base.results= single_point_base.results,
      dat = dat,
      single_point_cl_df =  dat.ss.obs,
      single_point_vd_df =  dat.fd.obs
    )
  )

}

#' Perform extended single-point pharmacokinetic calculations
#'
#' Extends the single-point pharmacokinetic calculations by incorporating additional logic to derive clearance (\eqn{CL}), volume of distribution (\eqn{V_d}), and absorption rate constant (\eqn{k_a}) based on the availability of data. The function evaluates which parameters were not successfully calculated in the \code{single_point_base} step and uses appropriate methods to estimate missing parameters.
#'
#' @param single_point_base.lst A list object returned by \code{single_point_base}, containing preprocessed data and initial calculations for \eqn{CL} and \eqn{V_d}.
#'
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
#'
#' dat <- theo_md
#' fdat <- processData(dat)$dat
#' half_life<-half_life_estimated(dat = fdat)$half_life_median
#' single_point_extra(single_point_base(fdat,half_life))
#'
#'
#' @export
#'

single_point_extra <- function(single_point_base.lst,
                               half_life) {
  dat <- single_point_base.lst$dat
  trimmed_mean_ka <- NA
  ka_single_point.out<-NA
  trimmed_mean_cl <-
    single_point_base.lst$single_point_base.results$cl
  trimmed_mean_vd <-
    single_point_base.lst$single_point_base.results$vd
  dat.ss.obs <- single_point_base.lst$single_point_cl_df
  dat.fd.obs <- single_point_base.lst$single_point_vd_df
  start.time <-
    single_point_base.lst$single_point_base.results$starttime

  approx.vc.out <- approx.vc(
    dat = dat,
    single_point_base.lst = single_point_base.lst,
    half_life = half_life
  )

  ##############################Single Point Extra#############################
  # Both can be calculated in the base part
  if (is.na(trimmed_mean_cl) == F & is.na(trimmed_mean_vd) == F) {
    single_point.message <-
      "Clearance was calculated using steady-state data and volume of distribution was calculated based on data within the dosing interval following the first dose "
  }

  # if half_life is available, Vd is not available
  if (is.na(trimmed_mean_cl) == F &
      is.na(trimmed_mean_vd) == T & is.na(half_life) == F) {
    message(black(
      paste0(
        "Insufficient data points to support Vd calculation (single-dose IV or oral); derived from clearance and estimated half-life instead.",
        strrep(".", 20)
      )
    ))

    # individual_mean_vd <- tryCatch( aggregate(cl* half_life/log(2) ~ ID, data = dat.ss.obs, FUN = trimmed_geom_mean),error=function(e) {NA})

    trimmed_mean_vd <-  signif(trimmed_mean_cl * half_life / log(2), 3)
    single_point.message <-
      "Clearance was calculated using steady-state data and volume of distribution was calculated based on clearance and estimated half_life"
  }

  # single-point method volume of distribution + half-life estimated
  if (is.na(trimmed_mean_cl) == T &
      is.na(trimmed_mean_vd) == F & is.na(half_life) == F) {
    message(black(
      paste0(
        "Insufficient steady-state data for CL calculation; derived from Vd and estimated half-life instead.",
        strrep(".", 20)
      )
    ))

    trimmed_mean_cl <-
      signif(trimmed_mean_vd * log(2) / half_life, 3)
    single_point.message <-
      "Vd was calculated based on data within the dosing interval following the first dose and clearance was calculated based on volume of distribution and estimated half_life"
  }

  if (is.na(trimmed_mean_cl) == T &
      is.na(trimmed_mean_vd) == T & is.na(half_life) == F) {
    message(black(
      paste0(
        "Neither single-dose (IV) data nor steady-state data supports the calculation of clearance (CL) and volume of distribution (Vd). Vd will be estimated using Dose/Cmax, and CL will be derived based on the estimated half-life and Vd.",
        strrep(".", 20)
      )
    ))

    trimmed_mean_vd <- approx.vc.out$approx.vc.value

    trimmed_mean_cl <-
      signif(trimmed_mean_vd * log(2) / half_life, 3)

    single_point.message <-
      "Vd was calculated based on Cmax in each dose interval and clearance was calculated based on volume of distribution and estimated half_life"
  }

  # extra ka for oral case
  if (unique(dat[dat$EVID == 1, ]$route) == "oral") {
    if (is.na(trimmed_mean_cl) == F & is.na(trimmed_mean_vd) == F) {
      datobs <- dat[dat$EVID == 0, ]

      cmax_by_group2 <- datobs %>%
        group_by(ID, dose_number) %>%
        mutate(Tmax = TIME[which.max(DV)]) %>%
        filter(TIME < Tmax |
                 (
                   n() == 1 & tad < 0.2 * log(2) * trimmed_mean_cl / trimmed_mean_vd
                 )) %>%
        slice_max(order_by = DV, with_ties = FALSE) %>%
        ungroup() %>%
        select(-Tmax)

      if (nrow(cmax_by_group2[cmax_by_group2$EVID == 0  &
                              cmax_by_group2$iiobs == 0, ]) > 0) {
        ka_single_point.out <- run_ka_solution(df = cmax_by_group2,
                                               cl = trimmed_mean_cl,
                                               ke = trimmed_mean_cl / trimmed_mean_vd)

        trimmed_mean_ka <-
          tryCatch(
            trimmed_geom_mean(
              ka_single_point.out$ka_calc_dat$ka_calcv,
              trim = 0.05,
              na.rm = TRUE
            ),
            error = function(e) {
              NA
            }
          )

      }
    }
  }

  end.time <- Sys.time()

  time.spent <- round(difftime(end.time, start.time), 4)

  # Only selected the key columns
  singlepoint.results <- data.frame(
    ka = signif(trimmed_mean_ka, 3),
    cl = signif(trimmed_mean_cl, 3),
    vd = signif(trimmed_mean_vd, 3),
    starttime = start.time,
    time.spent = time.spent,
    single_point.message = single_point.message
  )

  single.point.lst<- list(
      singlepoint.results = singlepoint.results,
      dat = dat,
      single_point_ka_df =   ka_single_point.out,
      single_point_cl_df =  dat.ss.obs,
      single_point_vd_df =  dat.fd.obs,
      approx.vc.out = approx.vc.out
    )

   class(single.point.lst) <- "single.point.lst"

  return(
    single.point.lst
  )

}






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
#' # Example 1: Bolus administration
#' dat<-Bolus_2CPT
#' dat<-processData(dat)$dat
#' approx.vc(dat,half_life = half_life_estimated(dat)$half_life_median)
#'
#' # Example 2: Oral administration
#' dat<-Oral_2CPT
#' dat<-processData(dat)$dat
#' approx.vc(dat,half_life = half_life_estimated(dat)$half_life_median)
#'
#' # Example 3: Theno_sd
#' dat<-theo_sd
#' dat<-processData(dat)$dat
#' approx.vc(dat,half_life = half_life_estimated(dat)$half_life_median)
#'
#' # Example 3: Theno_md
#' dat<-theo_md
#' dat<-processData(dat)$dat
#' approx.vc(dat,half_life = half_life_estimated(dat)$half_life_median)
#'
#'
#' @export

approx.vc <- function(dat,
                      half_life,
                      single_point_base.lst) {
  if (missing(half_life) || is.na(half_life) || is.null(half_life)) {
    stop("No half life provided for approximate central volume of distribution calculation")
  }

  if (missing(single_point_base.lst)) {
    single_point_base.lst <- single_point_base(dat, half_life)
  }

  cmax_by_group1 <- NULL
  cmax_by_group2 <- NULL
  approx.vc.value <- NA

  if (nrow(dat[dat$EVID == 0 & dat$dose_number == 1, ]) > 0) {
    datobs_fd <- dat[dat$EVID == 0 & dat$dose_number == 1, ]

    cmax_by_group1 <- datobs_fd %>%
      group_by(ID, dose_number) %>%
      slice_max(order_by = DV, with_ties = FALSE) %>%
      ungroup()

    cmax_by_group1 <- cmax_by_group1[cmax_by_group1$tad < half_life * 0.2, ]

    if (unique(dat[dat$EVID == 1, ]$route) == "infusion") {
      cmax_by_group1$vd <-
        signif(
          pmin(cmax_by_group1$TIME, cmax_by_group1$duration_obs) * cmax_by_group1$rateobs / cmax_by_group1$DV,
          3
        )
    }

    cmax_by_group1$vd <-
      signif(cmax_by_group1$dose / cmax_by_group1$DV, 3)
  }

  # Accumulation ratio was borrowed to approximately estimate the Cmax after single dose
  if (length(single_point_base.lst$single_point_cl_df) > 1) {
    dat.ss.obs <- single_point_base.lst$single_point_cl_df

    dat.ss.obs$Rac <-
      1 / (1 - exp(-(log(2) / half_life) * dat.ss.obs$dose_interval))

    cmax_by_group2 <-  dat.ss.obs %>%
      group_by(ID, dose_number) %>%
      slice_max(order_by = DV, with_ties = FALSE) %>%
      ungroup()

    cmax_by_group2 <-
      cmax_by_group2[cmax_by_group2$tad < half_life * 0.2, ]

    if (unique(dat[dat$EVID == 1, ]$route) == "infusion") {
      cmax_by_group2$vd <-
        signif(
          pmin(cmax_by_group2$TIME, cmax_by_group2$duration_obs) * cmax_by_group2$rateobs / (cmax_by_group2$DV /
                                                                                               cmax_by_group2$Rac),
          3
        )
    }

    cmax_by_group2$vd <-
      signif(cmax_by_group2$dose / (cmax_by_group2$DV / cmax_by_group2$Rac),
             3)
  }

  cmax_by_group <- rbind(
    data.frame(
      ID = cmax_by_group1$ID,
      dose_number = cmax_by_group1$dose_number,
      vd = cmax_by_group1$vd
    ),
    data.frame(
      ID = cmax_by_group2$ID,
      dose_number = cmax_by_group2$dose_number,
      vd = cmax_by_group2$vd
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
      approx.vc.dat.sd = cmax_by_group1,
      approx.vc.dat.md = cmax_by_group2
    )
  )

}



#' Calculate the trimmed geometric mean
#'
#' Computes the trimmed geometric mean of a numeric vector.
#' A specified proportion of the smallest and largest values is removed before calculating the geometric mean.
#'
#' @param x A numeric vector containing the values for which the trimmed geometric mean is to be calculated.
#' @param trim A numeric value between 0 and 0.5 indicating the proportion of observations
#'   to be trimmed from each end of the vector. Default is 0 (no trimming).
#' @param na.rm Logical. Should missing values (NA) be removed before computation? Defaults to TRUE.
#' @return The trimmed geometric mean of the input vector as a numeric value.
#'   If no values remain after trimming, returns `NA`.
#'
trimmed_geom_mean <- function(x, trim = 0, na.rm = TRUE) {
  if (na.rm) x <- na.omit(x)  # Remove NA values if na.rm is TRUE

  # Sort the vector and apply trimming
  x <- sort(x)
  n <- length(x)
  lower <- floor(n * trim) + 1  # Index of first value to keep
  upper <- n - floor(n * trim)  # Index of last value to keep

  # Keep only trimmed range
  x_trimmed <- x[lower:upper]

  # Calculate the geometric mean of the trimmed values
  exp(mean(log(x_trimmed), na.rm = FALSE))
}
