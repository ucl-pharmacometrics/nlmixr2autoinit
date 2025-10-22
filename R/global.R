# ================================================================
# Global variable declarations
# Each section corresponds to a specific function module.
# ================================================================

if (getRversion() >= "3.5.0") {

  # ------------------------ approx.vc ----------------------------
  utils::globalVariables(c(
    "ID", "dose_number", "DV"
  ))

  # ------------------------ bin.time -----------------------------
  utils::globalVariables(c(
    "EVID", "tad", "DVstd", "median"
  ))

  # ------------------------ calculate_cl -------------------------
  utils::globalVariables(c(
    "SteadyState", "ID", "dose_number", "DV", "tad",
    "max_value", "min_value", "max_interval", "min_interval",
    "TIME", "dose", "Css_avg_i", "recent_ii"
  ))

  # ------------------------ calculate_tad ------------------------
  utils::globalVariables(c(
    "DOSE", "dose", "ID", "resetflag", "TIME", "EVID", "dose_number",
    "AMT", "II", "RATE", "route", "last_dose_time", "last_dose_number",
    "last_dose", "last_ii", "last_rate", "last_route", "routeobs"
  ))

  # ------------------------ calculate_vd -------------------------
  utils::globalVariables(c(
    "EVID", "dose_number", "tad", "iiobs", "ID", "C_first_flag",
    "has_any_C_first", "TIME", "min_time", "dose", "DV",
    "durationobs", "rateobs", "vd"
  ))

  # ------------------------ getPPKinits --------------------------
  utils::globalVariables(c(
    "ID", "resetflag", "data_type",
    "Relative Root Mean Squared Error (rRMSE1)",
    "Relative Root Mean Squared Error (rRMSE2)",
    "Mean Absolute Percentage Error (MAPE)",
    "CL Method", "Vd Method", "min_count", "trim"
  ))

  # ------------------------ get_pooled_data ----------------------
  utils::globalVariables(c(
    "EVID", "AMT", "ID", "TIME", "interval", "tad",
    "dose_number", "tad_check", "iiobs", "resetflag"
  ))

  # ------------------------ getnca / getsigma / getsigmas --------
  utils::globalVariables(c(
    "EVID", "ID", "resetflag", "dose_number"
  ))

  # ------------------------ is_ss -------------------------------
  utils::globalVariables(c(
    "ID", "TIME", "dose_time_hist", "last_two_doses_interval",
    "doses_required", "dose_amt_hist", "doses_to_check", "median",
    "intervals", "amts_to_check", "tad", "dose_interval", "SSflag",
    "EVID", "dose_count_before_obs", "is_continuous", "is_same_dose",
    "is_within_last_dose_interval"
  ))

  # ------------------------ mark_dose_number ---------------------
  utils::globalVariables(c(
    "ID", "resetflag", "TIME", "CMT", "EVID"
  ))

  # ------------------------ processData --------------------------
  utils::globalVariables(c(
    ".data", ".", "ID", "reset_event", "EVID", "TIME", "SS", "II",
    "SS_events", "resetflag", "CMT", "obs_cmt", "n_unique_cmt",
    "rateobs", "dose", "dose_number", "DV", "n_post_Tmax", "Tmax",
    "iiobs", "AMT", "route"
  ))

  # ------------------------ run_ka_solution ----------------------
  utils::globalVariables(c(
    "ID", "TIME", "DV", "Tmax", "dose_number", "EVID", "SteadyState"
  ))

  # ------------------------ run_single_point_extra ---------------
  utils::globalVariables(c(
    "ID", "dose_number", "TIME", "DV", "Tmax", "tad"
  ))

  # ------------------------ sim_sens_1cmpt_mm --------------------
  utils::globalVariables(c(
    ".", "value", "prev", "rel_diff", "Vmax", "Km"
  ))

  # ------------------------ sim_sens_2cmpt -----------------------
  utils::globalVariables(c(
    ".", "value", "prev", "rel_diff", "vc", "vp",
    "Vc", "Vp", "Q", "CL", "Ka"
  ))

  # ------------------------ sim_sens_3cmpt -----------------------
  utils::globalVariables(c(
    ".", "value", "prev", "rel_diff", "estvc", "estvp", "estvp2",
    "vc", "ratio", "vp1", "vp2", "Vp1", "Vp2", "Q1", "Q2", "CL", "Ka"
  ))

  # ------------------------ Shared / general-use variables -------
  utils::globalVariables(c(
    "DVstd", "DOSE", "vd", "interval", "durationobs", "Rac",
    "Css_avg_i", "recent_ii", "SteadyState"
  ))

  # ------------------------ Fit function  -----------------------
  utils::globalVariables(c("/", "/<-", "dt", "centre", "depot", "A1", "A2", "A3"))

  # ------------------------ Model  -----------------------
  utils::globalVariables(c("ini", "model"))
}
