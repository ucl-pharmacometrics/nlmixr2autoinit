#' Control settings for fallback rules in parameter estimation
#'
#' @param enable_ka_fallback Logical value indicating whether to apply a fallback
#'   to ka = 1 if the estimated value is invalid.
#' @param sigma_method_additive Method for additive sigma. Options are
#'   "model" or "fixed_fraction".
#' @param sigma_method_proportional Method for proportional sigma. Options are
#'   "model" or "fixed_fraction".
#' @param sigma_fallback_fraction Numeric value specifying the fallback fraction,
#'   for example, 0.2 corresponds to 20 percent of the mean of observed concentrations.
#'
#' @return A list of fallback control parameters.
#'
#' @examples
#' fallback_control()
#'
#' @export

fallback_control <- function(enable_ka_fallback = FALSE,
                             sigma_method_additive = "model",
                             sigma_method_proportional = "model",
                             sigma_fallback_fraction = 0.2) {
  list(
    enable_ka_fallback = enable_ka_fallback,
    sigma_method_additive = sigma_method_additive,
    sigma_method_proportional = sigma_method_proportional,
    sigma_fallback_fraction = sigma_fallback_fraction

  )
}

#' Create full control list for initial parameter estimation
#'
#' Aggregates modular control functions into a structured list for use in
#' population pharmacokinetic parameter initialization.
#'
#' @param ss.control A control list consistent with the structure returned by
#'   ss_control().
#' @param pooled.control A control list consistent with the structure returned
#'   by pooled_control().
#' @param nca.control A control list consistent with the structure returned by
#'   nca_control().
#' @param fallback.control A control list consistent with the structure
#'   returned by fallback_control().
#'
#' @param selmetrics A character string or vector specifying model performance
#'   metrics to evaluate. Must be one or more of "APE", "MAE", "MAPE", "RMSE",
#'   "rRMSE1", or "rRMSE2". Default is "rRMSE2".
#'
#' @param hybrid.base Logical. If TRUE, enables hybrid evaluation mode in which
#'   model performance is assessed using mixed parameter combinations across
#'   methods. If FALSE, each method is evaluated independently. Default is TRUE.
#'
#' @param preferNCA Logical. If TRUE and selmetrics equals "rRMSE2", the lowest
#'   rRMSE2 is selected first. If the best-performing method is not NCA-based,
#'   the function then checks whether an NCA-based method offers a lower rRMSE1.
#'   If so, the NCA method is selected. Default is TRUE.
#'
#' @return A named list combining all control modules for parameter estimation.
#'
#' @examples
#' initsControl(
#'   pooled.control = pooled_control(nbins = 8),
#'   fallback.control = fallback_control(
#'     sigma_method_additive = "fixed_fraction"
#'   )
#' )
#'
#' @seealso \link{ss_control}, \link{pooled_control}, \link{nca_control},
#'   \link{fallback_control}
#'
#' @export

initsControl <- function(ss.control = ss_control(),
                         pooled.control = pooled_control(),
                         nca.control = nca_control(),
                         fallback.control = fallback_control(),
                         selmetrics = "rRMSE2",
                         hybrid.base = TRUE,
                         preferNCA = TRUE) {
  list(
    ss.control = ss.control,
    pooled.control = pooled.control,
    nca.control = nca.control,
    fallback.control = fallback.control,
    selmetrics = selmetrics,
    hybrid.base = hybrid.base,
    preferNCA = preferNCA
  )
}

