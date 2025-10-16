#' Control settings for fallback rules in parameter estimation
#'
#' @param enable_ka_fallback Logical. Whether to fallback to ka = 1 if invalid.
#' @param sigma_method_additive Method for additive sigma ("model" or "fixed_fraction").
#' @param sigma_method_proportional Method for proportional sigma ("model" or "fixed_fraction").
#' @param sigma_fallback_fraction Numeric. e.g. 0.2 (20% of DV mean).
#'
#' @return A list of fallback control parameters.
#' @export
#' @examples
#' fallback_control(sigma_method_additive = "fixed_fraction", sigma_fallback_fraction = 0.25)
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

#' Create full control list for initial PPK parameter estimation
#'
#' Aggregates all modular control functions into one structured list.
#'
#' @param ss.control Output of `ss_control()`.
#' @param pooled.control Output of `pooled_control()`.
#' @param nca.control Output of `nca_control()`.
#' @param fallback.control Output of `fallback_control()`.
#' @param selmetrics A character string (or vector) specifying which model performance
#' metric(s) to evaluate. Must be one or more of: \code{"APE"}, \code{"MAE"},
#' \code{"MAPE"}, \code{"RMSE"}, \code{"rRMSE1"}, \code{"rRMSE2"}. Default is
#' \code{"rRMSE2"}.
#' @param hybrid.base Logical. If \code{TRUE}, enables hybrid evaluation mode,
#' where model performance is assessed using mixed parameter combinations across
#' methods (e.g., combining \code{ka} from NCA, \code{cl} from graphic methods,
#' and \code{vd} from adaptive single-point method). If \code{FALSE}, evaluation
#' is performed separately for each parameter estimation method (e.g., using only
#' \code{simpcal}, only \code{nca}, etc.). This allows more flexible exploration
#' of parameter set synergies across estimation strategies. Default is \code{TRUE}.
#' @param preferNCA Logical.
#' If \code{TRUE} and \code{selmetrics == "rRMSE2"}, the function will first
#' select the best method by the lowest rRMSE2. If that method is not NCA-based,
#' it will then check whether a purely NCA-based method exists with a lower
#' rRMSE1 than the rRMSE2-best method. If such a candidate exists, the NCA
#' method will be selected instead. This option expresses a preference for NCA
#' methods without enforcing them. Default is \code{TRUE}.
#'
#' @return A named list combining all control modules.
#' @export
#' @examples
#' initsControl(
#'   pooled.control = pooled_control(nbins = 8),
#'   fallback.control = fallback_control(sigma_method_additive = "fixed_fraction")
#' )
#'
#' @seealso \code{\link{ss_control}}, \code{\link{pooled_control}}, \code{\link{nca_control}}, \code{\link{fallback_control}}
#'
initsControl <- function(ss.control = ss_control(),
                         pooled.control = pooled_control(),
                         nca.control = nca_control(),
                         fallback.control = fallback_control(),
                         selmetrics= "rRMSE2",
                         hybrid.base = TRUE,
                         preferNCA=TRUE) {
  list(
    ss.control = ss.control,
    pooled.control = pooled.control,
    nca.control = nca.control,
    fallback.control = fallback.control,
    selmetrics= selmetrics,
    hybrid.base =hybrid.base,
    preferNCA=preferNCA
  )
}

