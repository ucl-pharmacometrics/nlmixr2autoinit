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
#' @param pooled_control Output of `pooled_control()`.
#' @param nca_control Output of `nca_control()`.
#' @param selection_control Output of `selection_control()`.
#' @param fallback_control Output of `fallback_control()`.
#' @param selmetrics A character string (or vector) specifying which model performance
#' metric(s) to evaluate. Must be one or more of: \code{"APE"}, \code{"MAE"},
#' \code{"MAPE"}, \code{"RMSE"}, \code{"rRMSE1"}, \code{"rRMSE2"}. Default is
#' \code{"rRMSE2"}.
#' @param hybrid.base Logical. If \code{TRUE}, enables hybrid evaluation mode,
#' where model performance is assessed using mixed parameter combinations across
#' methods (e.g., combining \code{ka} from NCA, \code{cl} from graphical method,
#' and \code{vd} from steady-state). If \code{FALSE}, evaluation is performed
#' separately for each parameter estimation method (e.g., using only \code{simpcal},
#' only \code{nca}, etc.). This allows more flexible exploration of parameter set
#' synergies across estimation strategies. Default is \code{TRUE}.
#'
#' @return A named list combining all control modules.
#' @export
#' @examples
#' initsControl(
#'   pooled_control = pooled_control(nbins = 8),
#'   fallback_control = fallback_control(sigma_method_additive = "fixed_fraction")
#' )
#'
#' @seealso \code{\link{ss_control}}, \code{\link{pooled_control}}, \code{\link{nca_control}}, \code{\link{fallback_control}}}
#'
initsControl <- function(ss.control = ss_control(),
                         pooled.control = pooled_control(),
                         nca.control = nca_control(),
                         fallback.control = fallback_control(),
                         selmetrics= "rRMSE2",
                         hybrid.base= TRUE) {
  list(
    ss.control = ss.control,
    pooled.control = pooled.control,
    nca.control = nca.control,
    fallback.control = fallback.control,
    selmetrics= selmetrics,
    hybrid.base =hybrid.base
  )
}

