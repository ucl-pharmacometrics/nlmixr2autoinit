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
fallback_control <- function(enable_ka_fallback = TRUE,
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
#'
#' @return A named list combining all control modules.
#' @export
#' @examples
#' initsControl(
#'   pooled_control = pooled_control(nbins = 8),
#'   fallback_control = fallback_control(sigma_method_additive = "fixed_fraction")
#' )
initsControl <- function(ss.control = ss_control(),
                         pooled.control = pooled_control(),
                         nca.control = nca_control(),
                         fallback.control = fallback_control(),
                         selection.criterion= "All") {
  list(
    ss.control = ss.control,
    pooled.control = pooled.control,
    nca.control = nca.control,
    fallback.control = fallback.control,
    selection.criterion= selection.criterion
  )
}

