#' Evaluate model performance for one-compartment PK model
#'
#' Computes common predictive performance metrics (e.g., APE, MAE, MAPEM, RMSE, rRMSE)
#' for a pharmacokinetic (PK) model based on simulated and observed concentration data.
#'
#' @param dat A data.frame containing PK data with columns such as `EVID` and `DV`.
#' @param est.method Estimation method passed to the fitting function.
#' Defaults to using `rxSolve` for model simulation and parameter estimation.
#' Alternative methods can be specified as needed (e.g., "nls", "nlm","nlminb").
#' @param ka Absorption rate constant (only used for oral route).
#' @param cl Clearance value.
#' @param vd Volume of distribution.
#' @param route A character string indicating the route of administration.
#' Must be one of `"oral"`, `"infusion"`, or `"bolus"`. Defaults to `"bolus"`.
#'
#' @return A numeric vector of length 5 containing the following performance metrics (rounded to 3 decimal places):
#' \describe{
#'   \item{APE}{Absolute prediction error}
#'   \item{MAE}{Mean absolute error}
#'   \item{MAPE}{Mean absolute percentage error}
#'   \item{RMSE}{Root mean square error}
#'   \item{rRMSE}{Relative RMSE}
#' }
#' If simulation fails or input parameters are missing, returns a vector of NAs.
#'
#' @details
#' Internally selects the appropriate one-compartment model fitting function depending on the route of administration:
#' \itemize{
#'   \item Uses `Fit_1cmpt_oral()` for oral data
#'   \item Uses `Fit_1cmpt_iv()` for intravenous data
#' }
#'
#' @seealso \code{\link{Fit_1cmpt_oral}}, \code{\link{Fit_1cmpt_iv}}, \code{\link{metrics.}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' eval_perf_1cmpt(
#'   dat = Oral_1CPT,
#'   est.method = "rxSolve",
#'   ka = 1,
#'   cl = 4,
#'   vd = 70,
#'   route = "oral"
#' )
#' }
eval_perf_1cmpt <- function(dat,
                            est.method = "rxSolve",
                            ka = NULL,
                            cl = NULL,
                            vd = NULL,
                            route = c("bolus", "infusion", "oral")) {
  # safe check
  route <- tryCatch(
    match.arg(route, choices = c("bolus", "oral", "infusion")),
    error = function(e) {
      stop(
        sprintf(
          "Invalid value for `%s`: '%s'. Must be one of: %s.",
          "route",
          as.character(route),
          paste(shQuote(c(
            "bolus", "oral", "infusion"
          )), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  )

  # Defensive checks
  if (any(is.na(c(cl, vd))) || cl <= 0 || vd <= 0 ||
      (route == "oral" && (is.na(ka) || ka <= 0))) {
    return(rep(NA, 5))
  }

  model_func <- if (route == "oral")
    Fit_1cmpt_oral
  else
    Fit_1cmpt_iv
  if (any(is.na(c(cl, vd))) ||
      (route == "oral" && is.na(ka)))
    return(rep(NA, 5))
  sim <- tryCatch({
    model_func(
      data = dat[dat$EVID != 2, ],
      est.method = est.method,
      input.ka = if (route == "oral")
        ka
      else
        NA,
      input.cl = cl,
      input.vd = vd,
      input.add = 0
    )
  }, error = function(e)
    return(NULL))

  if (is.null(sim))
    return(rep(NA, 5))

  obs <- dat[dat$EVID == 0, ]$DV
  pred <- sim$cp

  metrics_vec <- metrics.(pred.x = pred, obs.y = obs)
  names(metrics_vec) <- c("APE", "MAE", "MAPE", "RMSE", "rRMSE")

  rm(sim)
  gc()

  return(round(metrics_vec, 3))
}


#' Calculate metrics for model predictive performance evaluation
#'
#' calculates several error metrics, including Absolute Percentage Error (APE),
#' Mean Absolute Error (MAE), Mean Absolute Percentage Error (MAPE), Root Mean Squared Error (RMSE),
#' and Relative Root Mean Squared Error (rRMSE) based on predicted and observed values.
#'
#' @param pred.x A numeric vector representing the predicted values or simulations.
#' @param obs.y A numeric vector representing the observed values.
#'
#' @details This function will stop with an error message if the length of `pred.x` and `obs.y` are not the same.
#' It computes the following metrics:
#' \itemize{
#'   \item \strong{APE}: Absolute Predicted Error (Sum of absolute differences between predicted and observed values)
#'   \item \strong{MAE}: Mean Absolute Error (Average absolute difference between predicted and observed values)
#'   \item \strong{MAPE}: Mean Absolute Percentage Error (Sum of absolute percentage differences)
#'   \item \strong{RMSE}: Root Mean Squared Error (Square root of the average squared differences between predicted and observed values)
#'   \item \strong{rRMSE}: Relative Root Mean Squared Error (Mean of relative squared differences scaled by the average of predicted and observed values)
#' }
#'
#' @return A numeric vector with the computed values of the following metrics:
#' \enumerate{
#'   \item APE (Absolute Predicted Error)
#'   \item MAE (Mean Absolute Error)
#'   \item MAPE (Mean Absolute Percentage Error)
#'   \item RMSE (Root Mean Squared Error)
#'   \item rRMSE (Relative Root Mean Squared Error)
#' }
#'
#' @examples
#' \dontrun{
#'   ## Define the one-compartment model
#'   one.compartment <- function() {
#'     ini({
#'       tka <- log(1.57); label("Ka")
#'       tcl <- log(2.72); label("Cl")
#'       tv <- log(31.5); label("V")
#'       eta.ka ~ 0.6
#'       eta.cl ~ 0.3
#'       eta.v ~ 0.1
#'       add.sd <- 0.7
#'     })
#'     # Model block with error and model specification
#'     model({
#'       ka <- exp(tka + eta.ka)
#'       cl <- exp(tcl + eta.cl)
#'       v <- exp(tv + eta.v)
#'       d/dt(depot) <- -ka * depot
#'       d/dt(center) <- ka * depot - cl / v * center
#'       cp <- center / v
#'       cp ~ add(add.sd)
#'     })
#'   }
#'
#'   ## Perform the model fitting using nlmixr2
#'   fit <- nlmixr2(one.compartment, theo_sd, est="saem", saemControl(print=0))
#'
#'   ## Calculate error metrics using the metrics function
#'   metrics.(pred.x = fit$PRED, obs.y = fit$DV)
#' }
#'
#' @export

metrics. <-function(pred.x,
                    obs.y){

  # Check.length of pred.x and length of y whether same
  if (length(pred.x)!=length(obs.y)){
    stop("Error, number of simuations are different with predictions")
  }
  # absolute percentage error  (APE)
  metrics.ape <- sum(abs(pred.x - obs.y),na.rm = T)

  # mean absolute error (MAE)
  metrics.mae <- mean(abs(pred.x - obs.y),na.rm = T)

  # mean absolute percentage error  (MAPE)
  metrics.mape <-mean(abs(pred.x - obs.y)/obs.y,na.rm = T)*100

  # Root mean squared error (RMSE)
  metrics.rmse <- sqrt(mean((pred.x - obs.y)^2,na.rm = T))

  # relative root mean squared error (rRMSE)
  metrics.rrmse<- sqrt(mean(((pred.x - obs.y)^2/ ((pred.x + obs.y)/2)^2),na.rm = T))*100

  return(c( metrics.ape=metrics.ape ,
            metrics.mae=metrics.mae,
            metrics.mape=metrics.mape ,
            metrics.rmse=  metrics.rmse,
            metrics.rrmse= metrics.rrmse))

}
