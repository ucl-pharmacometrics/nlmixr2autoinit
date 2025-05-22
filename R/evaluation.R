#' Calculate metrics for model predictive performance evaluation
#'
#' Calculates several error metrics, including absolute percentage error (APE),
#' mean Aabsolute error (MAE), mean absolute percentage error (MAPE), root mean squared error (RMSE),
#' and traditional relative root mean squared error (rRMSEtra) and Symmetric relative root mean squared error (rRMSEsym).
#'
#' @param pred.x A numeric vector representing the predicted values or simulations.
#' @param obs.y A numeric vector representing the observed values.
#'
#' @details This function will stop with an error message if the length of `pred.x` and `obs.y` are not the same.
#' It computes the following metrics:
#' It computes the following metrics:
#' \itemize{
#'   \item \strong{APE}: Sum of absolute differences between predicted and observed values.
#'   \item \strong{MAE}: Mean of absolute differences between predicted and observed values.
#'   \item \strong{MAPE}: Mean of absolute percentage differences between predicted and observed values.
#'   \item \strong{RMSE}: Square root of the mean of squared differences between predicted and observed values.
#'   \item \strong{rRMSEtype1}: relative RMSE, calculated as RMSE divided by the mean of observed values.
#'   \item \strong{rRMSEtype2}: relative RMSE, calculated by dividing the squared error by the square of the average of each predicted–observed pair.
#' }
#'
#' @return A numeric vector with the computed values of the following metrics:
#' \enumerate{
#'   \item APE (Absolute Predicted Error)
#'   \item MAE (Mean Absolute Error)
#'   \item MAPE (Mean Absolute Percentage Error)
#'   \item RMSE (Root Mean Squared Error)
#'   \item rRMSEtype1 (Relative RMSE, type 1)
#'   \item rRMSEtype2 ( Relative RMSE, type 2)
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' obs.y  <- rnorm(100, mean = 100, sd = 10)
#' pred.x <- obs.y + rnorm(100, mean = 0, sd = 5)
#'
#' ## Calculate error metrics
#' metrics.(pred.x = pred.x, obs.y = obs.y)
#' }
#'
#' @export

metrics. <- function(pred.x,
                     obs.y) {
  # Check.length of pred.x and length of y whether same
  if (length(pred.x) != length(obs.y)) {
    stop("Error, number of simuations are different with predictions")
  }
  # Absolute Predicted Error (APE)
  metrics.ape <- sum(abs(pred.x - obs.y), na.rm = T)

  # Mean Absolute Error (MAE)
  metrics.mae <- mean(abs(pred.x - obs.y), na.rm = T)

  # mean absolute percentage error  (MAPE)
  metrics.mape <- mean(abs(pred.x - obs.y) / obs.y, na.rm = T) * 100

  # Root mean squared error (RMSE)
  metrics.rmse <- sqrt(mean((pred.x - obs.y) ^ 2, na.rm = T))

  # relative root mean squared error (rRMSE) (type1, traditional Relative RMSE)
  metrics.rrmse1 <-
    sqrt(mean((pred.x - obs.y) ^ 2, na.rm = TRUE)) / mean(obs.y, na.rm = TRUE) * 100

  # relative root mean squared error (rRMSE) (type2, Symmetric Relative RMSE)
  metrics.rrmse2 <-
    sqrt(mean(((pred.x - obs.y) ^ 2 / ((pred.x + obs.y) / 2) ^ 2), na.rm = T)) * 100

  return(
    c(
      metrics.ape = metrics.ape ,
      metrics.mae = metrics.mae,
      metrics.mape = metrics.mape ,
      metrics.rmse =  metrics.rmse,
      metrics.rrmse1 = metrics.rrmse1,
      metrics.rrmse2 = metrics.rrmse2
    )
  )

}


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
  names(metrics_vec) <- c("APE", "MAE", "MAPE", "RMSE", "rRMSE1","rRMSE2")

  rm(sim)
  gc()

  return(round(metrics_vec, 3))
}


#' Generate Unique Mixture Parameter Grid (with Deduplication and NA Removal)
#'
#' Constructs a grid of all combinations of `ka`, `cl`, and `vd` parameters from
#' different sources (e.g., simpcal, graph, NCA methods), and removes combinations
#' that are redundant based on relative tolerance. Only oral routes consider `ka_value`
#' for deduplication. Any parameter value set that includes NA is removed up front.
#'
#' @param route Character; administration route, e.g., `"oral"` or `"iv"`.
#' @param sp_out_ka Numeric; ka from simpcal.
#' @param sp_out_cl Numeric; cl from simpcal.
#' @param sp_out_vd Numeric; vd from simpcal.
#' @param graph_out_ka Numeric; ka from graph-based method.
#' @param graph_out_cl Numeric; cl from graph.
#' @param graph_out_vd Numeric; vd from graph.
#' @param nca_fd_ka Numeric; ka from NCA FD.
#' @param nca_fd_cl Numeric; cl from NCA FD.
#' @param nca_fd_vd Numeric; vd from NCA FD.
#' @param nca_efd_ka Numeric; ka from NCA eFD.
#' @param nca_efd_cl Numeric; cl from NCA eFD.
#' @param nca_efd_vd Numeric; vd from NCA eFD.
#' @param nca_all_ka Numeric; ka from NCA All.
#' @param nca_all_cl Numeric; cl from NCA All.
#' @param nca_all_vd Numeric; vd from NCA All.
#'
#' @return A `data.frame` of unique parameter combinations with source labels and values.
#'
#' @keywords internal#
#' @export

hybrid_eval_perf_1cmpt <- function(route = "bolus",
                                   sp_out_ka,
                                   sp_out_cl,
                                   sp_out_vd,
                                   graph_out_ka,
                                   graph_out_cl,
                                   graph_out_vd,
                                   nca_fd_ka,
                                   nca_fd_cl,
                                   nca_fd_vd,
                                   nca_efd_ka,
                                   nca_efd_cl,
                                   nca_efd_vd,
                                   nca_all_ka,
                                   nca_all_cl,
                                   nca_all_vd) {
  # Source labels
  all_sources <-
    c("simpcal", "graph", "nca_fd", "nca_efd", "nca_all")

  # Precompiled parameter value vectors
  ka_values <- c(
    simpcal = sp_out_ka,
    graph = graph_out_ka,
    nca_fd = nca_fd_ka,
    nca_efd = nca_efd_ka,
    nca_all = nca_all_ka
  )

  cl_values <- c(
    simpcal = sp_out_cl,
    graph = graph_out_cl,
    nca_fd = nca_fd_cl,
    nca_efd = nca_efd_cl,
    nca_all = nca_all_cl
  )

  vd_values <- c(
    simpcal = sp_out_vd,
    graph = graph_out_vd,
    nca_fd = nca_fd_vd,
    nca_efd = nca_efd_vd,
    nca_all = nca_all_vd
  )

  valid_cl_sources <- names(cl_values)[!is.na(cl_values)]
  valid_vd_sources <- names(vd_values)[!is.na(vd_values)]

  if (route == "oral") {
    valid_ka_sources <- names(ka_values)[!is.na(ka_values)]

    param_grid <- expand.grid(
      ka_source = valid_ka_sources,
      cl_source = valid_cl_sources,
      vd_source = valid_vd_sources,
      stringsAsFactors = FALSE
    )

    # Separate base and hybrid
    base_combos <- param_grid[param_grid$ka_source == param_grid$cl_source &
                                param_grid$cl_source == param_grid$vd_source,]

    hybrid_combos <- param_grid[!(
      param_grid$ka_source == param_grid$cl_source &
        param_grid$cl_source == param_grid$vd_source
    ),]

    # Assign parameter values
    base_combos$ka_value <- ka_values[base_combos$ka_source]
    base_combos$cl_value <- cl_values[base_combos$cl_source]
    base_combos$vd_value <- vd_values[base_combos$vd_source]

    hybrid_combos$ka_value <- ka_values[hybrid_combos$ka_source]
    hybrid_combos$cl_value <- cl_values[hybrid_combos$cl_source]
    hybrid_combos$vd_value <- vd_values[hybrid_combos$vd_source]

    # Deduplication on hybrid combos
    keep_idx <- rep(TRUE, nrow(hybrid_combos))
    numeric_matrix <-
      as.matrix(hybrid_combos[, c("ka_value", "cl_value", "vd_value")])

    within_tol <- function(x1, x2, tol = 0.2) {
      all(abs((x1 - x2) / x2) <= tol)
    }

    for (i in 1:(nrow(numeric_matrix) - 1)) {
      if (!keep_idx[i])
        next
      for (j in (i + 1):nrow(numeric_matrix)) {
        if (!keep_idx[j])
          next
        if (within_tol(numeric_matrix[i,], numeric_matrix[j,])) {
          keep_idx[j] <- FALSE
        }
      }
    }

    param_grid_unique <-
      rbind(base_combos, hybrid_combos[keep_idx,])

  } else {
    # IV/bolus: no ka combinations
    param_grid <- expand.grid(
      cl_source = valid_cl_sources,
      vd_source = valid_vd_sources,
      stringsAsFactors = FALSE
    )
    param_grid$ka_source <- NA
    param_grid$ka_value <- NA
    param_grid$cl_value <- cl_values[param_grid$cl_source]
    param_grid$vd_value <- vd_values[param_grid$vd_source]
    # Separate base and hybrid combinations (by cl == vd)
    base_combos <- param_grid[param_grid$cl_source == param_grid$vd_source,]

    hybrid_combos <- param_grid[param_grid$cl_source != param_grid$vd_source,]

    # Deduplicate hybrid combos based on cl/vd similarity
    keep_idx <- rep(TRUE, nrow(hybrid_combos))
    numeric_matrix <-
      as.matrix(hybrid_combos[, c("cl_value", "vd_value")])

    within_tol <- function(x1, x2, tol = 0.2) {
      all(abs((x1 - x2) / x2) <= tol)
    }

    for (i in 1:(nrow(numeric_matrix) - 1)) {
      if (!keep_idx[i])
        next
      for (j in (i + 1):nrow(numeric_matrix)) {
        if (!keep_idx[j])
          next
        if (within_tol(numeric_matrix[i,], numeric_matrix[j,])) {
          keep_idx[j] <- FALSE
        }
      }
    }
    # Combine base and deduplicated hybrid
    param_grid_unique <-
      rbind(base_combos, hybrid_combos[keep_idx,])
  }

  # Reorder columns: ka_source, cl_source, vd_source, followed by other columns
  param_grid_unique <- param_grid_unique[c("ka_source", "cl_source", "vd_source",
                                           setdiff(
                                             names(param_grid_unique),
                                             c("ka_source", "cl_source", "vd_source")
                                           ))]

  progressr::handlers(
    progressr::handler_progress(format = ":message [:bar] :percent (:current/:total)", width = 80)
  )

  results <- progressr::with_progress({
    p <-
      progressr::progressor(steps = nrow(param_grid_unique))  # Initialize progress bar

    lapply(seq_len(nrow(param_grid_unique)), function(i) {
      row <- param_grid_unique[i, ]

      # Evaluate model performance
      perf <- eval_perf_1cmpt(dat,
                              "rxSolve",
                              as.numeric(row[["ka_value"]]),
                              as.numeric(row[["cl_value"]]),
                              as.numeric(row[["vd_value"]]),
                              route)

      p(sprintf("Evaluating combination %d", i))

      # Combine row info + performance results
      c(as.list(row), as.list(perf))
    })
  })

  # Convert to data.frame
  results_df <-
    do.call(rbind,
            lapply(results, as.data.frame.list, stringsAsFactors = FALSE))

  # Convert numeric columns
  numeric_cols <- c("ka_value",
                    "cl_value",
                    "vd_value",
                    grep("^performance", names(results_df), value = TRUE))

  results_df[numeric_cols] <-
    lapply(results_df[numeric_cols], as.numeric)

  return(results_df)
}
