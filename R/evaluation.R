#' Calculate metrics for model predictive performance evaluation
#'
#' Computes common error metrics that quantify the predictive performance
#' of pharmacometric models by comparing predicted (pred.x) and observed (obs.y)
#' concentration values.
#'
#' @param pred.x Numeric vector of model-predicted values.
#' @param obs.y Numeric vector of corresponding observed values.
#'
#' @details
#' The function stops with an error if pred.x and obs.y have unequal lengths.
#' The following metrics are calculated:
#'
#' \deqn{APE = \sum |pred.x - obs.y|}
#' Absolute prediction error (APE) is the sum of absolute differences.
#'
#' \deqn{MAE = \frac{1}{n} \sum |pred.x - obs.y|}
#' Mean absolute error (MAE) expresses the average absolute deviation.
#'
#' \deqn{MAPE = \frac{100}{n} \sum \left| \frac{pred.x - obs.y}{obs.y} \right|}
#' Mean absolute percentage error (MAPE) normalizes the error by observed values.
#'
#' \deqn{RMSE = \sqrt{\frac{1}{n} \sum (pred.x - obs.y)^2}}
#' Root mean squared error (RMSE) penalizes larger deviations.
#'
#' \deqn{rRMSE1 = \frac{RMSE}{\bar{obs.y}} \times 100}
#' Relative RMSE type 1 is the RMSE normalized by the mean observed value.
#'
#' \deqn{rRMSE2 = 100 \times \sqrt{\frac{1}{n} \sum \left(
#' \frac{pred.x - obs.y}{(pred.x + obs.y)/2} \right)^2}}
#' Relative RMSE type 2 is symmetric and normalizes by the mean of each
#' predicted–observed pair.
#'
#' @return
#' A numeric vector with named elements:
#'   - APE: absolute prediction error
#'   - MAE: mean absolute error
#'   - MAPE: mean absolute percentage error
#'   - RMSE: root mean squared error
#'   - rRMSE1: relative RMSE (type 1)
#'   - rRMSE2: relative RMSE (type 2)
#'
#' @examples
#'
#' obs.y  <- rnorm(100, mean = 100, sd = 10)
#' pred.x <- obs.y + rnorm(100, mean = 0, sd = 5)
#' metrics.(pred.x = pred.x, obs.y = obs.y)
#'
#' @export

metrics. <- function(pred.x,
                     obs.y) {
  # Check.length of pred.x and length of y whether same
  if (length(pred.x) != length(obs.y)) {
    stop("Error, number of simuations are different with predictions")
  }
  # Absolute Predicted Error (APE)
  metrics.ape <- sum(abs(pred.x - obs.y), na.rm = TRUE)

  # Mean Absolute Error (MAE)
  metrics.mae <- mean(abs(pred.x - obs.y), na.rm = TRUE)

  # mean absolute percentage error  (MAPE)
  metrics.mape <- mean(abs(pred.x - obs.y) / obs.y, na.rm = TRUE) * 100

  # Root mean squared error (RMSE)
  metrics.rmse <- sqrt(mean((pred.x - obs.y) ^ 2, na.rm = TRUE))

  # relative root mean squared error (rRMSE) (type1, traditional Relative RMSE)
  metrics.rrmse1 <-
    sqrt(mean((pred.x - obs.y) ^ 2, na.rm = TRUE)) / mean(obs.y, na.rm = TRUE) * 100

  # relative root mean squared error (rRMSE) (type2, Symmetric Relative RMSE)
  metrics.rrmse2 <-
    sqrt(mean(((pred.x - obs.y) ^ 2 / ((pred.x + obs.y) / 2) ^ 2), na.rm = TRUE)) * 100

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


#' Evaluates predictive performance of a one-compartment model
#'
#' Computes predictive error metrics by comparing simulated and observed
#' concentration–time data using specified pharmacokinetic parameters and dosing route.

#' @param dat A data frame containing raw time–concentration data in the
#'   standard nlmixr2 format.
#' @param est.method Estimation method passed to the fitting function.
#' Defaults to using `rxSolve` for model simulation and parameter estimation.
#' @param ka Absorption rate constant.
#' @param cl Clearance value.
#' @param vd Volume of distribution.
#' @param route A character string indicating the route of administration.
#' Must be one of `"oral"`, `"infusion"`, or `"bolus"`. Defaults to `"bolus"`.
#'
#' @details
#' Internally selects the appropriate one-compartment model fitting function, using
#' `Fit_1cmpt_oral()` for oral administration and `Fit_1cmpt_iv()` for intravenous administration.
#' Predictive performance is quantified using the `metrics.()` function.
#'
#' @return A numeric vector containing absolute prediction error, mean absolute error,
#' mean absolute percentage error, root mean square error, and relative root mean
#' square error.
#'
#' @seealso \code{\link{Fit_1cmpt_oral}}, \code{\link{Fit_1cmpt_iv}}, \code{\link{metrics.}}
#'
#' @export
#'
#' @examples
#'
#' eval_perf_1cmpt(
#'   dat = Oral_1CPT,
#'   est.method = "rxSolve",
#'   ka = 1,
#'   cl = 4,
#'   vd = 70,
#'   route = "oral"
#' )
#'
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

  # safely release the simulation object
  sim<-NULL

  return(round(metrics_vec, 3))
}


#' Generate Unique Mixture Parameter Grid (with Deduplication and NA Removal)
#'
#' Constructs a grid of all combinations of ka, cl, and vd parameters from
#' different sources (e.g., simpcal, graph, NCA methods), and removes combinations
#' that are redundant based on relative tolerance. Only oral routes consider ka value
#' for deduplication. Any parameter value set that includes NA is removed up front.
#'
#' @param route Route of administration. Must be one of bolus, oral, or infusion.
#' @param dat A data.frame containing PK data with columns such as EVID and DV.
#' @param sp_out_ka Numeric; ka estimated from adaptive single-point methods.
#' @param sp_out_cl Numeric; clearance estimated from adaptive single-point methods.
#' @param sp_out_vd Numeric; volume of distribution estimated from adaptive single-point methods.
#' @param graph_out_ka Numeric; ka estimated from naive pooled graphic methods.
#' @param graph_out_cl Numeric; clearance estimated from naive pooled graphic methods.
#' @param graph_out_vd Numeric; volume of distribution estimated from naive pooled graphic methods.
#' @param nca_fd_ka Numeric; ka estimated from naive pooled NCA using first-dose data.
#' @param nca_fd_cl Numeric; clearance estimated from naive pooled NCA using first-dose data.
#' @param nca_fd_vd Numeric; volume of distribution estimated from naive pooled NCA using first-dose data.
#' @param nca_efd_ka Numeric; ka estimated from naive pooled NCA using repeated-dose data.
#' @param nca_efd_cl Numeric; clearance estimated from naive pooled NCA using repeated-dose data.
#' @param nca_efd_vd Numeric; volume of distribution estimated from naive pooled NCA using repeated-dose data.
#' @param nca_all_ka Numeric; ka estimated from naive pooled NCA using combined first- and repeated-dose data.
#' @param nca_all_cl Numeric; clearance estimated from naive pooled NCA using combined first- and repeated-dose data.
#' @param nca_all_vd Numeric; volume of distribution estimated from naive pooled NCA using combined first- and repeated-dose data.
#' @param verbose Logical; if TRUE (default), displays a textual progress bar during model
#'   evaluation using the 'progressr' package. Set to FALSE to run
#'   silently without showing progress information, which is recommended for
#'   automated analyses or CRAN checks.
#'
#' @return A `data.frame` of unique parameter combinations with source labels and values.
#' @examples
#' dat <- Bolus_1CPT
#' # Example parameter estimates from different methods
#' sp_out_ka <- 1.2; sp_out_cl <- 3.5; sp_out_vd <- 50
#' graph_out_ka <- 1.1; graph_out_cl <- 3.6; graph_out_vd <- 52
#' nca_fd_ka <- 1.3; nca_fd_cl <- 3.4; nca_fd_vd <- 49
#' nca_efd_ka <- NA;  nca_efd_cl <- NA;  nca_efd_vd <- NA
#' nca_all_ka <- 1.25; nca_all_cl <- 3.55; nca_all_vd <- 51
#' # Run hybrid evaluation (silent)
#'  hybrid_eval_perf_1cmpt(
#'   route = "oral",
#'   dat = dat,
#'   sp_out_ka = sp_out_ka, sp_out_cl = sp_out_cl, sp_out_vd = sp_out_vd,
#'   graph_out_ka = graph_out_ka, graph_out_cl = graph_out_cl, graph_out_vd = graph_out_vd,
#'   nca_fd_ka = nca_fd_ka, nca_fd_cl = nca_fd_cl, nca_fd_vd = nca_fd_vd,
#'   nca_efd_ka = nca_efd_ka, nca_efd_cl = nca_efd_cl, nca_efd_vd = nca_efd_vd,
#'   nca_all_ka = nca_all_ka, nca_all_cl = nca_all_cl, nca_all_vd = nca_all_vd,
#'   verbose = FALSE
#' )
#' @export

hybrid_eval_perf_1cmpt <- function(route = "bolus",
                                   dat,
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
                                   nca_all_vd,
                                   verbose=TRUE) {
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

  within_tol <- function(x1, x2, tol = 0.2) {
    all(abs((x1 - x2) / x2) <= tol, na.rm = TRUE)
  }

  if (route == "oral") {
    valid_ka_sources <- names(ka_values)[!is.na(ka_values)]

    param_grid <- expand.grid(
      ka_source = valid_ka_sources,
      cl_source = valid_cl_sources,
      vd_source = valid_vd_sources,
      stringsAsFactors = FALSE
    )

  # Separate base and hybrid
    base_combos <-
      param_grid[param_grid$ka_source == param_grid$cl_source &
                   param_grid$cl_source == param_grid$vd_source,]

    hybrid_combos <-
      param_grid[!(
        param_grid$ka_source == param_grid$cl_source &
          param_grid$cl_source == param_grid$vd_source
      ), ]

    # Assign parameter values
    base_combos$ka_value <- ka_values[base_combos$ka_source]
    base_combos$cl_value <- cl_values[base_combos$cl_source]
    base_combos$vd_value <- vd_values[base_combos$vd_source]

    hybrid_combos$ka_value <- ka_values[hybrid_combos$ka_source]
    hybrid_combos$cl_value <- cl_values[hybrid_combos$cl_source]
    hybrid_combos$vd_value <- vd_values[hybrid_combos$vd_source]

    compare_cols <- c("ka_value", "cl_value", "vd_value")

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

    compare_cols <- c("cl_value", "vd_value")
  }


  # Deduplication
  # Step 1: Filter hybrid combinations similar to base combinations
  keep_base <- rep(TRUE, nrow(hybrid_combos))
  if (nrow(base_combos) > 0 && nrow(hybrid_combos) > 0) {
    for (i in seq_len(nrow(hybrid_combos))) {
      hybrid_params <- as.numeric(hybrid_combos[i, compare_cols])

      for (j in seq_len(nrow(base_combos))) {
        base_params <- as.numeric(base_combos[j, compare_cols])
        #relative error comparison (tolerance: 20%)
        rel_diff <- abs((hybrid_params - base_params) / base_params)
        if (all(rel_diff <= 0.2, na.rm = TRUE)) {
          keep_base[i] <- FALSE
          break
        }
      }
    }
  }
  hybrid_filtered <- hybrid_combos[keep_base, ]

  # Step 2: Remove similar combinations among the remaining hybrids
  if (nrow(hybrid_filtered) > 1) {
    keep_internal <- rep(TRUE, nrow(hybrid_filtered))
    param_matrix <- as.matrix(hybrid_filtered[, compare_cols])
    for (i in 1:(nrow(param_matrix) - 1)) {
      if (!keep_internal[i]) next

      for (j in (i + 1):nrow(param_matrix)) {
        if (!keep_internal[j]) next
        #relative error comparison (tolerance: 20%)
        rel_diff <- abs((param_matrix[i, ] - param_matrix[j, ]) / param_matrix[j, ])
        if (all(rel_diff <= 0.2, na.rm = TRUE)) {
          keep_internal[j] <- FALSE
        }
      }
    }
    hybrid_filtered <- hybrid_filtered[keep_internal, ]
  }

  param_grid_unique <- rbind(base_combos, hybrid_filtered)

  col_order <- c("ka_source", "cl_source", "vd_source",
                   setdiff(names(param_grid_unique),
                           c("ka_source", "cl_source", "vd_source")))
  param_grid_unique <- param_grid_unique[, col_order]

  if (isTRUE(verbose)) {
    handler.lists <- list(
      progressr::handler_progress(
        format = ":message [:bar] :percent (:current/:total)",
        width = 80
      )
    )
  } else {
    handler.lists <- list(progressr::handler_void())
  }

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
  },handlers =handler.lists  )

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
