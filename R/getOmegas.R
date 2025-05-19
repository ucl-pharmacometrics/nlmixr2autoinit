#' Generate ETA Variance and Covariance Table
#'
#' This function constructs a combined table containing:
#' \itemize{
#'   \item \strong{ETA variances} (e.g., \code{eta.cl}, \code{eta.vc}), which represent inter-individual variability (IIV) in pharmacokinetic parameters;
#'   \item \strong{Derived covariances} (e.g., \code{cov.eta_cl_vc}) computed from ETA variances and assumed pairwise correlations.
#' }
#'
#' ETA variances are initialized to 0.1 by default. Correlations within defined omega blocks
#' (block 1: \code{eta.vmax}, \code{eta.km}; block 2: \code{eta.cl}, \code{eta.vc}, \code{eta.vp}, \code{eta.vp2}, \code{eta.q}, \code{eta.q2})
#' are assumed to be 0.1 and used to compute covariances as:
#' \deqn{Cov(i, j) = sqrt(Var_i) * sqrt(Var_j) * Corr(i, j)}
#'
#' The resulting output format aligns with \code{Recommended_initial_estimates}.
#'
#' @return A \code{data.frame} with columns: \code{Parameters}, \code{Methods}, and \code{Values}.
#' @export

getOmegas <- function() {
  # Define ETA variance parameters (fixed)
  eta_names <- c(
    "eta.ka",
    "eta.cl",
    "eta.vc",
    "eta.vp",
    "eta.q",
    "eta.vp2",
    "eta.q2",
    "eta.vmax",
    "eta.km"
  )

  # Assign all ETA variances to 0.1
  eta_values <- setNames(rep(0.1, length(eta_names)), eta_names)

  # Omega blocks defining which correlations to include
  omega_block1 <- c("eta.vmax", "eta.km")
  omega_block2 <-
    c("eta.cl", "eta.vc", "eta.vp", "eta.vp2", "eta.q", "eta.q2")

  # Combine all correlation pairs from the blocks
  make_combinations <-
    function(block)
      combn(block, 2, simplify = FALSE)
  correlation_pairs <-
    c(make_combinations(omega_block1),
      make_combinations(omega_block2))

  # Assign default correlation values (e.g., 0.1)
  cor_eta_values <- sapply(correlation_pairs, function(x)
    0.1)
  names(cor_eta_values) <-
    sapply(correlation_pairs, function(combo) {
      paste0("cor.eta_",
             gsub("eta.", "", combo[1]),
             "_",
             gsub("eta.", "", combo[2]))
    })

  # Create ETA variance entries (diagonal)
  eta_table <- data.frame(
    Parameters = names(eta_values),
    Methods = rep("fixed_values", length(eta_values)),
    Values = format(eta_values, scientific = FALSE, digits = 3),
    stringsAsFactors = FALSE
  )

  # Compute covariances from correlations
  cov_table <- data.frame(
    Parameters = character(0),
    Methods = character(0),
    Values = character(0),
    stringsAsFactors = FALSE
  )

  for (pair in correlation_pairs) {
    i <- pair[1]
    j <- pair[2]
    cov_name <-
      paste0("cor.eta_", gsub("eta.", "", i), "_", gsub("eta.", "", j))

    cov_value <-
      sqrt(eta_values[[i]]) * sqrt(eta_values[[j]]) * cor_eta_values[[cov_name]]

    cov_table <- rbind(
      cov_table,
      data.frame(
        Parameters = cov_name,
        Methods = "derived_from_eta_and_correlation (0.1)",
        Values = format(cov_value, scientific = FALSE, digits = 3),
        stringsAsFactors = FALSE
      )
    )
  }

  # Combine into one unified table
  full_table <- rbind(eta_table, cov_table)
  rownames(full_table) <- NULL
  return(full_table)
}



