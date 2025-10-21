#' Bin time-concentration data using quantile or algorithmic binning
#'
#' Bins data by time using either equal-frequency (quantile) binning
#' or algorithmic binning methods.
#'
#' @param dat A data frame containing PK data. Must include:
#'   - tad: time after dose
#'   - DVstd: standardized concentration (DV/dose)
#'   - EVID: optional event ID column used to filter observations (EVID == 0)
#'
#' @param nbins Number of bins or "auto". If numeric with
#'   `bin.method = "quantile"`, specifies equal-frequency bins. If "auto",
#'   10 bins are used for quantile; otherwise binning is determined by
#'   `vpc::auto_bin`. Numeric nbins for non-quantile methods is passed to
#'   `vpc::auto_bin`.
#'
#' @param bin.method Binning strategy (default = "quantile").
#'   Available options are:
#'   - quantile: equal-frequency binning by empirical quantiles
#'   - jenks: natural breaks minimizing within-bin variance
#'   - kmeans, pretty, sd, equal, density: alternative binning
#'     methods from vpc::auto_bin
#'
#' @details
#' Supports quantile-based binning and other data-driven methods
#' (jenks, kmeans, pretty, sd, equal, density), with optional automatic
#' bin count selection.
#'
#' @return A list containing summary results of the time-concentration binning process.
#'
#' @seealso
#' \code{\link[vpc:auto_bin]{vpc::auto_bin}}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' dat <- Bolus_1CPT
#' dat <- nmpkconvert(dat)
#' dat <- calculate_tad(dat)
#' dat$DVstd <- dat$DV / dat$dose
#' bin.time(dat)
#'
#' @export

bin.time <- function(dat,
                     nbins   = "auto",
                     bin.method  = c("quantile",
                                     "jenks",
                                     "kmeans",
                                     "pretty",
                                     "sd",
                                     "equal",
                                     "density")) {
  # Defensive check for nbins
  if (!(is.numeric(nbins) || identical(nbins, "auto"))) {
    stop("`nbins` must be either a numeric value or the string \"auto\".",
         call. = FALSE)
  }

  if (is.numeric(nbins)) {
    if (length(nbins) != 1L || nbins <= 0 || !is.finite(nbins)) {
      stop("`nbins` must be a single positive finite number.", call. = FALSE)
    }
  }

  # get the default [1] if no input
  bin.method <- match.arg(bin.method)

  df <- subset(dat, EVID == 0, select = c(tad, DVstd))
  names(df) <- c("Time", "Conc")
  df <- df[order(df$Time), ]

  # 1. Determine bin breaks
  if (bin.method == "quantile") {
    if (identical(nbins, "auto")) {
      nbins <- 10
    }
    probs <- seq(0, 1, length.out = nbins + 1)

    brks  <-
      stats::quantile(
        unique(df$Time),
        probs = probs,
        na.rm = TRUE,
        type = 7
      )
    # Use type = 7: default quantile (linear interpolation between order statistics)
    brks  <- sort(unique(brks))
  } else {
    brks <- vpc::auto_bin(df$Time, type = bin.method, n_bins = nbins)
  }

  # 2. Assign groupings
  df$group <-
    cut(
      df$Time,
      breaks = brks,
      include.lowest = TRUE,
      labels = FALSE
    )

  # 3. Summary statistics per bin
  med.time <-
    stats::aggregate(Time ~ group, df, median, na.rm = TRUE)
  med.conc <-
    stats::aggregate(Conc ~ group, df, median, na.rm = TRUE)
  # max.time <- aggregate(Time ~ group, df, max,    na.rm = TRUE)

  bin_limits.df <-
    data.frame(
      Group = seq_along(brks[-1]),
      Lower = utils::head(brks, -1),
      Upper = utils::tail(brks, -1)
    )

  list(
    binned.df        = data.frame(Time = med.time$Time, Conc = med.conc$Conc),
    bin_limits.df    = bin_limits.df,
    breaks           = brks,
    method_used      = bin.method,
    nbins_final      = length(brks) - 1
  )
}
