#' Bin PK time-concentration data using quantile or algorithmic binning
#'
#' Bins pharmacokinetic (PK) data by time, using either equal-frequency (quantile)
#' binning or algorithmic methods supported by the `vpc` package, such as Jenks
#' natural breaks or k-means clustering.
#'
#' This function supports both quantile-based binning (used in NCA) and
#' adaptive binning strategies (used in VPC), with optional automatic
#' bin count determination.
#'
#' @param dat A data frame containing PK data. Must include:
#'   - `tad`: time after dose
#'   - `DVstd`: standardized concentration (e.g., DV/dose)
#'   - `EVID`: (optional) event ID; rows with `EVID == 0` will be used
#' @param nbins Number of bins or `"auto"`:
#'   - For `method = "quantile"`: `"auto"` defaults to 10 (deciles)
#'   - For other methods: `"auto"` uses Sturges' rule: `ceiling(log2(n)) + 1`
#' @param method Binning strategy (default = `"quantile"`). Options:
#'   - `"quantile"`: Equal-frequency binning by empirical quantiles
#'   - `"jenks"`: Natural breaks (minimizes intra-bin variance)
#'   - `"kmeans"`, `"pretty"`, `"sd"`, `"equal"`, `"density"`: other methods from `classInt`
#'
#' @return A list containing:
#' \describe{
#'   \item{binned.df}{Data frame of median time and concentration per bin}
#'   \item{bin_limits.df}{Data frame with upper time boundary for each bin}
#'   \item{breaks}{Breakpoints used to define time bins}
#'   \item{method_used}{Binning method used}
#'   \item{nbins_final}{Final number of bins (length of breaks - 1)}
#' }
#'
#' @examples
#'
# dat <- Bolus_1CPT
# # prepare data
# dat <- nmpkconvert(dat)
# dat <- calculate_tad(dat)
# dat$DVstd <- dat$DV / dat$dose
#'
#' # Default: Equal-frequency binning into deciles
#' bin.time(dat)
#'
#' \dontrun{
#' # Jenks binning with auto-selected bin count
#' bin.time(dat, method = "jenks")
#' }
#'
#' @importFrom stats quantile aggregate
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
  # get the default [1] if no input
  bin.method <- match.arg(bin.method)

  df <- subset(dat, EVID == 0, select = c(tad, DVstd))
  names(df) <- c("Time", "Conc")
  df <- df[order(df$Time), ]

  # 1. Determine bin breaks
  if (bin.method == "quantile") {
    if (identical(nbins, "auto"))
      nbins <- 10
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
      Lower = head(brks, -1),
      Upper = tail(brks, -1)
    )

  list(
    binned.df        = data.frame(Time = med.time$Time, Conc = med.conc$Conc),
    bin_limits.df    = bin_limits.df,
    breaks           = brks,
    method_used      = bin.method,
    nbins_final      = length(brks) - 1
  )
}

