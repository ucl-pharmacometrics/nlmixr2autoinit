#' Create time bins for pooling a pharmacokinetic data
#'
#' Creates time bins for the provided pharmacokinetic data and calculates the median time and concentration for each bin.
#' @param testdat A data frame containing the pharmacokinetic data with columns for time after dose (tad) and normalised concentration (DVnor).
#' @param nbins The number of bins to create.
#' @return A list containing two data frames: `test.pool.normalised` with median time and concentration for each bin, and `test.pool.group` with the maximum time in each bin.
#' @import nlmixr2
#' @examples
#' dat <- Bolus_1CPT
#' dat <- nmpkconvert(dat)
#' dat <- calculate_tad(dat)
#' dat$DVnor<-dat$DV/dat$dose
#' pk.time.binning(dat,8)
#' @export
#'
pk.time.binning <- function(testdat,
                            nbins) {

  test.pool <- testdat[testdat$EVID == 0, ]
  test.pool <- subset(test.pool, select = c(tad, DVnor))

  colnames(test.pool) <- c("Time", "Conc")
  test.pool <- test.pool[order(test.pool$Time), ]

  # Define the number of bins (current 8)
  times.quantiles <- quantile(sort(unique(test.pool$Time)),
                              probs = seq(0, 1,
                                          length.out = nbins))
  # print(times.quantiles)
  adjusted_times.quantiles <- sort(unique(times.quantiles))

  # Group the plasma concentration based on the calculated quantiles
  test.pool$group <- cut(
    test.pool$Time,
    breaks = adjusted_times.quantiles,
    include.lowest = TRUE,
    labels = FALSE
  )

  # Draw vertical split lines for time quantiles
  # find the max time in each groups
  test.pool.group <- aggregate(test.pool$Time,
                               list(test.pool$group),
                               FUN = max,
                               na.rm = T)

  # Calculate the median time for each group
  test.pool.subject.time <- aggregate(test.pool$Time,
                                      list(test.pool$group),
                                      FUN = median,
                                      na.rm = T)
  # Calculate the median conc for each group
  test.pool.subject.conc <- aggregate(test.pool$Conc,
                                      list(test.pool$group),
                                      FUN = median,
                                      na.rm = T)

  # Reason for collect the time and concentration separately, in case collect the extreme concentration.

  test.pool.normalised <-
    data.frame(time = test.pool.subject.time$x,
               test.pool.subject.conc$x)
  colnames(test.pool.normalised) <- c("Time", "Conc")

  colnames(test.pool.group) <- c("Group", "Upper Time Boundary ")

  return(
    list(
      test.pool.normalised = test.pool.normalised,
      test.pool.group = test.pool.group
    )
  )

}
