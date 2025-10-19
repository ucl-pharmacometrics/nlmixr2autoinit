#' Expand additional dosing (ADDL) records for pharmacokinetic analysis
#'
#' Expands dosing records that contain additional doses (ADDL) using the
#' specified interdose interval (II). Each additional dose is converted into
#' an explicit record to provide a complete dosing history suitable for
#' population pharmacokinetic modeling.
#'
#' @param dat A data frame containing raw time–concentration data in the standard nlmixr2 format.
#'
#' @return A data frame with expanded dosing records. The columns ADDL and II
#'   are reset to zero after expansion.
#'
#' @details
#' Dosing records with ADDL greater than zero are expanded using the formula:
#' TIME_new = TIME + n × II, where n ranges from 1 up to ADDL.
#' Observation records (EVID = 0) are not modified.
#'
#' @examples
#' # Dataset with a single subject and additional dosing
#' dat <- data.frame(
#'   ID   = rep(1, 6),
#'   EVID = c(1, 0, 0, 1, 0, 0),
#'   ADDL = c(2, 0, 0, 0, 0, 0),
#'   TIME = c(0, 1, 2, 3, 4, 5),
#'   II   = c(24, 0, 0, 0, 0, 0),
#'   AMT  = c(100, 0, 0, 0, 0, 0)
#' )
#' nmpkconvert(dat)
#'
#' @export

nmpkconvert <- function(dat) {
  # Extract the dose and concentration
  dat.dose <- dat[dat$EVID %in% c(1, 101, 4), ]
  dat.conc <- dat[dat$EVID == 0, ]
  suppressWarnings(suppressMessages(rownames(dat.dose) <-
                                      seq(1, nrow(dat.dose), 1)))
  # Identify rows with ADDL > 0, output as logical vector
  addl_rows <- dat.dose$ADDL > 0
  # Expand rows where ADDL > 0
  expanded_dose <- do.call(rbind,
                           lapply(which(addl_rows), function(di) {
                             df <- dat.dose[di, , drop = FALSE]
                             new.addl <-
                               df[rep(seq_len(nrow(df)), each = df$ADDL), ]
                             new.addl$TIME <-
                               seq(df$TIME + df$II, df$TIME + df$ADDL * df$II, by = df$II)
                             return(new.addl)
                           }))
  # Combine original and expanded dose data
  dat.dose <- rbind(dat.dose, expanded_dose)
  # Combine concentration and dose data, and sort
  dat.new <- rbind(dat.conc, dat.dose)
  dat.new <- dat.new[with(dat.new, order(ID, TIME, -AMT)), ]
  dat.new$ADDL <- 0
  dat.new$II <- 0
  return(dat.new)
}
