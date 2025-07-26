#' Convert pharmacokinetic dataset to depreciated format with additional dosing
#'
#' Process pharmacokinetic (PK) data to handle additional dosing (ADDL) by expanding the dose records accordingly.
#' @param dat A data frame containing the pharmacokinetic data with columns for EVID, ADDL, TIME, II, AMT, and other relevant variables.
#' @return A data frame with expanded dosing records and sorted by ID and TIME.
#' @examples
#' dat <- data.frame(ID = rep(1, 6), EVID = c(1, 0, 0, 1, 0, 0), ADDL = c(2, 0, 0, 0, 0, 0), TIME = c(0, 1, 2, 3, 4, 5), II = c(24, 0, 0, 0, 0, 0), AMT = c(100, 0, 0, 0, 0, 0))
#' nmpkconvert(dat)
#' @export
nmpkconvert <- function(dat) {
  # Extract the dose and concentration
  dat.dose <- dat[dat$EVID %in% c(1, 101, 4),]
  dat.conc <- dat[dat$EVID == 0,]
  rownames(dat.dose) <- suppressWarnings(suppressMessages(seq(1, nrow(dat.dose), 1)))
  # Identify rows with ADDL > 0, output as logical vector
  addl_rows <- dat.dose$ADDL > 0
  # Expand rows where ADDL > 0
  expanded_dose <- do.call(rbind,
                           lapply(which(addl_rows), function(di) {
                             df <- dat.dose[di, , drop = FALSE]
                             new.addl <-
                               df[rep(seq_len(nrow(df)), each = df$ADDL),]
                             new.addl$TIME <-
                               seq(df$TIME + df$II, df$TIME + df$ADDL * df$II, by = df$II)
                             return(new.addl)
                           }))
  # Combine original and expanded dose data
  dat.dose <- rbind(dat.dose, expanded_dose)
  # Combine concentration and dose data, and sort
  dat.new <- rbind(dat.conc, dat.dose)
  dat.new <- dat.new[with(dat.new, order(ID, TIME,-AMT)),]
  dat.new$ADDL <- 0
  dat.new$II <- 0
  return(dat.new)
}
