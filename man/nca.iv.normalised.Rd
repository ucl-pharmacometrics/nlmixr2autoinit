% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_nca.normalised.R
\name{nca.iv.normalised}
\alias{nca.iv.normalised}
\title{Non-compartmental analysis for intravenous pharmacokinetic data}
\usage{
nca.iv.normalised(dat, trap.rule.method, nlastpoints)
}
\arguments{
\item{dat}{A data frame with at least two columns: TIME and DV.}

\item{nlastpoints}{Number of last points to use for the linear regression on terminal slope (default is 4).}
}
\value{
A named vector containing the calculated clearance (cl), volume of distribution (vd), slope, and half-life.
}
\description{
Perform non-compartmental analysis (NCA) on intravenous data to calculate pharmacokinetic parameters (clearance and volume of distribution) from the provided data.
}
\examples{
dat <- data.frame(TIME = c(0.5, 1, 2, 4, 6, 8, 10), DV = c(12, 8, 5, 3, 2, 1.5, 1 ))
calc.nca.iv.normalised(dat, nlastpoints = 4)
}
