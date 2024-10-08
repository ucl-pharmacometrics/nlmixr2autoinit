% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_tad.R
\name{calculate_tad}
\alias{calculate_tad}
\title{Calculate time after dose (TAD) for a pharmacokinetic data}
\usage{
calculate_tad(dat, infusion_flag = 0)
}
\arguments{
\item{dat}{A data frame containing the pharmacokinetic data.}

\item{infusion_flag}{A flag indicating whether it is an infusion case}
}
\value{
A data frame with the calculated tad (time after the last dose).
}
\description{
Calculates the time after dose (TAD) for each observation in a pharmacokinetic dataset. The TAD is calculated as the time difference between each observation and the most recent dose administration. The function handles both bolus and infusion dosing scenarios, adjusting the calculations accordingly based on the infusion_flag parameter.
}
\examples{

Example 1:
dat <- Bolus_1CPT
calculate_tad(dat)

Example 2:
dat <- Infusion_1CPT
dat<-calculate_tad(dat,infusion_flag=1)
dat
}
