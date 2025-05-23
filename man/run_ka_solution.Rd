% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extra_ka.R
\name{run_ka_solution}
\alias{run_ka_solution}
\title{Calculate absorption rate constant (ka) for oral administration Data}
\usage{
run_ka_solution(df, cl, ke, Fbio = 1)
}
\arguments{
\item{df}{A data frame containing pharmacokinetic data with columns: \code{ID} (subject identifier), \code{TIME} (time point), \code{DV} (drug concentration), and \code{DOSE} (dose amount). The data frame is typically filtered to include data only from the first dosing interval (e.g., single dose and evid = 0).}

\item{cl}{A numeric value representing the clearance (CL) of the drug.}

\item{ke}{A numeric value representing the elimination rate constant (ke).}

\item{Fbio}{A numeric value representing the fraction of the drug absorbed (bioavailability, F).}
}
\value{
A list containing:
\item{ka_calc_median}{The median ka value across all individuals.}
\item{data_before_tmax}{The original data frame filtered for times before Tmax, with an added \code{ka_calc} column for individual ka values.}
}
\description{
Calculate the absorption rate constant (ka) using sampling points from the absorption phase (t<Tmax). The single-point method is applied, with each plasma concentration point used to calculate a corresponding ka value
}
\details{
The function first calculates the time to maximum concentration (Tmax) for each individual in the dataset. It then filters the data to include only the time points up to Tmax for each individual. The absorption rate constant (ka) is calculated for each time point using the \code{ka_calculation} function.
}
\examples{
# Example usage:
 df<-Oral_1CPT[Oral_1CPT$SD==1,]
 df<-processData(df)$dat
 result <- run_ka_solution(df = df, cl = 4, ke = 4/70, Fbio = 1)
 ka_median <- result[[1]]
 data_with_ka <- result[[2]]

}
