% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_hf.R
\name{get_hf}
\alias{get_hf}
\title{Estimate half-life from pooled pharmacokinetic data}
\usage{
get_hf(dat, dose_type = "first_dose", pooled = NULL, ...)
}
\arguments{
\item{dat}{A data frame containing the pharmacokinetic data. The required columns
depend on the specified \code{dose_type} and are validated within \code{\link{get_pooled_data}}.}

\item{dose_type}{Character string specifying the type of data to analyze. One of:
\itemize{
\item{\code{"first_dose"}}: Analysis based on first-dose data only.
\item{\code{"repeated_doses"}}: Analysis based on repeated-dose data only.
\item{\code{"combined_doses"}}: Combined analysis using both first and repeated doses.
}}

\item{pooled}{(Optional) A pooled data object as returned by \code{\link{get_pooled_data}}.
If not supplied, the function will internally generate pooled data using \code{dat}
and any applicable arguments to \code{bin.time}.}

\item{...}{Additional arguments passed to either \code{\link{bin.time}} (for pooling)
or \code{\link{find_best_lambdaz}} (for elimination slope calculation). The function
automatically splits these arguments based on their intended function.}
}
\value{
A named list containing:
\itemize{
\item{\code{half_life_median}}: The median of all positive half-life estimates from the subsets.
\item{\code{half_life_fd}}: The half-life estimate based on first-dose data, if available.
\item{\code{half_life_md}}: The half-life estimate based on repeated-dose data, if available.
\item{\code{half_life_all}}: The half-life estimate based on the full dataset, if available.
}
}
\description{
Estimates the terminal half-life of a drug using pooled and binned
pharmacokinetic data. It supports analysis based on first-dose data, repeated-dose data,
or a combined profile that includes both. The estimation is performed by applying
linear regression on the terminal phase of log-transformed concentration-time data.
}
\details{
The function works by:
\enumerate{
\item Generating (or using) pooled data using \code{\link{get_pooled_data}}, which applies
time binning to normalize the time-course concentration data.
\item Estimating the terminal elimination slope (\code{lambdaz}) using \code{\link{find_best_lambdaz}}.
\item Computing the terminal half-life as \code{log(2) / lambdaz} for each data subset.
\item Returning the median of all valid (positive, non-missing) half-life estimates.
}
}
\examples{

\dontrun{
# Example: half-life estimation from a combined profile (first + repeated doses)
dat <- Bolus_1CPT
dat <- processData(dat)$dat
get_hf(dat, dose_type = "combined_doses", nlastpoints = 3)

# Using externally pooled data
pooled <- get_pooled_data(dat, dose_type = "combined_doses")
get_hf(dat, dose_type = "combined_doses", pooled = pooled, nlastpoints = 2)
}

}
\seealso{
\code{\link{get_pooled_data}}, \code{\link{bin.time}}, \code{\link{find_best_lambdaz}}
}
