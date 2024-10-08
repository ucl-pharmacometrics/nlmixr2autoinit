% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_nca.normalised.R
\name{trap.rule_linear_up_log_down}
\alias{trap.rule_linear_up_log_down}
\title{Calculate AUC using linear up and log down trapezoidal rule}
\usage{
trap.rule_linear_up_log_down(x, y)
}
\arguments{
\item{x}{A numeric vector representing the time points.}

\item{y}{A numeric vector representing the corresponding concentration at each time point.}
}
\value{
A numeric value representing the estimated AUC using the linear up/log down method.
}
\description{
This function computes the area under the curve (AUC) by using the trapezoidal rule
for phases where concentration is increasing, and a logarithmic rule for phases
where concentration is decreasing.
}
\examples{
x <- c(0.5, 1, 2, 4, 6, 8, 10)
y <- c(12, 8, 5, 3, 2, 1.5, 1)
trap.rule_linear_up_log_down(x, y)

}
