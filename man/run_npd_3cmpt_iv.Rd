% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_modelfit.R
\name{run_npd_3cmpt_iv}
\alias{run_npd_3cmpt_iv}
\title{Conduct and evaluate three-compartment modelling}
\usage{
run_npd_3cmpt_iv(dat, est.method, input.cl, input.vd)
}
\arguments{
\item{dat}{A data frame containing the intravenous pharmacokinetic data.}

\item{est.method}{The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".}

\item{input.cl}{The initial estimate for clearance.}

\item{input.vd}{The initial estimate for volume of distribution.}
}
\value{
A list containing the results of the model fitting, including the fitted parameters (\code{npd.3cmpt_results}), APE (\code{npd.3cmpt.APE}), MAPE (\code{npd.3cmpt.MAPE}), and the model fitting list (\code{npd.list.3cmpt}).
}
\description{
Perform parameter estimation with naive pooled data approach analysis for a three-compartment intravenous model, and use provided clearance and volume of distribution estimates to fit the model and calculate the absolute prediction error (APE) and mean absolute prediction error (MAPE).
}
\examples{
dat <- Bolus_2CPT
run_npd_3cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vd=70)
}
