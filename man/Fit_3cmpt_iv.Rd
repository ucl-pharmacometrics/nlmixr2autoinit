% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npd_modelcodes.R
\name{Fit_3cmpt_iv}
\alias{Fit_3cmpt_iv}
\title{Fit intravenous pharmacokinetic data to a three-compartment model}
\usage{
Fit_3cmpt_iv(
  data,
  est.method,
  input.cl,
  input.vc3cmpt,
  input.vp3cmpt,
  input.vp23cmpt,
  input.q3cmpt,
  input.q23cmpt,
  input.add
)
}
\arguments{
\item{data}{Intravenous pharmacokinetic data in the nlmixr2 format.}

\item{est.method}{The estimation method in nlmixr2 to use (e.g., "nls", "nlm", focei"...). The default value is "nls".}

\item{input.cl}{Initial estimate for clearance.}

\item{input.vc3cmpt}{Initial estimate for the central volume of distribution in a three-compartment model.}

\item{input.vp3cmpt}{Initial estimate for the first peripheral volume of distribution in a three-compartment model.}

\item{input.vp23cmpt}{Initial estimate for the second peripheral volume of distribution in a three-compartment model.}

\item{input.q3cmpt}{Initial estimate for first inter-compartmental clearance in a three-compartment model.}

\item{input.q23cmpt}{Initial estimate for second inter-compartmental clearance in a three-compartment model.}

\item{input.add}{Initial estimate for additive error.}
}
\value{
An object of class 'nlmixr2' containing the results of the estimation.
}
\description{
Perform parameter estimation using the naive pooled data approach with a one-compartment model for the provided intravenous pharmacokinetic data.
}
\examples{
dat <- Bolus_2CPT
Fit_3cmpt_iv(dat=dat,est.method="nls",input.cl=4,input.vc3cmpt=70,input.vp3cmpt=35,input.vp23cmpt=35,input.q3cmpt=10,input.q23cmpt=10,input.add=1)
}
