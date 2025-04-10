% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPPKinits.R
\name{getPPKinits}
\alias{getPPKinits}
\title{Get initial estimates for a population pharmacokinetic modelling}
\usage{
getPPKinits(dat, run.option = 0, getinitsControl = initsControl())
}
\arguments{
\item{dat}{A data frame containing the pharmacokinetic data.}

\item{run.option}{Integer value indicating which methods to use. It has three possible values:
\itemize{
\item \code{0}: Default pipeline calculation, not involving parameter estimation or model fitting.
\item \code{1}: Applies both the pipeline calculation methods and naive pooled data compartmental analysis.
\item \code{2}: Uses only naive pooled data compartmental analysis.
}}

\item{getinit.settings}{A list or data frame containing calculation settings (optional). The following settings can be provided:
\itemize{
\item \code{half_life}: Numeric value for the half-life of the drug (default is \code{NA}). If not provided, it will be estimated based on the data.
\item \code{nlastpoints}: Numeric value specifying the number of last data points used for linear regression to obtain the slope in the terminal phase. (default is \code{4}).
\item \code{trapezoidal.rule}: Numeric value indicating the trapezoidal rule method to use (default is \code{0}).
\code{1} refers to the linear trapezoidal method, while \code{2} refers to the linear-up/log-down trapezoidal method.
\item \code{nbins}: Numeric value specifying the number of time windows for quantile-based partitioning of the time variable in the dataset when performing naive pooling of data (default is \code{8}).
\item \code{est.method}: Character string indicating the estimation method for naive pooled data compartmental analysis to use (default is \code{"nls"}) or other methods ((e.g., "nls", "nlm", "nlminb","foce", "focei") depending on the analysis.
\item \code{selection.criteria}: Character string indicating the selection criteria for method comparison (default is \code{"MAPE"}). The selection method can be set to either Mean Absolute Percentage Error (\code{"MAPE"}) or Absolute Percentage Error (\code{"APE"}), depending on the user requirement.
\item \code{npdcmpt.inits.strategy}: Numeric value indicating the strategy for setting initial estimates in the model (default is \code{0}).  \code{0} means all initial estimates are set to \code{1} in this step, while \code{1} means that the initial estimates are based on parameters derived from non-model-fitting calculation methods.

}
If any of these settings are not provided, default values will be used.}
}
\value{
A list containing data information, initial parameter estimates, messages, and run history.
}
\description{
Computes initial values of pharmacokinetic parameters using integrated pipeline which includes single-point method, non-compartmental analysis, graphical methods, simulation-based analysis and parameter estimation with naive pooled data approach compartmental analysis if specified.
}
\examples{
inits.out<-getPPKinits(dat = Infusion_1CPT,run.option = 2,getinitsControl = initsControl(est.method = "nls"))
inits.out




}
