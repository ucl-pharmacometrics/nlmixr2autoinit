% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_sens_params.R
\name{sim_sens_1cmpt_mm}
\alias{sim_sens_1cmpt_mm}
\title{Parameter sweeping for a one-compartment nonlinear (Michaelis-Menten) model}
\usage{
sim_sens_1cmpt_mm(
  dat,
  sim_vmax = list(mode = "auto", values = NULL, est.cl = NULL),
  sim_km = list(mode = "auto", values = NULL),
  sim_vd = list(mode = "manual", values = NULL),
  sim_ka = list(mode = "manual", values = NULL),
  route = c("iv", "oral")
)
}
\arguments{
\item{dat}{A data frame containing the pharmacokinetic dataset.}

\item{sim_vmax}{A list for Vmax simulation. Default: \code{list(mode = "auto", values = NULL)}.
\describe{
\item{mode}{\code{"manual"} or \code{"auto"}.}
\item{values}{Vector of candidate Vmax values if \code{"manual"}.}
\item{est.cl}{Estimated clearance, required if \code{"auto"}.}
}}

\item{sim_km}{A list for Km simulation. Default: \code{list(mode = "auto", values = NULL)}.
\describe{
\item{mode}{\code{"manual"} or \code{"auto"}.}
\item{values}{Vector of candidate Km values if \code{"manual"}.}
}}

\item{sim_vd}{A list for volume of distribution simulation. Default: \code{list(mode = "manual", values = NULL)}.
\describe{
\item{mode}{Must be \code{"manual"}.}
\item{values}{Vector of Vd values. Required.}
}}

\item{sim_ka}{A list for absorption rate constant (Ka). Used only if \code{route = "oral"}. Default: \code{list(mode = "manual", values = NULL)}.
\describe{
\item{mode}{Must be \code{"manual"}.}
\item{values}{Vector of Ka values. Required for oral.}
}}

\item{route}{A character string indicating administration route. One of \code{"iv"} or \code{"oral"}. Default is \code{"iv"}.}
}
\value{
A data frame with simulated parameter combinations and evaluation metrics:
\describe{
\item{Vmax, Km, Vd, Ka}{The parameter values used.}
\item{APE, MAE, MAPE, RMSE, rRMSE}{Model fit evaluation metrics.}
\item{Cumulative.Time.Sec}{Total elapsed time (in seconds) for each simulation.}
}
}
\description{
Performs sensitivity analysis via parameter sweeping in a one-compartment pharmacokinetic model
with Michaelis-Menten elimination. By systematically varying user-defined or automatically
generated values for key pharmacokinetic parameters (e.g., Vmax, Km, Vd, and Ka), the function
evaluates model sensitivity and fit across a defined grid of parameter combinations.
}
\examples{
\dontrun{
# Example 1: IV route
out <- sim_sens_1cmpt_mm(
  dat = Bolus_1CPTMM,
  sim_vmax = list(mode = "auto",est.cl=4),
  sim_km   = list(mode = "auto"),
  sim_vd   = list(mode = "manual", values = c(70)),
  sim_ka   = list(mode = "manual", values = NA),  # NA for IV route
  route = "iv"
)
head(out)

# Example 2: Oral route
#out <- sim_sens_1cmpt_mm(
#  dat = Oral_1CPTMM,
#  sim_vmax = list(mode = "auto",est.cl=1),
#  sim_km   = list(mode = "auto"),
#  sim_vd   = list(mode = "manual", values = c(48)),
#  sim_ka   = list(mode = "manual", values = c(1)),  # 1 for oral route
#  route = "oral"
#)
#head(out)
}

}
