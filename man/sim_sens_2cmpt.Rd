% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_sens_params.R
\name{sim_sens_2cmpt}
\alias{sim_sens_2cmpt}
\title{Simulation-based sensitivity analysis for  for a two-compartment model}
\usage{
sim_sens_2cmpt(
  dat,
  sim_vc_list,
  sim_vp_list = NA,
  sim_q_list = NA,
  estcl,
  estka = NA,
  noniv_flag = 0
)
}
\arguments{
\item{dat}{A data frame containing the pharmacokinetic data}

\item{sim_vc_list}{A list of central compartment volumes (Vc) to simulate (optional).}

\item{sim_vp_list}{A list of peripheral compartment volumes (Vp) to simulate (optional).}

\item{sim_q_list}{A list of inter-compartmental clearance (Q) to simulate (optional).}

\item{estcl}{Estimated clearance.}

\item{estka}{Estimated absorption rate if oral case.}

\item{estvc}{Estimated central volume of distribution (optional, required if \code{sim_vc_list} and \code{sim_vp_list} are not provided).}
}
\value{
A data frame containing the Vc, Vp, APE, MAPE, and time spent for each simulation.
}
\description{
Performs sensitivity analysis by testing a series of potential ratios of vc to vp.
}
\examples{

# Example 1. IV case
dat <- Bolus_2CPT
head(sim_sens_2cmpt(dat, estcl = 4, sim_vc_list = 70))

# Example 2. Oral case
dat <- Oral_2CPT
head(sim_sens_2cmpt(dat, estka=1, estcl = 4, sim_vc_list = 70))

}
