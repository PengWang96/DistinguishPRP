% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_posterior.R
\name{fixed_effect}
\alias{fixed_effect}
\title{Fit Fixed Effect Model using Monte Carlo Algorithm}
\usage{
fixed_effect(N, m, hat_beta, hat_sigma_sq)
}
\arguments{
\item{N}{Total number of iterations for the Monte Carlo simulation.}

\item{m}{Number of replications for the test statistic.}

\item{hat_beta}{The original or the simulated estimated effects.}

\item{hat_sigma_sq}{The squared standard errors of the estimated effects.}
}
\value{
A list containing the p-value, bar_beta estimates.
}
\description{
This function performs a fixed effect analysis.
It calculates the p-value for the fixed effect model by comparing the test statistic
of the observed data with that of the simulated data.
}
\examples{
fixed_effect(10000, 2, c(0.4, 0.6), c(0.1, 0.2))

}
