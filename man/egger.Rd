% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_posterior.R
\name{egger}
\alias{egger}
\title{Modified Egger test statistics}
\usage{
egger(hat_beta, bar_beta, k, hat_sigma_sq, m)
}
\arguments{
\item{hat_beta}{The original or the simulated estimated effects.}

\item{bar_beta}{The estimated true underlying effect.}

\item{k}{The level of heterogeneity.}

\item{hat_sigma_sq}{The squared standard errors of the estimated effects.}

\item{m}{The number of replications.}
}
\value{
The Egger test statistic value.
}
\description{
This function provides the calculation for Egger test quantities.
}
\examples{
egger(c(0.4, 0.6), 0.5, 1, c(0.1, 0.2), 2)
}
