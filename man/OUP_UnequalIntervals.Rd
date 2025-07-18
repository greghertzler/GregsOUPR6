% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MaximumLikelihoodHertz.R
\docType{data}
\name{OUP_UnequalIntervals}
\alias{OUP_UnequalIntervals}
\title{Unequal observation intervals and the Ornstein-Uhlenbeck Process}
\format{
csv file with 151 rows and 4 columns
}
\description{
Simulated data to demonstrate that observations are never 'missing'.
}
\details{
\itemize{
\item year 1: time variable in annual increments
\item z rho 0.5: sample paths with annual time increments or longer and a rate of convergence of 0.5
\item year 0.05: time variable in increments of 1/20 year
\item z rho 10: sample paths with time increments of 1/20 year or longer and a rate of convergence of 10
}

An observation consists of the initial state, the initial time, the terminal state
and the terminal time.  The terminal time minus the initial time can differ among
observations because each observation has its own mean and variance.  Longer
observation intervals have smaller means and larger variances but, otherwise,
parameters estimated from unequal observation intervals are statistically
equivalent to parameters estimated from equal intervals.

The usual time-series methods are special cases of estimating the Ornstein-Uhlenbeck
Process.  They estimate location, mu, and scale, sigma, but not the rate of
convergence, rho.  They eliminate rho by assuming weak stationarity.  The mean of
each observation is assumed to have converged and lost its connection to its
initial state.  The variance of each observation does not depend upon an initial
state and may still be converging.  Observations will have equal variances if
observation intervals are equal.

For the Ornstein-Uhlenbeck Process, variances converge twice as fast as means.
If variances are still converging, so are means.  Observations over equal time
intervals will have equal variances but different means.  In other words, weak
stationarity does not exist.  Any assumption of stationarity is strong stationarity
in which both the means and variances have converged.

Stationarity is an hypothesis to test, not an assumption to impose.  Goodness
of fit for the Ornstein-Uhlenbeck Process tests for stationarity, and also
for the other extreme of Brownian Motion.

If you have a time series with measurements at sporadic times, don't fill in for
observations that aren't missing.  If you must, first estimate the parameters.
Then use the parameters to calculate means for the times of your choosing.
Means are the maximum likelihood estimates of unobserved states of nature.
}
\keyword{datasets}
