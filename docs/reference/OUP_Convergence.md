# Rates of convergence for the Ornstein-Uhlenbeck Process

Monte-Carlo simulation to demonstrate different rates of convergence,
rho.

## Format

csv file with 177 rows and 16 columns

## Details

- year: time variable in annual increments for all sample paths

- z1-z5: sample paths in sets of three, each set with the same
  pseudo-random shocks but different rates of convergence

The rate of convergence, rho, determines the probability distribution of
the estimated parameters and the correlation between two sets of
parameters in hypothesis tests. Small rho tends toward Brownian Motion,
which does not converge. Large rho tends toward a stationary or ergodic
process which has converged every time it is observed. In between is an
Ornstein-Uhlenbeck Process which converges but has not yet converged.

Parameters for Browian Motion have an Erlang distribution. Parameters
for a stationary or ergodic process have a Chi^2 distribution. These
distributions are special cases of a Gamma distribution. In general,
parameters for the Ornstein-Uhlenbeck Process have a Gamma distribution.

The shape parameter, alpha, identifies the distribution with 0.5 \<=
alpha \<=1. Chi^2 has alpha = 0.5. Erlang has alpha = 1.
