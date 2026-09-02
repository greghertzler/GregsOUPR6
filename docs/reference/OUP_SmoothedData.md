# Smoothed sample paths for the Ornstein-Uhlenbeck Process

Monte Carlo simulation to demonstrate the effect of smoothing on
parameter estimates.

## Format

csv file with 177 rows and 11 columns

## Details

- year: time variable in annual increments

- Data: simulated sample path

- G uno-G nueve: nine successively smoother paths

Simulated data is smoothed ten times. First, the parameters are
estimated and the means calculated for each observation. Then the
calculated means are used, as if they are data, in a subsequent
estimation. Means are calculated again and used in the next estimation
and so on ten times. The true rate of convergence, rho, and location,
mu, are recovered from each estimation, but the scale, sigma, goes
toward zero.

By the eighth smoothing, the log of the likelihood becomes positive.
Hence, the likelihood, as the anti-log, becomes greater than one. In
other words, if the eighth and ninth smoothings were real samples, the
probability of observing them would be greater than 100%. Further
smoothings would increase this probability. A small sigma is a tell-tale
sign the data has been smoothed. A positive log likelihood is a sure
sign.

This is the best possible smoothing method. The model used for the
smoothing is consistent with the data. Surprisingly, hypothesis tests
and decision thresholds are not greatly affected. However, passage times
are calculated to be much larger and are completely unreliable.

The best possible smoothing is unlikely. Any model used to smooth the
data is probably not the Ornstein-Uhlenbeck Process. The estimates will
be wrong and the actual system will be much more uncertain. How much
more uncertain is uncertain.

Always use the raw data.
