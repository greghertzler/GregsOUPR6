# Brownian Motion

Brownian Motion is constructed from stationary increments.

## Format

csv file with 1001 rows and 4 columns

## Details

- tau: time variable

- e(s): stationary increments e(s)=eps(t-s)^0.5 where eps is a
  pseudo-random standard normal variable

- z(t)=z(s)+e(s): Brownian Motion from stationary increments

- z(t)-z(s)=e(s): residuals of Brownian Motion equal to the increments

Brownian Motion, also called a Wiener Process, does not converge and can
never be stationary. Paradoxically, it is constructed from stationary
increments. Beginning from a fixed state, say z(0), a stationary
increment is added to give the next state, z(1). That state becomes the
initial condition of the next observation to which another stationary
increment is added to give the next state, z(2). This can be solved to
show that z(2)=z(0)+e(0)+e(1). Continuing on, Brownian Motion at any
time is the initial fixed state plus a running total of the stationary
increments.

The test for Brownian Motion is the opposite of the test for stationary
increments. Parameter rho is expected to be zero. As rho goes to zero,
the mean of the Ornstein-Uhlenbeck Process goes to z(s). In other words,
the expected value of the next state equals the observed state at the
beginning of each observation. The variance goes to sigma^2(t-s), which
equals sigma^2 in this date with time intervals equal to one.

To conduct a test, first estimate the unrestricted parameters. Then
impose a restriction. Set rho=0. There is no need to restrict mu because
it cancels from all formulas and the estimation algorithm will set it to
zero for you.

Then do a Likelihood Ratio Test comparing the unrestricted and
restricted estimates. If the Log Likelihoods are not significantly
different, the hypothesis of Brownian Motion cannot be rejected.

Because the stationary increments are eps(t-s)^0.5, where eps~N(0,1) and
t-s=1, parameter sigma should not be significantly different from one.
