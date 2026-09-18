# MonteCarlo_Rcpp functions for simulating an Ornstein-Uhlenbeck Process

Calculations for the R6 class 'MonteCarlo', with parallel processing.

## Usage

``` r
RcppOUPMCMinMax(matPaths)

RcppOUPMCStandardNormal(m,skip,paths,seed)

RcppOUPMCForwardPaths(stdnorm,x,m,skip,dt,rho,mu,sigma,method)

RcppOUPMCBackwardPaths(stdnorm,y,m,skip,ds,rho,mu,sigma,method)

RcppOUPMCBoundedPaths(stdnorm,k,x,m,skip,dt,rho,mu,sigma,method)

RcppOUPMCForwardCountY(forward,y,psi)

RcppOUPMCBackwardCountX(backward,x,phi,rho,r,ds)

RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)

RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)

RcppOUPMCHeatCountZ(matPaths,z)
```

## Arguments

- matPaths:

  matrix of paths

- m:

  number of rows for states over time

- skip:

  subdivide time interval but report every ds or dt 0\<skip\<20

- paths:

  number of columns for paths

- seed:

  seed for reproducibility

- engine:

  random number generator

- stdnorm:

  matrix of standard normal shocks

- x:

  initial state or vector of backward states

- dt:

  time interval for initial value problems

- rho:

  rate parameter 0\<=rho\<inf

- mu:

  location parameter -inf\<mu\<inf

- sigma:

  scale parameter -inf\<sigma\<inf

- method:

  4 for 4th order Runge-Kutta, 5 for integral equation

- y:

  terminal state or vector of forward states

- ds:

  time interval for terminal value problems

- k:

  threshold -inf\<k\<inf

- forward:

  matrix of forward paths

- psi:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- backward:

  matrix of backward paths

- phi:

  \<=0 for integral -inf to x, \>0 for integral x to inf

- r:

  discount rate 0\<r

- Ppct:

  probability for a percentile 0.01\<pct\<0.99

- fpt:

  vector of first passage times

- z:

  vector of states

## Value

minmax(2) \<- RcppOUPMCMinMax()

stdnorm((m-1)\*skip,paths) \<- RcppOUPMCStandardNormal()

forward(m,paths) \<- RcppOUPMCForwardPaths()

backward(m,paths) \<- RcppOUPMCBackwardPaths()

bndfpt(m+1,paths) \<- RcppOUPMCBoundedPaths()

mvdpd(m,3\*n+2) \<- RcppOUPMCForwardCountY()

dpo(m,3\*n) \<- RcppOUPMCBackwardCountX()

pctdp(m,5) \<- RcppOUPMCForwardCountT()

pctdp(m,5) \<- RcppOUPMCBoundedCountT()

heat(m,n) \<- RcppOUPMCHeatCountZ()

## Notes on Values

Return values are vectors and matrices allocated in Rcpp. The dimensions
are shown for information. Of course, do not include them in R calls.
For example:

    stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)

The return values:

    stdnorm((m-1)*skip,paths)
    forward(m,paths)
    backward(m,paths)
    heat(m,n)

are matrices of pseudo-random standard normal variables, forward paths,
backward paths and heat maps.

The return value:

    minmax(2)

is a vector of minimum and maximum values, subset in R as:

    minmax <- RcppOUPMCMinMax(matPaths)
    min <- minmax[1]
    max <- minmax[2]

The return value:

    bndfpt(m+1,paths)

is a matrix of bounded paths with a vector of first passage times for
each path in row m+1. There are NA entries for paths which previously
hit the threshold and NA entries in row m+1 for paths which have yet to
hit the threshold. Subset in R as:

    bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma)
    bounded <- bndfpt[1:m,,drop=FALSE]
    fpt <- bndfpt[m+1,,drop=FALSE]

The return value:

    mvdpd(m,3*n+2)

is a composite matrix containing two vectors for means and variances and
three matrices for densities, probabilities and double integrals. The
vectors and matrices are subset in R as:

    mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
    means <- mvdpd[,1,drop=FALSE]
    variances <- mvdpd[,2,drop=FALSE]
    densities <- mvdpd[,3:(n+2),drop=FALSE]
    probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
    doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]

Similarly, the return value:

    dpo(m,3*n)

is a composite of three contiguous matrices for prior densities, prior
probabilities and options, subset in R as:

    dpo <- RcppOUPMCBackwardCountX(backward,x,phi,rho,r,ds)
    densities <- dpo[,1:n,drop=FALSE]
    probabilities <- dpo[,(n+1):(2*n),drop=FALSE]
    options <- dpo[,(2*n+1):(3*n),drop=FALSE]

The return values:

    pctdp(m,5)

are matrices of five columns. The first column has only five entries for
times at the mode, median, mean, lower percentile, and upper percentile.
The second and third columns each have five entries containing the
corresponding densities and probabilities. The fourth and fifth columns
contain contain m entries for either visiting time or first passage time
densities and probabilities. The columns are subset in R as:

    pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
    mode <- pctdp[1,1:3,drop=FALSE]
    median <- pctdp[2,1:3,drop=FALSE]
    mean <- pctdp[3,1:3,drop=FALSE]
    lowerpct <- pctdp[4,1:3,drop=FALSE]
    upperpct <- pctdp[5,1:3,drop=FALSE]
    densities <- pctdp[,4,drop=FALSE]
    probabilities <- pctdp[,5,drop=FALSE]

## Discussion

A single-threaded R6 object is fast enough for many calculations, but
not for Monte-Carlo simulations. Attempts at parallel processing using
parApply() and future_apply() failed. The whole R6 object is copied to
each thread, which locks up the computer. Rccp can be hundreds of times
faster and makes Monte-Carlo simulations practical for interactive
applications such as RStudio and RShiny. RcppParallel speeds the
calculations another five to eight times.

For Monte Carlo simulations, the stochastic integral equation is shocked
by Brownian Motion. Brownian Motion is a time transform of standard
normal variables. The results are forward, backward and bounded paths.

Forward, backward and bounded paths are binned and counted to
approximate several solutions. The approximations converge to analytical
solutions as the number of paths increases. Binning and counting
1,000,000 paths will be accurate to 3 or 4 significant digits. Here are
microbenchmark median times for 100,000 and 1,000,000 paths over 100
time intervals, as calculated by R6+RccpParallel:

    Unit: milliseconds     paths                paths
              function   100,000            1,000,000
    ------------------------------------------------------------
        StandardNormal   21.6702             218.0919
          ForwardPaths   19.6883  ________   189.7758  _________
              Subtotal             41.4642              407.8677
           Probability   54.8138  ________   965.1143  _________
                 Total             96.2780             1372.9820

The R6 object is reactive and will call the StandardNormal function only
once. After that the standard normal variables will be passed to
ForwardPaths and the forward paths will be passed to Probability. To
calculate a median time, the R6 object is tricked into recalculating by
changing an input, calculating, changing the input back to the original,
and calculating the original again. Do this 11 times and record the
sixth fastest time as the median. Times for ForwardPaths, the Subtotal
and the Total are calculating by trickery. Times for StandardNormal and
Probability are inferred by subtraction.

Times for StandardNormal and ForwardPaths go up approximately ten-fold
with a ten-fold increase in paths. Times for Probability go up almost
18-fold with a ten-fold increase in paths.

The function Probability calls the Rcpp function ForwardCountY which
bins and counts means, variances, transition densities, transition
probabilities and double integrals. So five sets of plots can be drawn
from one simulation followed by a count. For comparison, a 3D plot by
Plotly can take up to a second on an RTX 2070 GPU. So calculations are
only part of the job.

The R6 object manages inputs and outputs and draws plots. All
calculations are in Rcpp and RcppParallel functions. RccpParallel uses
Intel's Threading Building Blocks (TBB) on the CPU. Unlike parallel
processing on a GPU or accelerator, memory isn't copied and there is
less overhead. On trivially small problems, sequential versions
calculate faster. On large problems, parallel versions calculate much
faster.

RcppParallel is an optional package. If it is installed, it will be
used. Function RcppParallelInstalled() will enquire whether code is
compiled with RcppParallel or has fallen back to Rcpp. Optional packages
for random number generation are dqrng and sitmo. The functions
RcppdqrngInstalled() and RcppsitmoInstalled() will enquire whether they
are installed.

## From the Console

Rcpp and RcppParallel functions are available in R, the RStudio console
and RShiny apps. From the console, a simulation of 1,000,000 forward
paths over 100 time intervals would be:

     stdnorm <- RcppOUPMCStandardNormal(101,1,1000000,9999,1)
     fwd <- RcppOUPMCForwardPaths(stdnorm,15,101,1,0.1,0.5,-15,15,5)

The R6 object doesn't give users a choice, but from the console there
are four random number generators available: dqrng, std::mt19937, sitmo,
and rnorm. In the last argument of the function, these are requested as
engine 1, 2, 3 or 4, respectively. Microbenchmark median times for
100,000,000 standard normal variables are:

    Unit: milliseconds
              language           rng    transform  StandardNormal
    -------------------------------------------------------------
                  Rcpp         rnorm    inversion       4145.0140
                  Rcpp   sitmo::prng   Box-Muller       3937.8950
                  Rcpp  std::mt19937   Box-Muller       3231.4860
                  Rcpp  dqrng::pcg64     Ziggurat        727.3865
          RcppParallel   sitmo::prng   Box-Muller        585.1609
          RcppParallel  std::mt19937   Box-Muller        547.9659
          RcppParallel  dqrng::pcg64     Ziggurat        247.1802

The random number generators generate uniform random variables. More
time is spent transforming uniform to normal random variables. The
method of transform is also listed. These include inverting the normal
probability, the Polar Box-Muller transform and the Ziggurat transform.
The Ziggurat transform is a sophisticated lookup table and much faster.
Results on your computer may vary. On Unix-alike operating systems, the
Polar Box-Muller transform has been replaced with the Ziggurat
transform.

Both dqrng and sitmo have other random number generators, but
dqrng::pcg64 and sitmo::prng are the defaults. Both packages are
optional. If dqrng is not installed, the fall back is std::mt19937,
which is always available. Therefore, sitmo::prng and rnorm will only
used if requested as engines 3 and 4. If sitmo::prng is requested but
not installed the fallback is rnorm.

Another choice available from the console is the method of simulation,
either a 4th-order Runge-Kutta numerical integration or the stochastic
integral equation, itself. The last argument of the function is the
method, with 4 for 4th-order Runge-Kutta and 5 for stochastic integral
equation. Arguments 1, 2 and 3 are reserved for possible future
implementations of 1st-order Euler and 2nd and 3rd order Maryuma
methods. But this could be dangerous. Users might use them. The purpose
would be to demonstrate that low-order numerical methods only converge
with short time intervals.

In the function, the skip argument divides the time intervals. For
example, if the number of times is 101, there are 100 time intervals.
Argument skip=10 subdivides 100 into 1000 time intervals for the
calculations and reports results at the 101 times.

Here are microbenchmark times for 1,000,000 paths over 100 time
intervals with increasing skips. Also shown are the maximum and minimum
differences.

    Unit: milliseconds   Standard   Integral     Runge-
                  skip     Normal   Equation      Kutta  max dif  min dif
    ---------------------------------------------------------------------
                     1   208.6749   190.2418   360.3870  8.4e-03  -8.9e-03
                     2   431.8739   214.1760   608.2229  2.1e-03  -2.1e-03
                     4   978.7199   295.2814  1120.7290  5.3e-04  -5.3e-04
                     8  2132.0220   424.4842  2086.3920  1.3e-04  -1.3e-04

Even larger skips will calculate, but microbenchmark becomes pac man and
starts chomping memory. For skip=8, the paths are the same to within
four significant digits. But the Runge-Kutta method is much slower. The
times for the standard normal variables and the Runge-Kutta simulation
takes 4.2 seconds. The integral equation is not improved by larger
skips. For skip=1, the integral equation does the job in 0.4 seconds.

A microbenchmark comparison of indirectly calling RcppParallel functions
from R6 with directly calling them from the console is:

    Unit: milliseconds            R6+       Console
              function   RcppParallel  RcppParallel
    -----------------------------------------------
          ForwardPaths       540.7106      587.6297
         BackwardPaths       543.9790      584.8225
          BoundedPaths       845.6546      579.5391

These timings are from a standing start, generating the standard normal
variables before simulating the paths. For Forward and Backward Paths,
the R6 object is faster, but for Bounded Paths, it is much slower.
Bounded Paths hit the boundary and have NA values thereafter. We might
speculate that R6 is slow with NA values.

The R6 object has advantages over the console. All inputs are optional
and are coordinated across functions. Enter an input once and calculate
several outputs. The R6 object is reactive. In other words, it stores
the inputs and outputs and maps inputs to outputs. If an input changes,
dependent outputs are nullified and will be recalculated, as requested,
but nothing is calculated twice. The console stores outputs in the
global environment, but there is no map of inputs to outputs and outputs
can be stale. Another advantage of the R6 object are pre-programmed
plots with Plotly. The same simulation can be plotted different ways
without recalculation.

Parallel processing has more overhead and is slower on small problems.
Here are microbenchmark median times for simulating a small number of
Forward Paths by calling the Rcpp and RcppParallel functions from the
console:

    Unit: milliseconds  Console       Console
                 paths     Rcpp   RcppParallel
    ------------------------------------------
                   100  0.06970        0.07600
                 1,000  0.55425        0.23540

For 100 paths, the Rcpp sequential function takes less time, but for
1,000 paths it takes over twice as long. Most simulations will have more
than 1,000 paths. So users get no choice. RcppParallel functions are
compiled if RcppParallel is installed. Otherwise compilation falls back
to Rcpp.

More threads may be faster but also have more overhead. Here are
microbenchmark median times by number of threads for generating standard
normal variables and simulating 1,000,000 Forward Paths over 100 time
intervals:

    Unit: milliseconds
               threads   stdnorm  ForwardPaths      total
    -----------------------------------------------------
                     1  742.8735        497.8163    1240.6898
                     2  422.2265        307.3954     729.6219
                     3  396.3184        307.5473     703.8657
                     4  257.9281        214.4442     472.3723
                     5  260.6473        197.7865     458.4338
                     6  248.4076        197.5342     445.9418
                     7  246.5812        191.5214     438.1026
                     8  231.4124        191.8745     423.2869
                     9  232.0408        194.9968     427.0376
                    10  235.8238        191.1864     427.0102
                    11  216.9428        192.1811     409.1239
                    12  211.6724        191.6413     403.3137

More is better, but a few is pretty good. You could set fewer threads
using the RcppParallel commands:

     library(RcppParallel)
     defaultNumThreads()
     setThreadOptions(numThreads=8)

Potentially, the Rcpp functions could be imported into other packages.
