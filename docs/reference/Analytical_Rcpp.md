# Analytical_Rcpp functions for calculating an Ornstein-Uhlenbeck Process

Calculations for the R6 class 'Analytical', with parallel processing.

## Usage

``` r
RcppOUPADrift(z,rho,mu)

RcppOUPADiffusion(z,sigma)

RcppOUPAMean(t,s,x,rho,mu,eps)

RcppOUPAVariance(t,s,rho,sigma,eps)

RcppOUPADensity(t,y,s,x,rho,mu,sigma)

RcppOUPAProbability(t,y,s,x,rho,mu,sigma,psi)

RcppOUPADoubleIntegral(t,y,s,x,rho,mu,sigma,psi)

RcppOUPAOption(s,x,t,y,rho,mu,sigma,r,phi,b,c)

RcppOUPAOptionEnvelope(s,x,t,y,rho,mu,sigma,r,phi,b,c)

RcppOUPAdOOdsZero(s,x,t,y,rho,mu,sigma,r,phi,b,c)

RcppOUPADecisionThreshold(y,rho,mu,sigma,r,phi,b,c)

RcppOUPAObligation(s,x,t,y,rho,mu,r,phi,b,c)

RcppOUPAPassageTimeModeMedianMean(k,s,x,omega,rho,mu,sigma,z)

RcppOUPAPassageTimePercentiles(k,s,x,omega,Ppct,rho,mu,sigma,z)

RcppOUPAPassageTimeDensity(t,k,s,x,omega,rho,mu,sigma,z)

RcppOUPAPassageTimeProbability(t,k,s,x,omega,rho,mu,sigma,z)
```

## Arguments

- z:

  vector of states or optional vector of alternate initial states

- rho:

  rate parameter 0\<=rho\<inf

- mu:

  location parameter -inf\<mu\<inf

- sigma:

  scale parameter -inf\<sigma\<inf

- t:

  terminal time or vector of forward times

- s:

  initial time or vector of backward times

- x:

  initial state or vector of backward states

- eps:

  proportion remaining after convergence 0\<=eps\<=1

- y:

  terminal state or vector of forward states

- psi:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- r:

  discount rate -inf\<r\<inf

- phi:

  \<=0 for integral -inf to x, \>0 for integral x to inf

- b:

  lump-sum benefit for entry option

- c:

  lump-sum cost for exit option

- k:

  decision threshold

- omega:

  degree of irreversibility 0\<=omega\<=1

- Ppct:

  passage time probability for a percentile 0.01\<=Ppct\<=0.99

## Value

g(n) \<- RcppOUPADrift()

h2(n) \<- RcppOUPADiffusion()

Gt(m+1) \<- RcppOUPAMean()

H2t(m+1) \<- RcppOUPAVariance()

p(m,n) \<- RcppOUPADensity()

P(m,n) \<- RcppOUPAProbability()

PP(m,n) \<- RcppOUPADoubleIntegral()

OO(m,n) \<- RcppOUPAOption()

OOs(2,n) \<- RcppOUPAOptionEnvelope()

dOOdszero(4,n+3) \<- RcppOUPAdOOdsZero()

kOO(2) \<- RcppOUPADecisionThreshold()

BC(m,n) \<- RcppOUPAObligation()

tmmmtmmmx(3,n+1) \<- RcppOUPAPassageTimeModeMedianMean()

tpcttpctx(3,n+1) \<- RcppOUPAPassageTimePercentiles()

ptptx(m,n+1) \<- RcppOUPAPassageTimeDensity()

PtPtx(m,n+1) \<- RcppOUPAPassageTimeProbability()

## Note on Notation

The notation for times and states is confusing. Times can be either
fixed or variable. States can be either observed or stochastic. As
arguments to the functions, they can be either scalars or vectors.
Functions in Rcpp have local scope and don't get confused. But this
documentation is constructed in the R manner, assuming all arguments are
globally defined.

Simple arguments have globally unique names, but times and states have
different roles, depending upon the problem. Initial value problems,
like probabilities, have an initial time, an initial state, a time
variable and a stochastic state. Terminal value problems, like options,
have a terminal time, a terminal state, a time variable and a stochastic
state. But initial value and terminal value problems are not different
problems. They are governed by the same stochastic differential and
integral equations. Times and states cannot swap places in the formulas
without violating the differential and integral equations.

Hence, we adopt Kolmogorov's notation with backward variables s and x
and forward variables t and y. In initial value problems, s is fixed, x
is observed, t is variable and y is stochastic. In terminal value
problems, t and y are fixed, s is variable and x is stochastic. But s,
x, t and y occupy their same places in all formulas.

The Rcpp functions aren't confused by this, but how can people know how
to enter the arguments? The time and state arguments listed first are
vectors. The time and state arguments listed second are scalars. For
example in RcppOUPAProbability(t,y,s,x,...), t and y are vectors and s
and x are scalars. In RcppOUPAOption(s,x,t,y,...), s and x are vectors
and t and y are scalars.

Finally, the state z is also schizophrenic. In the stochastic
differential equation, it represents either state x or state y. In
passage times, it is an optional argument for alternate initial states
x.

## Notes on Values

Return values are vectors and matrices allocated in Rcpp. The dimensions
are shown for information. Of course, do not include them in R calls.
For example:

    g <- RcppOUPADrift(z,rho,mu)

The return values:

    g(n)
    h2(n)

are vectors for drift and diffusion of the same dimension as z.

The return values:

    Gt(m+1)
    H2t(m+1)

are vectors with means G and Variances H2 followed by a time, teps. For
means, teps is the time until convergence to within epsilon of the
location. For Variances, teps is the time until convergence to within
epsilon of the asymptotic variance.

The return values:

    p(m,n)
    P(m,n)
    PP(m,n)

are matrices of Transition Densities, Transition Probabilities and
Double Integrals with rows for t and columns for y.

The return value:

    OO(m,n)

is a matrix of Options with rows for s and columns for x.

The return value:

    OOs(2,n)

is a matrix with two row vectors for option prices and corresponding
times along the option envelope. It is subset in R as:

    OOhat <- OOs[1,,drop=FALSE]
    shat <- OOs[2,,drop=FALSE]

where t is the terminal time.

The return value:

    dOOdszero(4,n+3)

is a matrix of four row vectors followed by three column vectors. The
row vectors are for option prices and times where the derivatives of
option prices with respect to times equal zero. The first and second
rows are for option prices. and times where option prices are convex in
time. The third and fourth rows are where option prices are concave in
time. The three column vectors are a patch to connect the row vectors
where the surface of the option prices transform from convex to concave.
Within the patch, the first row has option prices, the second row has
times and the third row has states. The matrix is subset in R as:

    dOOdsconvex <- dOOdszero[1:2,1:n,drop=FALSE]
    dOOdsconcave <- dOOdszero[3:4,1:n,drop=FALSE]
    dOOdspatch <- dOOdszero[1:3,(n+1):(n+3),drop=FALSE]

This return value is mostly a curiosity. Plotting it with the option
envelope shows that, at the decision threshold, the option envelope
jumps over the convexity to reach the terminal value.

The return value:

    kOO(2)

is a vector with two elements for the state at the decision threshold
and the option price at that state. It is subset in R as:

    k <- kOO[1]
    OOhat <- kOO[2]

The return values:

    tmmmtmmmx(3,n+1)
    tpcttpctx(3,n+1)

are matrices of three row vectors. The rows are for the mode, median and
mean or for the lower, median and upper percentiles. For each row, the
first n elements are for the alternate initial starting values, z. The
last element is for the initial state x. The matrices are subset in R
as:

    tmode <- tmmmtmmmx[1,n+1,drop=FALSE]
    tmedian <- tmmmtmmmx[2,n+1,drop=FALSE]
    tmean <- tmmmtmmmx[3,n+1,drop=FALSE]
    tmodes <- tmmmtmmmx[1,1:n,drop=FALSE]
    tmedians <- tmmmtmmmx[2,1:n,drop=FALSE]
    tmeans <- tmmmtmmmx[3,1:n,drop=FALSE]

    tlower <- tpcttpctx[1,n+1,drop=FALSE]
    tmedian <- tpcttpctx[2,n+1,drop=FALSE]
    tupper <- tpcttpctx[3,n+1,drop=FALSE]
    tlowers <- tpcttpctx[1,1:n,drop=FALSE]
    tmedians <- tpcttpctx[2,1:n,drop=FALSE]
    tuppers <- tpcttpctx[3,1:n,drop=FALSE]

The alternate initial states, z, are optional. If omitted, the return
values are vectors for the initial state, x, subset in R as:

    tmode <- tmmmtmmmx[1]
    tmedian <- tmmmtmmmx[2]
    tmean <- tmmmtmmmx[3]

    tlower <- tpcttpctx[1]
    tmedian <- tpcttpctx[2]
    tupper <- tpcttpctx[3]

The return values:

    ptptx(m,n+1)
    PtPtx(m,n+1)

are matrices plus column vectors. A matrix contains the passage time
densities or probabilities for alternate initial states, z. The last
column is the passage time density or probability at the initial state,
x. The matrices are subset in R as:

    pt <- ptptx[,1:n,drop=FALSE]
    ptx <- ptptx[,n+1,drop=FALSE]
    Pt <- PtPtx[,1:n,drop=FALSE]
    Ptx <- PtPtx[,n+1,drop=FALSE]

The alternate initial states, z, are optional. If omitted, the return
values are vectors for the initial state, x:

    ptptx(m)
    PtPtx(m)

These do not need to be subset in R.

## Discussion

First, the analytical formulas were implemented in R6 as a
single-threaded application. Then the R6 code was translated into Rcpp
sequential code. Then it was translated into RcppParallel code. The R6
single-threaded code was archived and only the Rcpp and RcppParallel
codes remain.

Below are microbenchmark median times to calculate 40,000 results.
Calculations are on an i7 CPU with 12 threads running at a maximum speed
of 4.5 GHz.

    Unit: milliseconds             R6      R6+           R6+
              function  single-thread     Rcpp  RcppParallel
    --------------------------------------------------------
               Density        225.656   7.2685        1.1125
           Probability       5412.144  11.6550        1.5818
        DoubleIntegral       6909.977  21.6572        2.6955
                Option       5187.448  24.1113        2.9611
            Obligation        347.682   3.9241        0.7776

R6+Rcpp calculates from 31.0 to 464.4 times faster than R6
single-thread. R6+RcppParallel calculates from 5.0 to 8.1 times faster
than R6+Rcpp and from 202.8 to 3421.5 times faster than R6
single-thread.

For 40,000 Options, the Option Envelope requires 200 searches in the
time direction and the Decision Threshold is a sequential search in the
state direction with nested searches in the time direction. There is no
parallel code for the Decision Threshold. Below are the median times:

    Unit: milliseconds             R6      R6+           R6+
              function  single-thread     Rcpp  RcppParallel
    --------------------------------------------------------
        OptionEnvelope        651.259   3.7036        0.7437
     DecisionThreshold        236.470   1.3693

R6+Rcpp calculates 174.5 and 172.7 times faster than R6 single-thread.
R6+RcppParallel calculates 5.0 times faster than R6+Rcpp and 875.7 times
faster than R6 single-thread.

Passage Time calculations are expensive. The mode, median and
percentiles are searches over densities and probabilities. The mean is a
Gaussian quadrature. Below are the median times for 10,000 Passage
Times:

    Unit: milliseconds                    R6      R6+           R6+
                     function  single-thread     Rcpp  RcppParallel
    ---------------------------------------------------------------
    PassageTimeModeMedianMean       4576.543  22.3153        2.9293
       PassageTimePercentiles       4246.912  17.4183        2.3781
           PassageTimeDensity       2091.740  11.1987        1.6769
       PassageTimeProbability       3217.157  11.2643        1.5363

R6+Rcpp calculates from 186.8 to 285.6 times faster than R6
single-thread. R6+RcppParallel calculates from 6.7 to 7.6 times faster
than R6+Rcpp and from 1247.4 to 2094.1 times faster than R6
single-thread.

After Rcpp versions of the functions were coded, all but
DecisionThreshold were translated into RcppParallel versions.
RccpParallel uses Intel's Threading Building Blocks (TBB) on the CPU.
Unlike parallel processing on a GPU or accelerator, memory isn't copied
and there is less overhead. On trivially small problems, sequential
versions calculate faster. On large problems, parallel versions
calculate several times faster.

RcppParallel is an optional package. If it is installed, it will be
used. Function RcppParallelInstalled() will enquire whether code is
compiled with RcppParallel or has fallen back to Rcpp.

## From the Console

These functions are available in R, the RStudio console and RShiny apps.
As an example, a calculation of option prices would be:

    s <- seq(from=30,to=10,by=-0.1)
    x <- seq(from=-40,to=60,by=0.5)
    OO <- RcppOUPAOption(s,x,30,0,0.5,15,15,0.05,1,0,0)

Calling functions directly from the console is slightly faster than
calling them indirectly through R6 objects. Here are microbenchmark
median times for comparison with some of the previous results:

    Unit: milliseconds           R6+        Console
              function  RcppParallel   RcppParallel
    -----------------------------------------------
               Density        1.1125         0.9925
           Probability        1.5818         1.4358
       Double Integral        2.6955         2.5753
                Option        2.9611         2.8299
            Obligation        0.7776         0.5291

The extra times taken by the R6 object are fractions of a millisecond.
Although slower, the R6 object has advantages over the console. All
inputs are optional and are coordinated across functions. Enter an input
once and calculate several outputs. The R6 object is reactive. In other
words, it stores the inputs and outputs and maps inputs to outputs. If
an input changes, dependent outputs are nullified and will be
recalculated, as requested, but nothing is calculated twice. The console
stores outputs in the global environment, but there is no map of inputs
to outputs. Outputs can be stale. Another advantage of the R6 object is
predefined plots with Plotly. The same simulation can plotted different
ways without recalculation.

Potentially, the Rcpp functions could be imported into other packages.
