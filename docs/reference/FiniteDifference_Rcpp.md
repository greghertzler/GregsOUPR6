# FiniteDifference_Rcpp functions for numerical solutions of an Ornstein-Uhlenbeck Process

Calculations for the R6 class 'FiniteDifference', with sequential
processing.

## Usage

``` r
RcppOUPFDDrift(x,rho,mu)

RcppOUPFDDiffusion(x,sigma)

RcppOUPFDTerminalValue_Linear(x,xo,vs)

RcppOUPFDTerminalValue_Degenerate(x,xo,Vmax,Vmin)

RcppOUPFDTerminalValue_Stepped(x,xo,vs,Vmax,Vmin)

RcppOUPFDTerminalValue_Kinked(x,xo,vs,Vmax,Vmin)

RcppOUPFDTerminalValue_Butterfly(x,xo,xm,vs,Vmax,Vmin)

RcppOUPFDTerminalValue_Mitscherlich(x,xo,vr,Vmax,Vmin)

RcppOUPFDTerminalValue_Gompertz(x,xi,vr,Vmax,Vmin)

RcppOUPFDTerminalValue_Logistic(x,xi,vr,Vmax,Vmin)

RcppOUPFDTerminalValue_Transcendental(x,xo,xi,xm,Vmax,Vmin)

RcppOUPFDTerminalValue_YieldIndex(x,xo,xi,xm,Vmax,Vmin)

RcppOUPFDOption(s,x,V,r,theta,skip,rho,mu,sigma)

RcppOUPFDOptionEnvelope(s,x,V,r,theta,skip,rho,mu,sigma)

RcppOUPFDDecisionThreshold(x,V,OOenv,phi)
```

## Arguments

- x:

  vector of states

- rho:

  rate parameter 0\<=rho\<inf

- mu:

  location parameter -inf\<mu\<inf

- sigma:

  scale parameter -inf\<sigma\<inf

- xo:

  state at the intercept, step, kink or left wing

- vs:

  slope

- Vmax:

  maximum terminal value

- Vmin:

  minimum terminal value

- xm:

  state at the maximum or right wing

- vr:

  rate of change

- xi:

  state at the inflection point

- s:

  vector of times

- V:

  vector of terminal values

- r:

  discount rate 0\<=r\<inf

- theta:

  weight of current time in time stepping 0.5\<=theta\<=1

- skip:

  subdivide time intervals but report at times s 1\<=skip\<=1000

- OOenv:

  envelope of option values

- phi:

  search direction for exit or entry options

## Value

g(n) \<- RcppOUPFDDrift()

h2(n) \<- RcppOUPFDDiffusion()

V(n) \<- RcppOUPFDTerminalValue_Linear()

V(n) \<- RcppOUPFDTerminalValue_Degenerate()

V(n) \<- RcppOUPFDTerminalValue_Stepped()

V(n) \<- RcppOUPFDTerminalValue_Kinked()

V(n) \<- RcppOUPFDTerminalValue_Butterfly()

V(n) \<- RcppOUPFDTerminalValue_Mitscherlich()

V(n) \<- RcppOUPFDTerminalValue_Gompertz()

V(n) \<- RcppOUPFDTerminalValue_Logistic()

V(n) \<- RcppOUPFDTerminalValue_Transcendental()

V(n) \<- RcppOUPFDTerminalValue_YieldIndex()

c(m,n) \<- RcppOUPFDOption()

env(2,n) \<- RcppOUPFDOptionEnvelope()

dec(2) \<- RcppOUPFDDecisionThreshold()

## Notes on Values

Return values are vectors and matrices allocated in Rcpp. The dimensions
are shown for information. Of course, do not include them in R calls.
For example:

    g <- RcppOUPFDDrift(x,rho,mu)

The return values:

    g(n)
    h2(n)
    V(n)

are vectors of the same dimension as x.

The return value:

    c(m,n)

is a matrix with row dimension from s and column dimension from x.

The return value:

    env(2,n)

is a matrix with two row vectors for the option prices along the
envelope and the corresponding times. It is subset in R as:

    OOenv <- env[1,]
    tsenv <- env[2,]

The return value:

    dec(2)

is a vector with two elements for the state at the decision threshold
and the option price at that state. It is subset in R as:

    k <- dec[1]
    OOhat <- dec[2]

## Discussion

First, the finite difference method was implemented in R6 as a
single-threaded application. Then the R6 code was translated into Rcpp
sequential code. The finite difference method steps backward in time.
For each time, it uses LU decomposition, followed by forward and
backward substitution to solve matrices. The finite difference method is
inherently sequential. Global searches for envelopes of option prices
and decision thresholds are also sequential. There is no parallel
processing in this module. The R6 single-threaded code has been
archived, leaving the Rcpp code.

Below are microbenchmark median times to calculate 400,000 Option
Prices. Parameter skip=10 and only 40,000 Option Prices are reported.
Other rows show median times for the Option Envelope and the Decision
Threshold. Calculations are on an i7 CPU with a maximum speed of 4.5
GHz.

    Unit: milliseconds             R6     R6+
             function   single-thread    Rcpp
    -----------------------------------------
               Option        342.1143  7.7743
       OptionEnvelope        608.9655  4.6913
    DecisionThreshold          6.0607  0.1089

The Option Envelope is an input into the calculation of the Decision
Threshold. If the Decision Threshold is calculated from a standing
start, add the two times. R6+Rcpp calculates from 44.0 to 129.8 times
faster than R6 single-thread. The algorithms are not identical. R6+Rcpp
takes advantage of void subroutines and call-by-reference to manipulate
vectors and matrices without copying.

## From the Console

These functions are available in R, the RStudio console and RShiny apps.
For example, a calculation of option prices would be:

    s <- seq(from=30,to=10,by=-0.1)
    x <- seq(from=-40,to=60,by=0.5)
    n <- length(x)
    V <- vector("double",n)
    for(j in 1:n) { V[j] <- max(0,x[j]) }
    OO <- RcppOUPFDOption(s,x,V,0.05,0.5,10,0.5,15,15)

Calling functions from the console is slightly faster than calling them
through R6 objects. The table below compares the R6+Rcpp median times
above with the median times from the console:

    Unit: milliseconds     R6+  Console
             function     Rcpp     Rcpp
    -----------------------------------
               Option   7.7743   7.3196
       OptionEnvelope   4.6913   4.5448
    DecisionThreshold   0.1089   0.0174

There is less than half a millisecond time penalty for using the R6
object and it is much more convenient. It is reactive. In other words,
it stores inputs and outputs and maps inputs to outputs. Changing an
input will nullify dependent outputs, eliminating any danger of
reporting a stale output. Outputs are calculated only as needed and only
once. Then they are reused. In plots using Plotly, for example.

Potentially, the functions could be imported into other packages.
