#include <Rcpp.h>
using namespace Rcpp;
#ifdef USE_PARALLEL
#include <RcppParallel.h>
using namespace RcppParallel;
#endif
#include <cmath>
#include <limits>
#include "gammahertz.h"

// roxygen (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @title Analytical_Rcpp functions for calculating an Ornstein-Uhlenbeck Process
//'
//' @description
//' Calculations for the R6 class 'Analytical', with parallel processing.
//'
//' @details # Note on Notation
//' The notation for times and states is confusing.  Times can be either
//'  fixed or variable.  States can be either observed or stochastic.  As
//'  arguments to the functions, they can be either scalars or vectors.
//'  Functions in Rcpp have local scope and don't get confused.  But this
//'  documentation is constructed in the R manner, assuming all arguments are
//'  globally defined.
//'
//' Simple arguments have globally unique names, but times and states have
//'  different roles, depending upon the problem.  Initial value problems, like
//'  probabilities, have an initial time, an initial state, a time variable and
//'  a stochastic state.  Terminal value problems, like options, have a terminal
//'  time, a terminal state, a time variable and a stochastic state.  But initial
//'  value and terminal value problems are not different problems.  They are
//'  governed by the same stochastic differential and integral equations. Times
//'  and states cannot swap places in the formulas without violating the differential
//'  and integral equations.
//'
//' Hence, we adopt Kolmogorov's notation with backward variables s and x and
//'  forward variables t and y.  In initial value problems, s is fixed, x is
//'  observed, t is variable and y is stochastic.  In terminal value problems,
//'  t and y are fixed, s is variable and x is stochastic.  But s, x, t and y
//'  occupy their same places in all formulas.
//'
//' The Rcpp functions aren't confused by this, but how can people know how to
//'  enter the arguments?  The time and state arguments listed first are vectors.
//'  The time and state arguments listed second are scalars.  For example in
//'  RcppOUPAProbability(t,y,s,x,...), t and y are vectors and s and x are scalars.
//'  In RcppOUPAOption(s,x,t,y,...), s and x are vectors and t and y are scalars.
//'
//' Finally, the state z is also schizophrenic.  In the stochastic differential
//'  equation, it represents either state x or state y.  In passage times, it
//'  is an optional argument for alternate initial states x.
//'
//' @details # Notes on Values
//' Return values are vectors and matrices allocated in Rcpp.  The dimensions are
//'  shown for information.  Of course, do not include them in R calls.  For example:
//'
//'     g <- RcppOUPADrift(z,rho,mu)
//'
//' The return values:
//'
//'     g(n)
//'     h2(n)
//'
//'  are vectors for drift and diffusion of the same dimension as z.
//'
//' The return values:
//'
//'     Gt(m+1)
//'     H2t(m+1)
//'
//'  are vectors with means G and Variances H2 followed by a time, teps.  For means,
//'  teps is the time until convergence to within epsilon of the location.  For Variances,
//'  teps is the time until convergence to within epsilon of the asymptotic variance.
//'
//' The return values:
//'
//'     p(m,n)
//'     P(m,n)
//'     PP(m,n)
//'
//'  are matrices of Transition Densities, Transition Probabilities and Double
//'  Integrals  with rows for t and columns for y.
//'
//' The return value:
//'
//'     OO(m,n)
//'
//'  is a matrix of Options with rows for s and columns for x.
//'
//' The return value:
//'
//'     OOs(2,n)
//'
//'  is a matrix with two row vectors for option prices and corresponding times
//'  along the option envelope.  It is subset in R as:
//'
//'     OOhat <- OOs[1,,drop=FALSE]
//'     shat <- OOs[2,,drop=FALSE]
//'
//'  where t is the terminal time.
//'
//' The return value:
//'
//'     dOOdszero(4,n+3)
//'
//'  is a matrix of four row vectors followed by three column vectors.  The row vectors
//'  are for option prices and times where the derivatives of option prices with
//'  respect to times equal zero.  The first and second rows are for option prices.
//'  and times where option prices are convex in time. The third and fourth rows are
//'  where option prices are concave in time.  The three column vectors are a patch to
//'  connect the row vectors where the surface of the option prices transform from
//'  convex to concave.  Within the patch, the first row has option prices, the second
//'  row has times and the third row has states.  The matrix is subset in R as:
//'
//'     dOOdsconvex <- dOOdszero[1:2,1:n,drop=FALSE]
//'     dOOdsconcave <- dOOdszero[3:4,1:n,drop=FALSE]
//'     dOOdspatch <- dOOdszero[1:3,(n+1):(n+3),drop=FALSE]
//'
//' This return value is mostly a curiosity. Plotting it with the option envelope shows
//'  that, at the decision threshold, the option envelope jumps over the convexity to
//'  reach the terminal value.
//'
//' The return value:
//'
//'     kOO(2)
//'
//'  is a vector with two elements for the state at the decision threshold and
//'  the option price at that state.  It is subset in R as:
//'
//'     k <- kOO[1]
//'     OOhat <- kOO[2]
//'
//' The return values:
//'
//'     tmmmtmmmx(3,n+1)
//'     tpcttpctx(3,n+1)
//'
//'  are matrices of three row vectors.  The rows are for the mode, median and mean
//'  or for the lower, median and upper percentiles.  For each row, the first n
//'  elements are for the alternate initial starting values, z.  The last element
//'  is for the initial state x.  The matrices are subset in R as:
//'
//'     tmode <- tmmmtmmmx[1,n+1,drop=FALSE]
//'     tmedian <- tmmmtmmmx[2,n+1,drop=FALSE]
//'     tmean <- tmmmtmmmx[3,n+1,drop=FALSE]
//'     tmodes <- tmmmtmmmx[1,1:n,drop=FALSE]
//'     tmedians <- tmmmtmmmx[2,1:n,drop=FALSE]
//'     tmeans <- tmmmtmmmx[3,1:n,drop=FALSE]
//'
//'     tlower <- tpcttpctx[1,n+1,drop=FALSE]
//'     tmedian <- tpcttpctx[2,n+1,drop=FALSE]
//'     tupper <- tpcttpctx[3,n+1,drop=FALSE]
//'     tlowers <- tpcttpctx[1,1:n,drop=FALSE]
//'     tmedians <- tpcttpctx[2,1:n,drop=FALSE]
//'     tuppers <- tpcttpctx[3,1:n,drop=FALSE]
//'
//' The alternate initial states, z, are optional.  If omitted, the return values
//'  are vectors for the initial state, x, subset in R as:
//'
//'     tmode <- tmmmtmmmx[1]
//'     tmedian <- tmmmtmmmx[2]
//'     tmean <- tmmmtmmmx[3]
//'
//'     tlower <- tpcttpctx[1]
//'     tmedian <- tpcttpctx[2]
//'     tupper <- tpcttpctx[3]
//'
//' The return values:
//'
//'     ptptx(m,n+1)
//'     PtPtx(m,n+1)
//'
//'  are matrices plus column vectors.  A matrix contains the passage time
//'  densities or probabilities for alternate initial states, z.  The last
//'  column is the passage time density or probability at the initial state, x.
//'  The matrices are subset in R as:
//'
//'     pt <- ptptx[,1:n,drop=FALSE]
//'     ptx <- ptptx[,n+1,drop=FALSE]
//'     Pt <- PtPtx[,1:n,drop=FALSE]
//'     Ptx <- PtPtx[,n+1,drop=FALSE]
//'
//' The alternate initial states, z, are optional.  If omitted, the return values
//'  are vectors for the initial state, x:
//'
//'     ptptx(m)
//'     PtPtx(m)
//'
//' These do not need to be subset in R.
//'
//' @details # Discussion
//' First, the analytical formulas were implemented in R6 as a single-threaded
//'  application.  Then the R6 code was translated into Rcpp sequential code.
//'  Then it was translated into RcppParallel code.  The R6 single-threaded
//'  code was archived and only the Rcpp and RcppParallel codes remain.
//'
//' Below are microbenchmark median times to calculate 40,000 results.
//'  Calculations are on an i7 CPU with 12 threads running at a maximum speed
//'  of 4.5 GHz.
//'
//'     Unit: milliseconds             R6      R6+           R6+
//'               function  single-thread     Rcpp  RcppParallel
//'     --------------------------------------------------------
//'                Density        225.656   7.2685        1.1125
//'            Probability       5412.144  11.6550        1.5818
//'         DoubleIntegral       6909.977  21.6572        2.6955
//'                 Option       5187.448  24.1113        2.9611
//'             Obligation        347.682   3.9241        0.7776
//'
//' R6+Rcpp calculates from 31.0 to 464.4 times faster than R6 single-thread.
//'  R6+RcppParallel calculates from 5.0 to 8.1 times faster than R6+Rcpp and
//'  from 202.8 to 3421.5 times faster than R6 single-thread.
//'
//' For 40,000 Options, the Option Envelope requires 200 searches in the time
//'  direction and the Decision Threshold is a sequential search in the state
//'  direction with nested searches in the time direction.  There is no parallel
//'  code for the Decision Threshold.  Below are the median times:
//'
//'     Unit: milliseconds             R6      R6+           R6+
//'               function  single-thread     Rcpp  RcppParallel
//'     --------------------------------------------------------
//'         OptionEnvelope        651.259   3.7036        0.7437
//'      DecisionThreshold        236.470   1.3693
//'
//' R6+Rcpp calculates 174.5 and 172.7 times faster than R6 single-thread.
//'  R6+RcppParallel calculates 5.0 times faster than R6+Rcpp and 875.7 times
//'  faster than R6 single-thread.
//'
//' Passage Time calculations are expensive.  The mode, median and percentiles
//'  are searches over densities and probabilities.  The mean is a Gaussian
//'  quadrature.  Below are the median times for 10,000 Passage Times:
//'
//'     Unit: milliseconds                    R6      R6+           R6+
//'                      function  single-thread     Rcpp  RcppParallel
//'     ---------------------------------------------------------------
//'     PassageTimeModeMedianMean       4576.543  22.3153        2.9293
//'        PassageTimePercentiles       4246.912  17.4183        2.3781
//'            PassageTimeDensity       2091.740  11.1987        1.6769
//'        PassageTimeProbability       3217.157  11.2643        1.5363
//'
//' R6+Rcpp calculates from 186.8 to 285.6 times faster than R6 single-thread.
//'  R6+RcppParallel calculates from 6.7 to 7.6 times faster than R6+Rcpp and
//'  from 1247.4 to 2094.1 times faster than R6 single-thread.
//'
//' After Rcpp versions of the functions were coded, all but DecisionThreshold
//'  were translated into RcppParallel versions.  RccpParallel uses Intel's
//'  Threading Building Blocks (TBB) on the CPU.  Unlike parallel processing on
//'  a GPU or accelerator, memory isn't copied and there is less overhead.  On
//'  trivially small problems, sequential versions calculate faster.  On large
//'  problems, parallel versions calculate several times faster.
//'
//' RcppParallel is an optional package.  If it is installed, it will be used.
//'  Function RcppParallelInstalled() will enquire whether code is compiled with
//'  RcppParallel or has fallen back to Rcpp.
//'
//' @details # From the Console
//' These functions are available in R, the RStudio console and RShiny apps. As
//'  an example, a calculation of option prices would be:
//'
//'     s <- seq(from=30,to=10,by=-0.1)
//'     x <- seq(from=-40,to=60,by=0.5)
//'     OO <- RcppOUPAOption(s,x,30,0,0.5,15,15,0.05,1,0,0)
//'
//' Calling functions directly from the console is slightly faster than
//'  calling them indirectly through R6 objects.  Here are microbenchmark
//'  median times for comparison with some of the previous results:
//'
//'     Unit: milliseconds           R6+        Console
//'               function  RcppParallel   RcppParallel
//'     -----------------------------------------------
//'                Density        1.1125         0.9925
//'            Probability        1.5818         1.4358
//'        Double Integral        2.6955         2.5753
//'                 Option        2.9611         2.8299
//'             Obligation        0.7776         0.5291
//'
//' The extra times taken by the R6 object are fractions of a millisecond.
//'  Although slower, the R6 object has advantages over the console. All inputs
//'  are optional and are coordinated across functions.  Enter an input once
//'  and calculate several outputs.  The R6 object is reactive.  In other words,
//'  it stores the inputs and outputs and maps inputs to outputs.  If an input
//'  changes, dependent outputs are nullified and will be recalculated, as
//'  requested, but nothing is calculated twice.  The console stores outputs in
//'  the global environment, but there is no map of inputs to outputs.  Outputs
//'  can be stale.  Another advantage of the R6 object is predefined plots with
//'  Plotly. The same simulation can plotted different ways without recalculation.
//'
//' Potentially, the Rcpp functions could be imported into other packages.
//'
//' @name Analytical_Rcpp

// Helpers (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

std::vector<double> Laguerret()
{
  std::vector<double> t(41);
  t[0] = 0.034840064125239;
  t[1] = 0.183625087016639;
  t[2] = 0.451525058905516;
  t[3] = 0.838985610907935;
  t[4] = 1.34657594330787;
  t[5] = 1.97503995935016;
  t[6] = 2.72530671207309;
  t[7] = 3.59849894479794;
  t[8] = 4.59594284949542;
  t[9] = 5.71917969983466;
  t[10] = 6.96997971310658;
  t[11] = 8.3503584847076;
  t[12] = 9.86259639570966;
  t[13] = 11.5092614819208;
  t[14] = 13.2932363681025;
  t[15] = 15.2177500186414;
  t[16] = 17.2864152453923;
  t[17] = 19.5032731583782;
  t[18] = 21.872846064999;
  t[19] = 24.4002007458958;
  t[20] = 27.0910246000436;
  t[21] = 29.9517179152153;
  t[22] = 32.9895065670353;
  t[23] = 36.21258090697;
  t[24] = 39.6302686603258;
  t[25] = 43.2532526217178;
  t[26] = 47.0938482901412;
  t[27] = 51.1663631195495;
  t[28] = 55.4875691072365;
  t[29] = 60.0773363233172;
  t[30] = 64.9595008936776;
  t[31] = 70.1630847824147;
  t[32] = 75.7240620856209;
  t[33] = 81.6880101000193;
  t[34] = 88.1142662803787;
  t[35] = 95.0828121482822;
  t[36] = 102.706501661439;
  t[37] = 111.154921958725;
  t[38] = 120.707582561042;
  t[39] = 131.899754362277;
  t[40] = 146.110597447904;

  return t;
}

std::vector<double> Laguerrew()
{
  std::vector<double> w(41);
  w[0] = 8.63554570194484E-02;
  w[1] = 0.173332645592903;
  w[2] = 0.208568103938104;
  w[3] = 0.193350430207399;
  w[4] = 0.147724394333778;
  w[5] = 9.56297394958616E-02;
  w[6] = 5.31761202739599E-02;
  w[7] = 2.55882342423587E-02;
  w[8] = 1.06989176609361E-02;
  w[9] = 3.89523849800191E-03;
  w[10] = 1.2358867321078E-03;
  w[11] = 3.41684428543532E-04;
  w[12] = 8.22451802053959E-05;
  w[13] = 1.72110639995837E-05;
  w[14] = 3.12503224715398E-06;
  w[15] = 4.91094984483102E-07;
  w[16] = 6.65939575826955E-08;
  w[17] = 7.76479538935683E-09;
  w[18] = 7.75316169614965E-10;
  w[19] = 6.59856132000171E-11;
  w[20] = 4.76126275521546E-12;
  w[21] = 2.8950494277635E-13;
  w[22] = 1.47313265756315E-14;
  w[23] = 6.22367295054339E-16;
  w[24] = 2.1634252311101E-17;
  w[25] = 6.12366679384877E-19;
  w[26] = 1.39455738111393E-20;
  w[27] = 2.51965913863939E-22;
  w[28] = 3.5530244291388E-24;
  w[29] = 3.8348925676424E-26;
  w[30] = 3.09502899341052E-28;
  w[31] = 1.81545615594484E-30;
  w[32] = 7.47169937604288E-33;
  w[33] = 2.06348976691315E-35;
  w[34] = 3.60876174544819E-38;
  w[35] = 3.69668912238723E-41;
  w[36] = 1.98741198787653E-44;
  w[37] = 4.75752031803417E-48;
  w[38] = 3.87248620144704E-52;
  w[39] = 6.42131713584797E-57;
  w[40] = 5.91215183605371E-63;

  return w;
}

double OUPDensity(double t, double y, double s, double x, double rho, double mu, double sigma, double dy)
{
  double pij;
  double variance;
  if(rho < 0.0000000001) { variance = sigma*sigma*(t-s); }
  else { variance = sigma*sigma/(2*rho)*(1.0-exp(-2.0*rho*(t-s))); }
  double mean = mu+(x-mu)*exp(-rho*(t-s));
  double ymean = y-mean;
  if(variance < dy*dy)
  {
    if(ymean >= -0.5*dy && ymean < 0.5*dy)
    {
      if(dy > 0) { pij = 1/dy; }
      else { pij = std::numeric_limits<double>::infinity(); }
    }
    else { pij = 0.0; }
  }
  else
  {
    double u = pow(6.28318530717958*variance,-0.5);
    double v2 = 0.5*ymean*ymean/variance;
    pij = u*exp(-v2);
  }
  return pij;
}

double OUPProbability(double t, double y, double s, double x, double rho, double mu, double sigma, double psi)
{
  double Pij;
  double variance;
  if(rho < 0.0000000001) { variance = sigma*sigma*(t-s); }
  else { variance = sigma*sigma/(2*rho)*(1.0-exp(-2.0*rho*(t-s))); }
  double mean = mu+(x-mu)*exp(-rho*(t-s));
  double ymean = y-mean;
  if(variance < 0.0000000001)
  {
    if(ymean < -0.0000000001) { Pij = 0.0; }
    else if(ymean < 0.0000000001) { Pij = 0.5; }
    else { Pij = 1.0; }
    if(psi > 0) { Pij = 1.0-Pij; }

  }
  else if(psi <= 0)
  {
    double v2 = 0.5*ymean*ymean/variance;
    if(y < mean) { Pij = 0.5*GammaBigOneHalf(v2)/1.77245385090552; }
    else if(y == mean) { Pij = 0.5; }
    else { Pij = 0.5*(1+GammaSmallOneHalf(v2)/1.77245385090552); }
  }
  else
  {
    double v2 = 0.5*ymean*ymean/variance;
    if(y > mean) { Pij = 0.5*GammaBigOneHalf(v2)/1.77245385090552; }
    else if(y == mean) { Pij = 0.5; }
    else { Pij = 0.5*(1+GammaSmallOneHalf(v2)/1.77245385090552); }
  }
  return Pij;
}

double OUPDoubleIntegral(double t, double y, double s, double x, double rho, double mu, double sigma, double dy, double psi)
{
  double PPij;
  int ltgt = 1;
  if(psi > 0) { ltgt = -1; }
  if(t <= s) { PPij = ltgt*(y-x)*OUPProbability(t,y,s,x,rho,mu,sigma,psi); }
  else if(sigma == 0) { PPij = ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*OUPProbability(t,y,s,x,rho,mu,sigma,psi); }
  else if(rho < 0.0000000001) { PPij = ltgt*(y-x)*OUPProbability(t,y,s,x,rho,mu,sigma,psi)+sigma*sigma*(t-s)*OUPDensity(t,y,s,x,rho,mu,sigma,dy); }
  else { PPij = ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*OUPProbability(t,y,s,x,rho,mu,sigma,psi)+sigma*sigma/(2*rho)*(1-exp(-2*rho*(t-s)))*OUPDensity(t,y,s,x,rho,mu,sigma,dy); }

  return PPij;
}

double OUPOption(double s, double x, double t, double y, double rho, double mu, double sigma, double dy, double r, double phi, double b, double c)
{
  double OOij;
  int ltgt = 1;
  double bc = -c;
  if(phi > 0)
  {
    ltgt = -1;
    bc = b;
  }
  if(t <= s) { OOij = ltgt*(y-x)*OUPProbability(t,y,s,x,rho,mu,sigma,phi);  }
  else if(sigma == 0) { OOij = ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*OUPProbability(t,y,s,x,rho,mu,sigma,phi); }
  else if(rho < 0.0000000001) { OOij = ltgt*(y-x)*OUPProbability(t,y,s,x,rho,mu,sigma,phi)+sigma*sigma*(t-s)*OUPDensity(t,y,s,x,rho,mu,sigma,dy); }
  else { OOij = ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*OUPProbability(t,y,s,x,rho,mu,sigma,phi)+sigma*sigma/(2*rho)*(1-exp(-2*rho*(t-s)))*OUPDensity(t,y,s,x,rho,mu,sigma,dy); }

  return exp(-r*(t-s))*(OOij+bc);
}

double OUPObligation(double s, double x, double t, double y, double rho, double mu, double r, double phi, double b, double c)
{
  double BCij;
  if(phi <= 0) { BCij = (mu+(x-mu)*exp(-rho*(t-s))-y+b+c); }
  else { BCij = -(mu+(x-mu)*exp(-rho*(t-s))-y+b+c); }

  return exp(-r*(t-s))*BCij;
}

double OUPPassageTimeDensity(double s, double x, double t, double k, double omega, double rho, double mu, double sigma, double dt)
{
  double density = 0.0;
  if(k < std::numeric_limits<double>::infinity() && k > -std::numeric_limits<double>::infinity())
  {
    double varx;
    if(rho < 0.0000000001) { varx = sigma*sigma*(t-s); }
    else { varx = sigma*sigma*(1-exp(-2*rho*(t-s)))/(2*rho); }
    if(varx < 0.0000000001)
    {
      density = 0.0;
      if(x < mu)
      {
        if(k > mu+(x-mu)*exp(-rho*(t-0.5*dt-s)) && k <= mu+(x-mu)*exp(-rho*(t+0.5*dt-s))) { density = 10.0; }
      }
      else
      {
        if(k > mu+(x-mu)*exp(-rho*(t+0.5*dt-s)) && k <= mu+(x-mu)*exp(-rho*(t-0.5*dt-s))) { density = 10.0; }
      }
    }
    else if(x == k)
    {
      double u = abs(k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))*sigma*sigma/pow(varx,1.5);
      double v2 = (k-mu)*(1-exp(-rho*(t-s)));
      v2 = 0.5*v2*v2/varx;
      density = (1-omega)*u*exp(-v2)/(2*2.50662827431001);
    }
    else if(k == mu)
    {
      double u = abs(x-k)*exp(-rho*(t-s))*sigma*sigma/pow(varx,1.5);
      double v2 = (k-x)*exp(-rho*(t-s));
      v2 = 0.5*v2*v2/varx;
      density = (1+omega)*u*exp(-v2)/(2*2.50662827431001);
    }
    else
    {
      double meanx = mu+(x-mu)*exp(-rho*(t-s));
      double boundx = mu+(k-mu)*exp(-rho*(t-s))-(x-k)*exp(-rho*(t-s));
      double vm2 = 0.5*(k-meanx)*(k-meanx)/varx;
      double vb2 = 0.5*(k-boundx)*(k-boundx)/varx;
      double smallgammab = GammaSmallOneHalf(vb2);
      double biggammab = GammaBigOneHalf(vb2);
      double lnlambda = 2*(k-mu)*(1-exp(-rho*(t-s)))*(x-k)*exp(-rho*(t-s))/varx;
      double dlnlambda = -(k-mu)*(x-k)*sigma*sigma*(exp(-rho*(t-s))-2*exp(-2*rho*(t-s))+exp(-3*rho*(t-s)))/(varx*varx);
      if(x > k)
      {
        double u = ((x-mu)*exp(-rho*(t-s))-(k-mu)*exp(-2*rho*(t-s))+omega*((x-k)*exp(-rho*(t-s))-(k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))))*sigma*sigma/pow(varx,1.5);
        if(k < mu && k < boundx) { density = u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552); }
        else
        {
          if(lnlambda > 709.782712893384) { density = u*exp(-vm2)/(2*2.50662827431001); }
          else { density = u*exp(-vm2)/(2*2.50662827431001)+omega*biggammab*dlnlambda*exp(lnlambda)/(2*1.77245385090552); }
        }
      }
      else
      {
        double u = ((k-mu)*exp(-2*rho*(t-s))-(x-mu)*exp(-rho*(t-s))+omega*((k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))-(x-k)*exp(-rho*(t-s))))*sigma*sigma/pow(varx,1.5);
        if(k > mu && k > boundx) { density = u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552); }
        else
        {
          if(lnlambda > 709.782712893384) { density = u*exp(-vm2)/(2*2.50662827431001); }
          else { density = u*exp(-vm2)/(2*2.50662827431001)+omega*biggammab*dlnlambda*exp(lnlambda)/(2*1.77245385090552); }
        }
      }
    }
  }
  return density;
}

double OUPPassageTimeProbability(double s, double x, double t, double k, double omega, double rho, double mu, double sigma)
{
  double probability = 0.0;
  if(k < std::numeric_limits<double>::infinity() && k > -std::numeric_limits<double>::infinity())
  {
    double varx;
    if(x == k)
    {
      if(k == mu) { probability = 0.5*(1+omega); }
      else
      {
        if(rho < 0.0000000001) { varx = sigma*sigma*(t-s); }
        else { varx = sigma*sigma*(1-exp(-2*rho*(t-s)))/(2*rho); }
        if(varx < 0.0000000001) { probability = 0.5*(1+omega); }
        else
        {
          double v2 = (k-mu)*(1-exp(-rho*(t-s)));
          v2 = 0.5*v2*v2/varx;
          double smallgamma = GammaSmallOneHalf(v2);
          double biggamma = GammaBigOneHalf(v2);
          probability = (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552);
        }
      }
    }
    else if(k == mu)
    {
      if(rho < 0.0000000001) { varx = sigma*sigma*(t-s); }
      else { varx = sigma*sigma*(1-exp(-2*rho*(t-s)))/(2*rho); }
      if(varx < 0.0000000001) { probability = 0.0; }
      else
      {
        double v2 = (k-x)*exp(-rho*(t-s));
        v2 = 0.5*v2*v2/varx;
        double biggamma = GammaBigOneHalf(v2);
        probability = (1+omega)*biggamma/(2*1.77245385090552);
      }
    }
    else
    {
      double meanx = mu+(x-mu)*exp(-rho*(t-s));
      double boundx = mu+(k-mu)*exp(-rho*(t-s))-(x-k)*exp(-rho*(t-s));
      if(rho < 0.0000000001) { varx = sigma*sigma*(t-s); }
      else { varx = sigma*sigma*(1-exp(-2*rho*(t-s)))/(2*rho); }
      if(varx < 0.0000000001)
      {
        if(x > k && k > mu)
        {
          if(k < meanx) { probability = 0.0; }
          else if(k == meanx) { probability = 0.5; }
          else { probability = 1.0; }
        }
        else if(x > k && mu > k) { probability = 0.0; }
        else if(x < k && k < mu)
        {
          if(k > meanx) { probability = 0.0; }
          else if(k == meanx) { probability = 0.5; }
          else { probability = 1.0; }
        }
        else { probability = 0.0; }
      }
      else
      {
        double vm2 = k-meanx;
        vm2 = 0.5*vm2*vm2/varx;
        double vb2 = k-boundx;
        vb2 = 0.5*vb2*vb2/varx;
        double lnlambda = 2*(k-mu)*(1-exp(-rho*(t-s)))*(x-k)*exp(-rho*(t-s))/varx;
        double smallgammam = GammaSmallOneHalf(vm2);
        double biggammam = GammaBigOneHalf(vm2);
        double smallgammab = GammaSmallOneHalf(vb2);
        double biggammab = GammaBigOneHalf(vb2);
        if(x > k && k > mu)
        {
          if(k < meanx)
          {
            if(lnlambda > 709.782712893384) { probability = biggammam/(2*1.77245385090552); }
            else { probability = (biggammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552); }
          }
          else
          {
            if(lnlambda > 709.782712893384) { probability = (1.77245385090552+smallgammam)/(2*1.77245385090552); }
            else { probability = (1.77245385090552+smallgammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552); }
          }
        }
        else if(x > k && mu > k)
        {
          if(k > boundx) { probability = (biggammam+omega*exp(lnlambda)*biggammab)/(2*1.77245385090552); }
          else { probability = (biggammam+omega*exp(lnlambda)*(1.77245385090552+smallgammab))/(2*1.77245385090552); }
        }
        else if(x < k && k < mu)
        {
          if(k > meanx)
          {
            if(lnlambda > 709.782712893384) { probability = biggammam/(2*1.77245385090552); }
            else { probability = (biggammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552); }
          }
          else
          {
            if(lnlambda > 709.782712893384) { probability = (1.77245385090552+smallgammam)/(2*1.77245385090552); }
            else { probability = (1.77245385090552+smallgammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552); }
          }
        }
        else
        {
          if(k < boundx) { probability = (biggammam+omega*exp(lnlambda)*biggammab)/(2*1.77245385090552); }
          else { probability = (biggammam+omega*exp(lnlambda)*(1.77245385090552+smallgammab))/(2*1.77245385090552); }
        }
      }
    }
  }
  return probability;
}

double OUPPassageTimeProbabilityInf(double x, double k, double omega, double rho, double mu, double sigma)
{
  double PInf = 0.0;
  if(k < std::numeric_limits<double>::infinity() || k > -std::numeric_limits<double>::infinity())
  {
    if(omega == 1) { PInf = 1.0; }
    else if(k == mu) { PInf = 0.5*(1+omega); }
    else if(sigma*sigma < 0.0000000001)
    {
      if(x == k) { PInf = 1.0; }
      else if(x > k)
      {
        if(k > mu) { PInf = 1.0; }
        else { PInf = omega; }
      }
      else
      {
        if(k < mu) { PInf = 1.0; }
        else { PInf = omega; }
      }
    }
    else
    {
      double v2 = (k-mu)/sigma;
      v2 = rho*v2*v2;
      double smallgamma = GammaSmallOneHalf(v2);
      double biggamma = GammaBigOneHalf(v2);
      if(x == k) { PInf = (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552); }
      else if(x > k)
      {
        if(k > mu) { PInf = (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552); }
        else { PInf = (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552); }
      }
      else if(x < k)
      {
        if(k < mu) { PInf = (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552); }
        else { PInf = (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552); }
      }
    }
  }
  return PInf;
}

std::vector<double> OUPOptionMaxMin(double x, double y, double rho, double mu, double sigma, double dy, double r, double phi, double b, double c, double tsguess, double tsmax, double tsmin, double maxmin)
{
  std::vector<double> OOshat(2);
  if(rho < 0.0000000001 && r < 0.0000000001)
  {
    if(maxmin > 0 && sigma > 0)
    {
      OOshat[0] = std::numeric_limits<double>::infinity();
      OOshat[1] = -std::numeric_limits<double>::infinity();
    }
    else if(phi > 0)
    {
      if(x > y) { OOshat[0] = x-y+b; }
      else { OOshat[0] = b; }
      OOshat[1] = 0.0;
    }
    else
    {
      if(y > x) { OOshat[0] = y-x-c; }
      else { OOshat[0] = -c; }
      OOshat[1] = 0.0;
    }
  }
  else
  {
    double tensigdig = 10000000000000;
    std::size_t n = 1000;
    std::vector<double> Optn(5);
    std::vector<double> ts(5);
    if(tsguess > tsmin && tsguess < tsmax) { ts[2] = tsguess; }
    else { ts[2] = (tsmax+tsmin)/2.0; }
    Optn[2] = OUPOption(-ts[2],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
    double dt = (ts[2]-tsmin)/4.0;
    if(ts[2] > (tsmax+tsmin)/2.0) { dt = (tsmax-ts[2])/4.0; }
    ts[1] = ts[2]-dt;
    ts[0] = ts[1]-dt;
    ts[3] = ts[2]+dt;
    ts[4] = ts[3]+dt;
    Optn[0] = OUPOption(-ts[0],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
    Optn[1] = OUPOption(-ts[1],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
    Optn[3] = OUPOption(-ts[3],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
    Optn[4] = OUPOption(-ts[4],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
    // Rcout << ts[0] << "," << Optn[0] << ", " << ts[1] << "," << Optn[1] << ", " << ts[2] << "," << Optn[2] << ", " << ts[3] << "," << Optn[3] << ", " << ts[4] << "," << Optn[4] << std::endl;
    std::size_t m = 0;
    for(std::size_t i = 1; i < 5; i++)
    {
      if(maxmin*(Optn[i]-Optn[m]) > 0) { m = i; }
    }
    std::size_t j = 0;
    while(tensigdig*abs(Optn[3]-Optn[1]) >= abs(Optn[2]) && ts[2] < 10.0*tsmax && ts[2] > tsmin && j < n)
    {
      j += 1;
      if(m == 1 || m == 2 || m == 3)
      {
        if(m == 1)
        {
          ts[4] = ts[2];
          Optn[4] = Optn[2];
          ts[2] = ts[1];
          Optn[2] = Optn[1];
        }
        else if(m == 2)
        {
          ts[0] = ts[1];
          Optn[0] = Optn[1];
          ts[4] = ts[3];
          Optn[4] = Optn[3];
        }
        else
        {
          ts[0] = ts[2];
          Optn[0] = Optn[2];
          ts[2] = ts[3];
          Optn[2] = Optn[3];
        }
        ts[1] = (ts[2]+ts[0])/2;
        Optn[1] = OUPOption(-ts[1],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
        ts[3] = (ts[2]+ts[4])/2;
        Optn[3] = OUPOption(-ts[3],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
        dt *= 0.5;
        if(maxmin*(Optn[1]-Optn[2]) > 0) { m = 1; }
        else if(maxmin*(Optn[3]-Optn[2]) > 0) { m = 3; }
        else { m = 2; }
      }
      else if(m == 4)
      {
        for(std::size_t i = 1; i < 5; i++)
        {
          ts[i-1] = ts[i];
          Optn[i-1] = Optn[i];
        }
        ts[4] = ts[4]+dt;
        Optn[4] = OUPOption(-ts[4],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
        if(maxmin*(Optn[4]-Optn[3]) > 0)
        {
          m = 4;
          dt *= 2;
        }
        else { m = 3; }
      }
      else
      {
        for(std::size_t i = 4; i > 1; i--)
        {
          ts[i] = ts[i-1];
          Optn[i] = Optn[i-1];
        }
        if(ts[0] > tsmin)
        {
          ts[1] = ts[0];
          Optn[1] = Optn[0];
          ts[0] = ts[0]-dt;
          if(ts[0] < tsmin) { ts[0] = tsmin; }
          Optn[0] = OUPOption(-ts[0],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
        }
        else
        {
          ts[1] = (ts[1]+ts[0])/2;
          Optn[1] = OUPOption(-ts[1],x,0,y,rho,mu,sigma,dy,r,phi,b,c);
        }
        if(maxmin*(Optn[0]-Optn[1]) > 0)
        {
          m = 0;
          dt *= 2;
        }
        else { m = 1; }
      }
    }
    // Rcout << ts[0] << "," << Optn[0] << ", " << ts[1] << "," << Optn[1] << ", " << ts[2] << "," << Optn[2] << ", " << ts[3] << "," << Optn[3] << ", " << ts[4] << "," << Optn[4] << std::endl << std::endl;
    if(maxmin*(Optn[2]-Optn[0]) <= 0)
    {
      OOshat[0] = Optn[0];
      OOshat[1] = ts[0];
    }
    else
    {
      OOshat[0] = Optn[2];
      OOshat[1] = ts[2];
    }
  }
  return OOshat;
}

double OUPPassageTimeModeSearch(double s, double x, double k, double omega, double rho, double mu, double sigma, double median)
{
  double mode = median;
  if(median < std::numeric_limits<double>::infinity() && abs(x - k) >= 0.0000000001 && sigma*sigma >= 0.0000000001)
  {
    std::vector<int> m(2);
    std::vector<double> dt(2);
    std::vector<double> Dens(10);
    std::vector<double> t(10);
    double fivesigdig = 1000000;
    std::size_t n = 1000;
    dt[0] = 0.01*(median-s);
    dt[1] = 0.5*(median-s);
    t[0] = 0.0;
    t[0+5] = 0.0;
    Dens[0] = 0.0;
    Dens[0+5] = 0.0;
    m[0] = 0;
    m[1] = 0;
    for(std::size_t r = 1; r < 5; r++)
    {
      for(std::size_t c = 0; c < 2; c++)
      {
        std::size_t co = 5*c;
        t[r+co] = t[r-1+co]+dt[c];
        Dens[r+co] = OUPPassageTimeDensity(0,x,t[r+co],k,omega,rho,mu,sigma,0.05);
        if(Dens[r+co] > Dens[m[c]+co]) { m[c] = r; }
      }
    }
    std::size_t c = 0;
    if(Dens[m[0]] < 0.000001) { c = 1; }
    while(c < 2)
    {
      std::size_t co = 5*c;
      std::size_t j = 0;
      while(fivesigdig*(t[3+co]-t[1+co]) >= t[2+co] && j < n)
      {
        j = j+1;
        if(m[c] == 1 || m[c] == 2 || m[c] == 3)
        {
          if(m[c] == 1)
          {
            t[4+co] = t[2+co];
            Dens[4+co] = Dens[2+co];
            t[2+co] = t[1+co];
            Dens[2+co] = Dens[1+co];
          }
          else if(m[c] == 2)
          {
            t[0+co] = t[1+co];
            Dens[0+co] = Dens[1+co];
            t[4+co] = t[3+co];
            Dens[4+co] = Dens[3+co];
          }
          else
          {
            t[0+co] = t[2+co];
            Dens[0+co] = Dens[2+co];
            t[2+co] = t[3+co];
            Dens[2+co] = Dens[3+co];
          }
          t[1+co] = (t[2+co]+t[0+co])/2;
          Dens[1+co] = OUPPassageTimeDensity(0,x,t[1+co],k,omega,rho,mu,sigma,0.05);
          t[3+co] = (t[2+co]+t[4+co])/2;
          Dens[3+co] = OUPPassageTimeDensity(0,x,t[3+co],k,omega,rho,mu,sigma,0.05);
          if(Dens[1+co] > Dens[2+co]) { m[c] = 1; }
          else if(Dens[3+co] > Dens[2+co]) { m[c] = 3; }
          else { m[c] = 2; }
        }
        else if(m[c] == 4)
        {
          for(std::size_t r = 1; r < 5; r++)
          {
            t[r-1+co] = t[r+co];
            Dens[r-1+co] = Dens[r+co];
          }
          t[4+co] = t[4+co]+dt[c];
          Dens[4+co] = OUPPassageTimeDensity(0,x,t[4+co],k,omega,rho,mu,sigma,0.05);
          if(Dens[4+co] > Dens[3+co])
          {
            m[c] = 4;
            dt[c] = 2*dt[c];
          }
          else { m[c] = 3; }
        }
        else
        {
          for(std::size_t r = 4; r > 0; r--)
          {
            t[r+co] = t[r-1+co];
            Dens[r+co] = Dens[r-1+co];
          }
          t[0+co] = t[0+co]-dt[c];
          Dens[0+co] = OUPPassageTimeDensity(0,x,t[0+co],k,omega,rho,mu,sigma,0.05);
          if(Dens[0+co] > Dens[1+co])
          {
            m[c] = 0;
            dt[c] = 2*dt[c];
          }
          else { m[c] = 1; }
        }
      }
      c += 1;
    }
    if(Dens[2] > Dens[2+5]) { mode = s+t[2]; }
    else { mode = s+t[2+5]; }
  }
  return mode;
}

double OUPPassageTimePctSearch(double s, double x, double k, double omega, double rho, double mu, double sigma, double Ppct)
{
  double tpct = std::numeric_limits<double>::infinity();
  if(k < std::numeric_limits<double>::infinity() && k > -std::numeric_limits<double>::infinity())
  {
    if(abs(x - k) < 0.0000000001) { tpct = s; }
    else if(sigma*sigma < 0.0000000001)
    {
      if(((x < k && k < mu) || (x > k && k > mu)) && rho > 0) { tpct = s + log((x-mu)/(k-mu))/rho; }
      else if(omega == 0) { tpct = s; }
    }
    else
    {
      double fivesigdig = 1000000;
      std::size_t n = 1000;
      double PInf = OUPPassageTimeProbabilityInf(x,k,omega,rho,mu,sigma);
      std::vector<double> Prob(5);
      std::vector<double> t(5);
      t[0] = 0.0;
      Prob[0] = 0.0;
      double dt = 0.1;
      for(std::size_t i = 1; i < 5; i++)
      {
        t[i] = dt;
        Prob[i] = OUPPassageTimeProbability(0,x,t[i],k,omega,rho,mu,sigma);
        dt = dt*10.0;
      }
      std::size_t m = 0;
      while(Prob[m] < Ppct*PInf && m < 4) { m +=1; }
      if(m > 0)
      {
        std::size_t j = 0;
        while(fivesigdig*(t[3]-t[1]) >= t[2] && j < n)
        {
          j = j+1;
          if(m == 1 || m == 2 || m == 3)
          {
            if(m == 1)
            {
              t[4] = t[2];
              Prob[4] = Prob[2];
              t[2] = t[1];
              Prob[2] = Prob[1];
            }
            else if(m == 2)
            {
              t[0] = t[1];
              Prob[0] = Prob[1];
              t[4] = t[3];
              Prob[4] = Prob[3];
            }
            else
            {
              t[0] = t[2];
              Prob[0] = Prob[2];
              t[2] = t[3];
              Prob[2] = Prob[3];
            }
            t[1] = (t[2]+t[0])/2;
            Prob[1] = OUPPassageTimeProbability(0,x,t[1],k,omega,rho,mu,sigma);
            t[3] = (t[2]+t[4])/2;
            Prob[3] = OUPPassageTimeProbability(0,x,t[3],k,omega,rho,mu,sigma);
            if(Prob[1] > Ppct*PInf) { m = 1; }
            else if(Prob[3] < Ppct*PInf) { m = 3; }
            else { m = 2; }
          }
          else if(m == 4)
          {
            for(std::size_t i = 1; i < 5; i++)
            {
              t[i-1] = t[i];
              Prob[i-1] = Prob[i];
            }
            t[4] = t[4]+dt;
            Prob[4] = OUPPassageTimeProbability(0,x,t[4],k,omega,rho,mu,sigma);
            if(Prob[4] <= Ppct*PInf)
            {
              m = 4;
              dt = 2*dt;
            }
            else
            {
              m = 3;
              dt = 5;
            }
          }
          else
          {
            for(std::size_t i = 4; i > 0; i--)
            {
              t[i] = t[i-1];
              Prob[i] = Prob[i-1];
            }
            t[0] = t[0]-dt;
            Prob[0] = OUPPassageTimeProbability(0,x,t[0],k,omega,rho,mu,sigma);
            if(Prob[0] >= Ppct*PInf)
            {
              m = 0;
              dt = 2*dt;
            }
            else
            {
              m = 1;
              dt = 5;
            }
          }
        }
        tpct = s+t[2];
      }
      else { tpct = s; }
    }
  }
  return tpct;
}

double OUPPassageTimeMeanIntegrate(double s, double x, double k, double omega, double rho, double mu, double sigma, double median)
{
  double mean = median;
  if(rho < 0.0000000001) { mean = std::numeric_limits<double>::infinity(); }
  else if(median < std::numeric_limits<double>::infinity() && abs(x - k) >= 0.0000000001 && sigma*sigma >= 0.0000000001)
  {
    double PInf = OUPPassageTimeProbabilityInf(x,k,omega,rho,mu,sigma);
    double dss = (median-s)/100.0;
    double ss = 0;
    double dnsty = 0;
    while(dnsty < 0.000001*PInf && ss < median-s)
    {
      ss += dss;
      dnsty = OUPPassageTimeDensity(0,x,ss,k,omega,rho,mu,sigma,0.05);
    }
    ss -= dss;
    std::vector<double> t = Laguerret();
    std::vector<double> w = Laguerrew();
    std::vector<double> wghtfnct(41);
    std::vector<double> density(41);
    double tscale = t[5];
    double wscale = 0;
    mean = 0;
    for(std::size_t i = 0; i < 41; i++)
    {
      wghtfnct[i] = exp(-t[i]);
      t[i] = t[i]*(median-s-ss)/tscale+ss;
      w[i] = w[i]*(median-s-ss)/tscale;
      density[i] = OUPPassageTimeDensity(0,x,t[i],k,omega,rho,mu,sigma,0.05);
      wscale += w[i]*density[i]/wghtfnct[i];
      mean += w[i]*t[i]*density[i]/wghtfnct[i];
    }
    mean = s+mean/wscale;
  }
  return mean;
}

// Exports (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPADrift(z,rho,mu)
//' @param  z   vector of states or optional vector of alternate initial states
//' @param  rho rate parameter 0<=rho<inf
//' @param  mu  location parameter -inf<mu<inf
//' @return g(n) <- RcppOUPADrift()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPADrift(NumericVector z, double rho, double mu)
{
  int n = z.size();
  NumericVector g(n);
  for(int j = 0; j < n; j++) { g[j] = -rho*(z[j]-mu); }

  return g;
}

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPADiffusion(z,sigma)
//' @param  z     vector of states or optional vector of alternate initial states
//' @param  sigma scale parameter -inf<sigma<inf
//' @return h2(n) <- RcppOUPADiffusion()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPADiffusion(NumericVector z, double sigma)
{
  int n = z.size();
  NumericVector h2(n);
  for(int j = 0; j < n; j++) { h2[j] = sigma*sigma; }

  return h2;
}

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAMean(t,s,x,rho,mu,eps)
//' @param  t   terminal time or vector of forward times
//' @param  s   initial time or vector of backward times
//' @param  x   initial state or vector of backward states
//' @param  rho rate parameter 0<=rho<inf
//' @param  mu  location parameter -inf<mu<inf
//' @param  eps proportion remaining after convergence 0<=eps<=1
//' @return Gt(m+1) <- RcppOUPAMean()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPAMean(NumericVector t, double s, double x, double rho, double mu, double eps)
{
  int m = t.size();
  NumericVector Gt(m+1);
  for(int i = 0; i < m; i++) { Gt[i] = mu+(x-mu)*exp(-rho*(t[i]-s)); }
  Gt[m] = s;
  if(eps <= 0 || rho < 0.0000000001) { Gt[m] = std::numeric_limits<double>::infinity(); }
  else { Gt[m] = s-log(eps)/rho; }

  return Gt;
}

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAVariance(t,s,rho,sigma,eps)
//' @param  t     terminal time or vector of forward times
//' @param  s     initial time or vector of backward times
//' @param  rho   rate parameter 0<=rho<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  eps   proportion remaining after convergence 0<=eps<=1
//' @return H2t(m+1) <- RcppOUPAVariance()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPAVariance(NumericVector t, double s, double rho, double sigma, double eps)
{
  int m = t.size();
  NumericVector H2t(m+1);
  if(rho < 0.0000000001)
  {
    for(int i = 0; i < m; i++) { H2t[i] = sigma*sigma*(t[i]-s); }
  }
  else
  {
    for(int i = 0; i < m; i++) { H2t[i] = sigma*sigma/(2*rho)*(1-exp(-2*rho*(t[i]-s))); }
  }
  H2t[m] = s;
  if(eps <= 0 || rho < 0.0000000001) { H2t[m] = std::numeric_limits<double>::infinity(); }
  else { H2t[m] = s-log(eps)/(2*rho); }

  return H2t;
}

#ifdef USE_PARALLEL
struct ROAD : public Worker
{
  RMatrix<double> p;
  const RVector<double> t;
  const RVector<double> y;
  double s;
  double x;
  double rho;
  double mu;
  double sigma;
  double dy;
  std::size_t m;

  ROAD(NumericMatrix& p, const NumericVector& t, const NumericVector& y, double s, double x, double rho, double mu, double sigma, double dy, std::size_t m)
    : p(p), t(t), y(y), s(s), x(x), rho(rho), mu(mu), sigma(sigma), dy(dy), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      p(i,j) = OUPDensity(t[i],y[j],s,x,rho,mu,sigma,dy);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPADensity(t,y,s,x,rho,mu,sigma)
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @return p(m,n) <- RcppOUPADensity()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPADensity(NumericVector t, NumericVector y, double s, double x, double rho, double mu, double sigma)
{
  std::size_t m = t.size();
  std::size_t n = y.size();
  double dy = 0.1;
  if(n > 1) { dy = (y[n-1]-y[0])/(n-1); }
  NumericMatrix p(m,n);
#ifdef USE_PARALLEL
  ROAD worker(p,t,y,s,x,rho,mu,sigma,dy,m);
  parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++)
    {
        p(i,j) = OUPDensity(t[i],y[j],s,x,rho,mu,sigma,dy);
    }
  }
#endif
  return p;
}

#ifdef USE_PARALLEL
struct ROAP : public Worker
{
  RMatrix<double> P;
  const RVector<double> t;
  const RVector<double> y;
  double s;
  double x;
  double rho;
  double mu;
  double sigma;
  double psi;
  std::size_t m;

  ROAP(NumericMatrix& P, const NumericVector& t, const NumericVector& y, double s, double x, double rho, double mu, double sigma, double psi, std::size_t m)
    : P(P), t(t), y(y), s(s), x(x), rho(rho), mu(mu), sigma(sigma), psi(psi), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      P(i,j) = OUPProbability(t[i],y[j],s,x,rho,mu,sigma,psi);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAProbability(t,y,s,x,rho,mu,sigma,psi)
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  psi   <=0 for integral -inf to y, >0 for integral y to inf
//' @return P(m,n) <- RcppOUPAProbability()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAProbability(NumericVector t, NumericVector y, double s, double x, double rho, double mu, double sigma, double psi)
{
  std::size_t m = t.size();
  std::size_t n = y.size();
  NumericMatrix P(m,n);
#ifdef USE_PARALLEL
  ROAP worker(P,t,y,s,x,rho,mu,sigma,psi,m);
  parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++)
    {
      P(i,j) = OUPProbability(t[i],y[j],s,x,rho,mu,sigma,psi);
    }
  }
#endif
  return P;
}

#ifdef USE_PARALLEL
struct ROAPP : public Worker
{
  RMatrix<double> PP;
  const RVector<double> t;
  const RVector<double> y;
  double s;
  double x;
  double rho;
  double mu;
  double sigma;
  double dy;
  double psi;
  std::size_t m;

  ROAPP(NumericMatrix& PP, const NumericVector& t, const NumericVector& y, double s, double x, double rho, double mu, double sigma, double dy, double psi, std::size_t m)
    : PP(PP), t(t), y(y), s(s), x(x), rho(rho), mu(mu), sigma(sigma), dy(dy), psi(psi), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      PP(i,j) = OUPDoubleIntegral(t[i],y[j],s,x,rho,mu,sigma,dy,psi);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPADoubleIntegral(t,y,s,x,rho,mu,sigma,psi)
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  psi   <=0 for integral -inf to y, >0 for integral y to inf
//' @return PP(m,n) <- RcppOUPADoubleIntegral()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPADoubleIntegral(NumericVector t, NumericVector y, double s, double x, double rho, double mu, double sigma, double psi)
{
  std::size_t m = t.size();
  std::size_t n = y.size();
  double dy = 0.1;
  if(n > 1) { dy = (y[n-1]-y[0])/(n-1); }
  NumericMatrix PP(m,n);
#ifdef USE_PARALLEL
  ROAPP worker(PP,t,y,s,x,rho,mu,sigma,dy,psi,m);
  parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++)
    {
      PP(i,j) = OUPDoubleIntegral(t[i],y[j],s,x,rho,mu,sigma,dy,psi);
    }
  }
#endif
  return PP;
}

#ifdef USE_PARALLEL
struct ROAOO : public Worker
{
  RMatrix<double> OO;
  const RVector<double> s;
  const RVector<double> x;
  double t;
  double y;
  double rho;
  double mu;
  double sigma;
  double dy;
  double r;
  double phi;
  double b;
  double c;
  std::size_t m;

  ROAOO(NumericMatrix& OO, const NumericVector& s, const NumericVector& x, double t, double y, double rho, double mu, double sigma, double dy, double r, double phi, double b, double c, std::size_t m)
    : OO(OO), s(s), x(x), t(t), y(y), rho(rho), mu(mu), sigma(sigma), dy(dy), r(r), phi(phi), b(b), c(c), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      OO(i,j) =  OUPOption(s[i],x[j],t,y,rho,mu,sigma,dy,r,phi,b,c);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAOption(s,x,t,y,rho,mu,sigma,r,phi,b,c)
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  r     discount rate -inf<r<inf
//' @param  phi   <=0 for integral -inf to x, >0 for integral x to inf
//' @param  b     lump-sum benefit for entry option
//' @param  c     lump-sum cost for exit option
//' @return OO(m,n) <- RcppOUPAOption()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAOption(NumericVector s, NumericVector x, double t, double y, double rho, double mu, double sigma, double r, double phi, double b, double c)
{
  std::size_t m = s.size();
  std::size_t n = x.size();
  double dy = 0.1;
  if(n > 1) { dy = (x[n-1]-x[0])/(n-1); }
  NumericMatrix OO(m,n);
#ifdef USE_PARALLEL
  ROAOO worker(OO,s,x,t,y,rho,mu,sigma,dy,r,phi,b,c,m);
  parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++)
    {
      OO(i,j) = OUPOption(s[i],x[j],t,y,rho,mu,sigma,dy,r,phi,b,c);
    }
  }
#endif
  return OO;
}

#ifdef USE_PARALLEL
struct ROAOOs : public Worker
{
  RMatrix<double> OOs;
  const RVector<double> x;
  double t;
  double y;
  double rho;
  double mu;
  double sigma;
  double dy;
  double r;
  double phi;
  double b;
  double c;
  double tsguess;
  double tsmax;
  double tsmin;
  double maxmin;

  ROAOOs(NumericMatrix& OOs, const NumericVector& x,double t, double y, double rho, double mu, double sigma, double dy, double r, double phi, double b, double c, double tsguess, double tsmax, double tsmin, double maxmin)
    : OOs(OOs), x(x), t(t), y(y), rho(rho), mu(mu), sigma(sigma), dy(dy), r(r), phi(phi), b(b), c(c), tsguess(tsguess), tsmax(tsmax), tsmin(tsmin), maxmin(maxmin) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      std::vector<double> OOshat = OUPOptionMaxMin(x[j],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,maxmin);
      OOs(0,j) = OOshat[0];
      OOs(1,j) = t-OOshat[1];
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAOptionEnvelope(s,x,t,y,rho,mu,sigma,r,phi,b,c)
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  r     discount rate -inf<r<inf
//' @param  phi   <=0 for integral -inf to x, >0 for integral x to inf
//' @param  b     lump-sum benefit for entry option
//' @param  c     lump-sum cost for exit option
//' @return OOs(2,n) <- RcppOUPAOptionEnvelope()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAOptionEnvelope(NumericVector s, NumericVector x, double t, double y, double rho, double mu, double sigma, double r, double phi, double b, double c)
{
  std::size_t m = s.size();
  std::size_t n = x.size();
  double dy = 0.1;
  if(n > 1) { dy = (x[n-1]-x[0])/(n-1); }
  NumericMatrix OOs(2,n);
  double tsmax = s[0];
  double tsmin = s[m-1];
  double tsguess = tsmin;
  tsguess = OUPOptionMaxMin(y,y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1)[1];
#ifdef USE_PARALLEL
  ROAOOs worker(OOs,x,y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1);
  parallelFor(0,n,worker);
#else
  if(phi > 0)
  {
    std::size_t j = 0;
    while(j < n)
    {
      std::vector<double> OOshat = OUPOptionMaxMin(x[j],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1);
      OOs(0,j) = OOshat[0];
      OOs(1,j) = t-OOshat[1];
      tsguess = OOshat[1];
      j += 1;
    }
  }
  else
  {
    std::size_t j = n;
    while(j > 0)
    {
      j -= 1;
      std::vector<double> OOshat = OUPOptionMaxMin(x[j],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1);
      OOs(0,j) = OOshat[0];
      OOs(1,j) = t-OOshat[1];
      tsguess = OOshat[1];
    }
  }
#endif
  return OOs;
}

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAdOOdsZero(s,x,t,y,rho,mu,sigma,r,phi,b,c)
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  r     discount rate -inf<r<inf
//' @param  phi   <=0 for integral -inf to x, >0 for integral x to inf
//' @param  b     lump-sum benefit for entry option
//' @param  c     lump-sum cost for exit option
//' @return dOOdszero(4,n+3) <- RcppOUPAdOOdsZero()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAdOOdsZero(NumericVector s, NumericVector x, double t, double y, double rho, double mu, double sigma, double r, double phi, double b, double c)
{
  int m = s.size();
  int n = x.size();
  double dy = 0.1;
  if(n > 1) { dy = (x[n-1]-x[0])/(n-1); }
  NumericMatrix dOOdszero(4,n+3);
  if(phi > 0 )
  {
    std::vector<double> env(2);
    double tsmax = s[0];
    double tsmin = s[m-1];
    double tsguess = tsmax;
    int j = -1;
    while(tsguess > tsmin && j < n-1)
    {
      j += 1;
      env = OUPOptionMaxMin(x[j],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1);
      dOOdszero(2,j) = env[0];
      dOOdszero(3,j) = t-env[1];
      tsguess = env[1];
    }
    int k = j;
    if(tsguess <= tsmin)
    {
      while(j < n)
      {
        dOOdszero(0,j) = NA_REAL;
        dOOdszero(1,j) = t-tsmax;
        dOOdszero(2,j) = NA_REAL;
        dOOdszero(3,j) = t-tsmin;
        j += 1;
      }
    }
    if(k > 0) { tsmax = t-dOOdszero(3,k-1); }
    tsguess = tsmax;
    int i = k;
    while(x[i] > y && i > 0)
    {
      i -= 1;
      env = OUPOptionMaxMin(x[i],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,-1);
      dOOdszero(0,i) = env[0];
      dOOdszero(1,i) = t-env[1];
      tsguess = env[1];
    }
    while(i > 0)
    {
      i -= 1;
      dOOdszero(0,i) = 0.0;
      dOOdszero(1,i) = t-tsmin;
    }
    if(k < n && k > 0)
    {
      dOOdszero(2,n) = x[k-1];
      dOOdszero(2,n+1) = x[k-1];
      dOOdszero(2,n+2) = x[k-1];
      dOOdszero(1,n) = dOOdszero(1,k-1);
      dOOdszero(1,n+1) = (dOOdszero(1,k-1)+dOOdszero(3,k-1))/2.0;
      dOOdszero(1,n+2) = dOOdszero(3,k-1);
      dOOdszero(0,n) = dOOdszero(0,k-1);
      dOOdszero(0,n+1) = OUPOption(-t+dOOdszero(1,n+1),x[k-1],0,y,rho,mu,sigma,dy,r,phi,b,c);
      dOOdszero(0,n+2) = dOOdszero(2,k-1);
    }
  }
  else
  {
    std::vector<double> env(2);
    double tsmax = s[0];
    double tsmin = s[m-1];
    double tsguess = tsmax;
    int j = n;
    while(tsguess > tsmin && j > 0)
    {
      j -= 1;
      env = OUPOptionMaxMin(x[j],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,1);
      dOOdszero(2,j) = env[0];
      dOOdszero(3,j) = t-env[1];
      tsguess = env[1];
    }
    int k = j;
    if(tsguess <= tsmin)
    {
      while(j > -1)
      {
        dOOdszero(0,j) = NA_REAL;
        dOOdszero(1,j) = t-tsmax;
        dOOdszero(2,j) = NA_REAL;
        dOOdszero(3,j) = t-tsmin;
        j -= 1;
      }
    }
    if(k < n) { tsmax = t-dOOdszero(3,k+1); }
    tsguess = tsmax;
    int i = k;
    while(x[i] < y && i < n)
    {
      i += 1;
      env = OUPOptionMaxMin(x[i],y,rho,mu,sigma,dy,r,phi,b,c,tsguess,tsmax,tsmin,-1);
      dOOdszero(0,i) = env[0];
      dOOdszero(1,i) = t-env[1];
      tsguess = env[1];
    }
    while(i < n-1)
    {
      i += 1;
      dOOdszero(0,i) = 0.0;
      dOOdszero(1,i) = t-tsmin;
    }
    if(k < n-1 && k > -1)
    {
      dOOdszero(2,n) = x[k+1];
      dOOdszero(2,n+1) = x[k+1];
      dOOdszero(2,n+2) = x[k+1];
      dOOdszero(1,n) = dOOdszero(1,k+1);
      dOOdszero(1,n+1) = (dOOdszero(1,k+1)+dOOdszero(3,k+1))/2.0;
      dOOdszero(1,n+2) = dOOdszero(3,k+1);
      dOOdszero(0,n) = dOOdszero(0,k+1);
      dOOdszero(0,n+1) = OUPOption(-t+dOOdszero(1,n+1),x[k+1],0,y,rho,mu,sigma,dy,r,phi,b,c);
      dOOdszero(0,n+2) = dOOdszero(2,k+1);
    }
  }
  return dOOdszero;
}

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPADecisionThreshold(y,rho,mu,sigma,r,phi,b,c)
//' @param  y     terminal state or vector of forward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  r     discount rate -inf<r<inf
//' @param  phi   <=0 for integral -inf to x, >0 for integral x to inf
//' @param  b     lump-sum benefit for entry option
//' @param  c     lump-sum cost for exit option
//' @return kOO(2) <- RcppOUPADecisionThreshold()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPADecisionThreshold(double y, double rho, double mu, double sigma, double r, double phi, double b, double c)
{
  NumericVector kOO(2);
  if(rho < 0.0000000001 && r < 0.0000000001)
  {
    if(phi > 0) { kOO[0] = std::numeric_limits<double>::infinity(); }
    else { kOO[0] = -std::numeric_limits<double>::infinity(); }
    kOO[1] = std::numeric_limits<double>::infinity();
  }
  else if(sigma*sigma < 0.0000000001 && rho < 0.0000000001)
  {
    kOO[0] = y;
    if(phi > 0) { kOO[1] = b; }
    else { kOO[1] = -c; }
  }
  else
  {
    double tensigdig = 10000000000000;
    int n = 1000;
    std::vector<double> x(5);
    std::vector<double> t(5);
    std::vector<double> Optn(5);
    std::vector<double> env(2);
    double tsguess = 0;
    int ltgt = 1;
    if(phi > 0) { ltgt = -1; }
    double dx;
    int m;
    if(ltgt < 0)
    {
      x[0] = y;
      env = OUPOptionMaxMin(x[0],y,rho,mu,sigma,0.1,r,phi,b,c,tsguess,1100,0,1);
      Optn[0] = env[0];
      t[0] = env[1];
      dx = (Optn[0]-b)/2;
      m = 0;
      for(int i = 1; i < 5; i++)
      {
        x[i] = x[i-1]+dx;
        env = OUPOptionMaxMin(x[i],y,rho,mu,sigma,0.1,r,phi,b,c,env[1],1100,0,1);
        Optn[i] = env[0];
        t[i] = env[1];
        if(Optn[i] > ltgt*(y-x[i])+b) { m = i; }
      }
    }
    else
    {
      x[4] = y;
      env = OUPOptionMaxMin(x[4],y,rho,mu,sigma,0.1,r,phi,b,c,tsguess,1100,0,1);
      Optn[4] = env[0];
      t[4] = env[1];
      dx = (Optn[4]+c)/2;
      m = 4;
      for(int i = 4; i > 0; i--)
      {
        x[i-1] = x[i]-dx;
        env = OUPOptionMaxMin(x[i-1],y,rho,mu,sigma,0.1,r,phi,b,c,env[1],1100,0,1);
        Optn[i-1] = env[0];
        t[i-1] = env[1];
        if(Optn[i-1] > ltgt*(y-x[i-1])-c) { m = i-1; }
      }
    }
    int j = 0;
    while(tensigdig*(x[3]-x[1]) >= abs(x[2]) && j < n)
    {
      j += 1;
      if(m == 1 || m == 2 || m == 3)
      {
        if(m == 1)
        {
          x[4] = x[2];
          Optn[4] = Optn[2];
          x[2] = x[1];
          Optn[2] = Optn[1];
        }
        else if(m == 2)
        {
          x[0] = x[1];
          Optn[0] = Optn[1];
          x[4] = x[3];
          Optn[4] = Optn[3];
        }
        else
        {
          x[0] = x[2];
          Optn[0] = Optn[2];
          x[2] = x[3];
          Optn[2] = Optn[3];
        }
        x[1] = (x[2]+x[0])/2;
        tsguess = (t[2]+t[0])/2;
        env = OUPOptionMaxMin(x[1],y,rho,mu,sigma,0.1,r,phi,b,c,tsguess,1100,0,1);
        Optn[1] = env[0];
        t[1] = env[1];
        x[3] = (x[2]+x[4])/2;
        tsguess = (t[2]+t[4])/2;
        env = OUPOptionMaxMin(x[3],y,rho,mu,sigma,0.1,r,phi,b,c,tsguess,1100,0,1);
        Optn[3] = env[0];
        t[3] = env[1];
        dx *= 0.5;
        if(ltgt < 0)
        {
          if(Optn[3] > ltgt*(y-x[3])+b) { m = 3; }
          else if(Optn[2] > ltgt*(y-x[2])+b) { m = 2; }
          else { m = 1; }
        }
        else
        {
          if(Optn[1] > ltgt*(y-x[1])-c) { m = 1; }
          else if(Optn[2] > ltgt*(y-x[2])-c) { m = 2; }
          else { m = 3; }
        }
      }
      else if(m == 4)
      {
        for(int i = 1; i < 5; i++)
        {
          x[i-1] = x[i];
          Optn[i-1] = Optn[i];
        }
        x[4] = x[4]+dx;
        env = OUPOptionMaxMin(x[4],y,rho,mu,sigma,0.1,r,phi,b,c,t[4],1100,0,1);
        Optn[4] = env[0];
        t[4] = env[1];
        if(ltgt < 0)
        {
          if(Optn[4] > ltgt*(y-x[4])+b)
          {
            m = 4;
            dx *= 2;
          }
          else { m = 3; }
        }
        else
        {
          if(Optn[4] > ltgt*(y-x[4])-c) { m = 3; }
          else
          {
            m = 4;
            dx *= 2;
          }
        }
      }
      else
      {
        for(int i = 4; i > 0; i--)
        {
          x[i] = x[i-1];
          Optn[i] = Optn[i-1];
        }
        x[0] = x[0]-dx;
        env = OUPOptionMaxMin(x[0],y,rho,mu,sigma,0.1,r,phi,b,c,t[0],1100,0,1);
        Optn[0] = env[0];
        t[0] = env[1];
        if(ltgt < 0)
        {
          if(Optn[0] > ltgt*(y-x[0])+b) { m = 1; }
          else
          {
            m = 0;
            dx *= 2;
          }
        }
        else
        {
          if(Optn[0] > ltgt*(y-x[0])-c)
          {
            m = 0;
            dx *= 2;
          }
          else { m = 1; }
        }
      }
      kOO[0] = x[2];
      if(ltgt < 0) { kOO[1] = ltgt*(y-x[2])+b; }
      else { kOO[1] = ltgt*(y-x[2])-c; }
    }
  }
  return kOO;
}

#ifdef USE_PARALLEL
struct ROABC : public Worker
{
  RMatrix<double> BC;
  const RVector<double> s;
  const RVector<double> x;
  double t;
  double y;
  double rho;
  double mu;
  double r;
  double phi;
  double b;
  double c;
  std::size_t m;

  ROABC(NumericMatrix& BC, const NumericVector& s, const NumericVector& x, double t, double y, double rho, double mu, double r, double phi, double b, double c, std::size_t m)
    : BC(BC), s(s), x(x), t(t), y(y), rho(rho), mu(mu), r(r), phi(phi), b(b), c(c), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      BC(i,j) =  OUPObligation(s[i],x[j],t,y,rho,mu,r,phi,b,c);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAObligation(s,x,t,y,rho,mu,r,phi,b,c)
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  t     terminal time or vector of forward times
//' @param  y     terminal state or vector of forward states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  r     discount rate -inf<r<inf
//' @param  phi   <=0 for integral -inf to x, >0 for integral x to inf
//' @param  b     lump-sum benefit for entry option
//' @param  c     lump-sum cost for exit option
//' @return BC(m,n) <- RcppOUPAObligation()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAObligation(NumericVector s, NumericVector x, double t, double y, double rho, double mu, double r, double phi, double b, double c)
{
  std::size_t m = s.size();
  std::size_t n = x.size();
  NumericMatrix BC(m,n);
#ifdef USE_PARALLEL
  ROABC worker(BC,s,x,t,y,rho,mu,r,phi,b,c,m);
  parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++)
    {
      BC(i,j) = OUPObligation(s[i],x[j],t,y,rho,mu,r,phi,b,c);
    }
  }
#endif
  return BC;
}

#ifdef USE_PARALLEL
struct ROAPTMMMzx : public Worker
{
  RMatrix<double> tmmmtmmmx;
  const RVector<double> zx;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;

  ROAPTMMMzx(NumericMatrix& tmmmtmmmx, const NumericVector& zx, double s, double k, double omega, double rho, double mu, double sigma)
    : tmmmtmmmx(tmmmtmmmx), zx(zx), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      tmmmtmmmx(1,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,0.5);
      tmmmtmmmx(0,j) = OUPPassageTimeModeSearch(s,zx[j],k,omega,rho,mu,sigma,tmmmtmmmx(1,j));
      tmmmtmmmx(2,j) = OUPPassageTimeMeanIntegrate(s,zx[j],k,omega,rho,mu,sigma,tmmmtmmmx(1,j));
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAPassageTimeModeMedianMean(k,s,x,omega,rho,mu,sigma,z)
//' @param  k     decision threshold
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  omega degree of irreversibility 0<=omega<=1
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  z     vector of states or optional vector of alternate initial states
//' @return tmmmtmmmx(3,n+1) <- RcppOUPAPassageTimeModeMedianMean()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAPassageTimeModeMedianMean(double k, double s, double x, double omega, double rho, double mu, double sigma, Nullable<NumericVector> z=R_NilValue)
{
  if(z.isNull())
  {
    NumericMatrix tmmmtmmmx(3,1);
    tmmmtmmmx(1,0) = OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,0.5);
    tmmmtmmmx(0,0) = OUPPassageTimeModeSearch(s,x,k,omega,rho,mu,sigma,tmmmtmmmx(1,0));
    tmmmtmmmx(2,0) = OUPPassageTimeMeanIntegrate(s,x,k,omega,rho,mu,sigma,tmmmtmmmx(1,0));
    return tmmmtmmmx;
  }
  else
  {
    NumericVector zx(z);
    zx.push_back(x);
    std::size_t n = zx.size();
    NumericMatrix tmmmtmmmx(3,n);
#ifdef USE_PARALLEL
    ROAPTMMMzx worker(tmmmtmmmx,zx,s,k,omega,rho,mu,sigma);
    parallelFor(0,n,worker);
#else
    for(std::size_t j = 0; j < n; j++)
    {
      tmmmtmmmx(1,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,0.5);
      tmmmtmmmx(0,j) = OUPPassageTimeModeSearch(s,zx[j],k,omega,rho,mu,sigma,tmmmtmmmx(1,j));
      tmmmtmmmx(2,j) = OUPPassageTimeMeanIntegrate(s,zx[j],k,omega,rho,mu,sigma,tmmmtmmmx(1,j));
    }
#endif
    return tmmmtmmmx;
  }
}

#ifdef USE_PARALLEL
struct ROAPTPctzx : public Worker
{
  RMatrix<double> tpcttpctx;
  const RVector<double> zx;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;
  double Ppct;

  ROAPTPctzx(NumericMatrix& tpcttpctx, const NumericVector& zx, double s, double k, double omega, double rho, double mu, double sigma, double Ppct)
    : tpcttpctx(tpcttpctx), zx(zx), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma), Ppct(Ppct) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      tpcttpctx(0,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,1-Ppct);
      tpcttpctx(1,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,0.5);
      tpcttpctx(2,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,Ppct);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAPassageTimePercentiles(k,s,x,omega,Ppct,rho,mu,sigma,z)
//' @param  k     decision threshold
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  omega degree of irreversibility 0<=omega<=1
//' @param  Ppct  passage time probability for a percentile 0.01<=Ppct<=0.99
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  z     vector of states or optional vector of alternate initial states
//' @return tpcttpctx(3,n+1) <- RcppOUPAPassageTimePercentiles()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAPassageTimePercentiles(double k, double s, double x, double omega, double Ppct, double rho, double mu, double sigma, Nullable<NumericVector> z=R_NilValue)
{
  if(z.isNull())
  {
    NumericMatrix tpcttpctx(3,1);
    double pct = Ppct;
    if(pct < 0.5) { pct = 1-pct; }
    tpcttpctx(0,0) = OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,1-pct);
    tpcttpctx(1,0) = OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,0.5);
    tpcttpctx(2,0) = OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,pct);
    return tpcttpctx;
  }
  else
  {
    NumericVector zx(z);
    zx.push_back(x);
    std::size_t n = zx.size();
    NumericMatrix tpcttpctx(3,n);
    double pct = Ppct;
    if(pct < 0.5) { pct = 1-pct; }
#ifdef USE_PARALLEL
    ROAPTPctzx worker(tpcttpctx,zx,s,k,omega,rho,mu,sigma,pct);
    parallelFor(0,n,worker);
#else
    for(std::size_t j = 0; j < n; j++)
    {
      tpcttpctx(0,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,1-pct);
      tpcttpctx(1,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,0.5);
      tpcttpctx(2,j) = OUPPassageTimePctSearch(s,zx[j],k,omega,rho,mu,sigma,pct);
    }
#endif
    return tpcttpctx;
  }
}

#ifdef USE_PARALLEL
struct ROAPTDx : public Worker
{
  RMatrix<double> ptptx;
  const RVector<double> t;
  double x;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;
  double dt;
  std::size_t m;

  ROAPTDx(NumericMatrix& ptptx, const NumericVector& t, double x, double s, double k, double omega, double rho, double mu, double sigma, double dt, std::size_t m)
    : ptptx(ptptx), t(t), x(x), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma), dt(dt), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) { ptptx(i,0) = OUPPassageTimeDensity(s,x,t[i],k,omega,rho,mu,sigma,dt); }
  }
};

struct ROAPTDzx : public Worker
{
  RMatrix<double> ptptx;
  const RVector<double> t;
  const RVector<double> zx;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;
  double dt;
  std::size_t m;

  ROAPTDzx(NumericMatrix& ptptx, const NumericVector& t, const NumericVector& zx, double s, double k, double omega, double rho, double mu, double sigma, double dt, std::size_t m)
    : ptptx(ptptx), t(t), zx(zx), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma), dt(dt), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      ptptx(i,j) = OUPPassageTimeDensity(s,zx[j],t[i],k,omega,rho,mu,sigma,dt);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAPassageTimeDensity(t,k,s,x,omega,rho,mu,sigma,z)
//' @param  t     terminal time or vector of forward times
//' @param  k     decision threshold
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  omega degree of irreversibility 0<=omega<=1
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  z     vector of states or optional vector of alternate initial states
//' @return ptptx(m,n+1) <- RcppOUPAPassageTimeDensity()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAPassageTimeDensity(NumericVector t, double k, double s, double x, double omega, double rho, double mu, double sigma, Nullable<NumericVector> z=R_NilValue)
{
  std::size_t m = t.size();
  double dt = 0.5;
  if(m > 1) { dt = (t[m-1]-t[0])/(m-1); }
  if(z.isNull())
  {
    NumericMatrix ptptx(m,1);
#ifdef USE_PARALLEL
    ROAPTDx worker(ptptx,t,x,s,k,omega,rho,mu,sigma,dt,m);
    parallelFor(0,m,worker);
#else
    for(std::size_t i = 0; i < m; i++) {  ptptx(i,0) = OUPPassageTimeDensity(s,x,t[i],k,omega,rho,mu,sigma,dt); }
#endif
    return ptptx;
  }
  else
  {
    NumericVector zx(z);
    zx.push_back(x);
    std::size_t n = zx.size();
    NumericMatrix ptptx(m,n);
#ifdef USE_PARALLEL
    ROAPTDzx worker(ptptx,t,zx,s,k,omega,rho,mu,sigma,dt,m);
    parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++) { ptptx(i,j) = OUPPassageTimeDensity(s,zx[j],t[i],k,omega,rho,mu,sigma,dt); }
  }
#endif
    return ptptx;
  }
}

#ifdef USE_PARALLEL
struct ROAPTPx : public Worker
{
  RMatrix<double> PtPtx;
  const RVector<double> t;
  double x;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;
  std::size_t m;

  ROAPTPx(NumericMatrix& PtPtx, const NumericVector& t, double x, double s, double k, double omega, double rho, double mu, double sigma, std::size_t m)
    : PtPtx(PtPtx), t(t), x(x), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++)
    {
      PtPtx(i,0) = OUPPassageTimeProbability(s,x,t[i],k,omega,rho,mu,sigma);
    }
  }
};

struct ROAPTPzx : public Worker
{
  RMatrix<double> PtPtx;
  const RVector<double> t;
  const RVector<double> zx;
  double s;
  double k;
  double omega;
  double rho;
  double mu;
  double sigma;
  std::size_t m;

  ROAPTPzx(NumericMatrix& PtPtx, const NumericVector& t, const NumericVector& zx, double s, double k, double omega, double rho, double mu, double sigma, std::size_t m)
    : PtPtx(PtPtx), t(t), zx(zx), s(s), k(k), omega(omega), rho(rho), mu(mu), sigma(sigma), m(m) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t ij = begin; ij < end; ij++)
    {
      std::size_t i = ij%m;
      std::size_t j = ij/m;
      PtPtx(i,j) = OUPPassageTimeProbability(s,zx[j],t[i],k,omega,rho,mu,sigma);
    }
  }
};
#endif

//' @rdname Analytical_Rcpp
//' @usage  RcppOUPAPassageTimeProbability(t,k,s,x,omega,rho,mu,sigma,z)
//' @param  t     terminal time or vector of forward times
//' @param  k     decision threshold
//' @param  s     initial time or vector of backward times
//' @param  x     initial state or vector of backward states
//' @param  omega degree of irreversibility 0<=omega<=1
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @param  z     vector of states or optional vector of alternate initial states
//' @return PtPtx(m,n+1) <- RcppOUPAPassageTimeProbability()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPAPassageTimeProbability(NumericVector t, double k, double s, double x, double omega, double rho, double mu, double sigma, Nullable<NumericVector> z=R_NilValue)
{
  std::size_t m = t.size();
  if(z.isNull())
  {
    NumericMatrix PtPtx(m,1);
#ifdef USE_PARALLEL
    ROAPTPx worker(PtPtx,t,x,s,k,omega,rho,mu,sigma,m);
    parallelFor(0,m,worker);
#else
    for(std::size_t i = 0; i < m; i++) { PtPtx(i,0) = OUPPassageTimeProbability(s,x,t[i],k,omega,rho,mu,sigma); }
#endif
    return PtPtx;
  }
  else
  {
    NumericVector zx(z);
    zx.push_back(x);
    std::size_t n = zx.size();
    NumericMatrix PtPtx(m,n);
#ifdef USE_PARALLEL
    ROAPTPzx worker(PtPtx,t,zx,s,k,omega,rho,mu,sigma,m);
    parallelFor(0,m*n,worker);
#else
  for(std::size_t i = 0; i < m; i++)
  {
    for(std::size_t j = 0; j < n; j++) { PtPtx(i,j) = OUPPassageTimeProbability(s,zx[j],t[i],k,omega,rho,mu,sigma); }
  }
#endif
    return PtPtx;
  }
}
