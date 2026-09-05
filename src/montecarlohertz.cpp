#include <Rcpp.h>
using namespace Rcpp;
#ifdef USE_PARALLEL
#include <RcppParallel.h>
using namespace RcppParallel;
#ifdef USE_SITMO
#include <sitmo.h>
#endif
#endif
#include <random>
#include <cmath>
#include <limits>
#include "gammahertz.h"

// roxygen (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @title MonteCarlo_Rcpp functions for simulating an Ornstein-Uhlenbeck Process
//'
//' @description
//' Calculations for the R6 class 'MonteCarlo', with parallel processing.
//'
//' @details # Notes on Values
//' Return values are vectors and matrices allocated in Rcpp.  The dimensions are
//'  shown for information.  Of course, do not include them in R calls.  For example:
//'
//'     stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
//'
//' The return values:
//'
//'     stdnorm((m-1)*skip,paths)
//'     forward(m,paths)
//'     backward(m,paths)
//'     heat(m,n)
//'
//'  are matrices of pseudo-random standard normal variables, forward paths,
//'  backward paths and heat maps.
//'
//' The return value:
//'
//'     minmax(2)
//'
//'  is a vector of minimum and maximum values, subset in R as:
//'
//'     minmax <- RcppOUPMCMinMax(matPaths)
//'     min <- minmax[1]
//'     max <- minmax[2]
//'
//' The return value:
//'
//'     bndfpt(m+1,paths)
//'
//'  is a matrix of bounded paths with a vector of first passage times for each
//'  path in row m+1.  There are NA entries for paths which previously hit the
//'  threshold and NA entries in row m+1 for paths which have yet to hit the
//'  threshold.  Subset in R as:
//'
//'     bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma)
//'     bounded <- bndfpt[1:m,,drop=FALSE]
//'     fpt <- bndfpt[m+1,,drop=FALSE]
//'
//' The return value:
//'
//'     mvdpd(m,3*n+2)
//'
//'  is a composite matrix containing two vectors for means and variances and three
//'  matrices for densities, probabilities and double integrals.  The vectors and
//'  matrices are subset in R as:
//'
//'     mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
//'     means <- mvdpd[,1,drop=FALSE]
//'     variances <- mvdpd[,2,drop=FALSE]
//'     densities <- mvdpd[,3:(n+2),drop=FALSE]
//'     probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
//'     doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
//'
//' Similarly, the return value:
//'
//'     dpo(m,3*n)
//'
//'  is a composite of three contiguous matrices for prior densities, prior
//'  probabilities and options, subset in R as:
//'
//'     dpo <- RcppOUPMCBackwardCountX(backward,x,phi,rho,r,ds)
//'     densities <- dpo[,1:n,drop=FALSE]
//'     probabilities <- dpo[,(n+1):(2*n),drop=FALSE]
//'     options <- dpo[,(2*n+1):(3*n),drop=FALSE]
//'
//' The return values:
//'
//'     pctdp(m,5)
//'
//'  are matrices of five columns.  The first column has only five entries for
//'  times at the mode, median, mean, lower percentile, and upper percentile. The
//'  second and third columns each have five entries containing the corresponding
//'  densities and probabilities.  The fourth and fifth columns contain contain m
//'  entries for either visiting time or first passage time densities and
//'  probabilities.  The columns are subset in R as:
//'
//'     pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
//'     mode <- pctdp[1,1:3,drop=FALSE]
//'     median <- pctdp[2,1:3,drop=FALSE]
//'     mean <- pctdp[3,1:3,drop=FALSE]
//'     lowerpct <- pctdp[4,1:3,drop=FALSE]
//'     upperpct <- pctdp[5,1:3,drop=FALSE]
//'     densities <- pctdp[,4,drop=FALSE]
//'     probabilities <- pctdp[,5,drop=FALSE]
//'
//' @details # Discussion
//' A single-threaded R6 object is fast enough for many calculations, but not
//'  for Monte-Carlo simulations.  Attempts at parallel processing using parApply()
//'  and future_apply() failed.  The whole R6 object is copied to each thread,
//'  which locks up the computer.  Rccp is often hundreds of times faster and makes
//'  Monte-Carlo simulations practical for interactive applications such as RStudio
//'  and RShiny. RcppParallel speeds the calculations another five to eight times.
//'
//' For the simulations, both a 4th order Runge-Kutta method and the stochastic
//'  integral equation are implemented.  With argument skip set to 10, and the
//'  same seed in the random number generators, they give the same paths to about six
//'  significant digits, but the stochastic integral equation calculates faster.
//'
//' A typical simulation might require 100,000,000 standard normal variables.  On
//'  an i7 processor with 12 threads running at a maximum speed of 4.5 GHz,
//'  microbenchmark median times to generate the standard normal variables are:
//'
//'     Unit: milliseconds             R6       R6+           R6+
//'               function  single-thread      Rcpp  RcppParallel
//'     ---------------------------------------------------------
//'         StandardNormal       5779.756  4002.566       660.010
//'
//' R6 single-thread and R6+Rcpp use rnorm(), the standard random number
//'  generator in R.  The R6+RcppParallel uses sitmo::prng_engine with a Box-Muller
//'  transform for uniform to normal random variables.  R6+Rcpp is 1.4 times
//'  faster than R6 single-thread.  R6+RcppParallel is 6.1 times faster than
//'  R6+Rcpp and 8.8 times faster than R6 single-thread.
//'
//' Once the standard normal variables are calculated, microbenchmark median
//'  times to simulate the 4th order Runge-Kutta and stochastic integral
//'  equation for 100,000 paths over 100 time intervals with skip=10 are:
//'
//'     Unit: milliseconds                      R6       R6+           R6+
//'                        function  single-thread      Rcpp  RcppParallel
//'     ------------------------------------------------------------------
//'           ForwardPathRungeKutta      207798.80  2688.932       242.135
//'     ForwardPathIntegralEquation       51298.03   416.793        47.316
//'
//' The stochastic integral equation calculates from 4.1 to 6.5 times faster than
//'  the 4th order Runge-Kutta method.  R6+Rcpp calculates from 77.3 to 123.1 times
//'  faster than R6 single-thread.  R6+RcppParallel calculates from 8.8 to 11.1
//'  times faster than R6+Rcpp and from 858.2 to 1084.2 times faster than R6
//'  single-thread.
//'
//' The skip parameter can increase the accuracy of the Runge-Kutta method and
//'  enable a better count of First Passage Times.  The stochastic integral
//'  equation is not improved by skip=10 and is penalized in the timings above.
//'
//' Forward, backward and bounded paths can be binned and counted to approximate
//'  solutions. The approximations converge to analytical solutions as the number
//'  of paths increases.  Binning and counting 1,000,000 paths will be accurate to
//'  3 or 4 significant digits.  Here are microbenchmark median times for the
//'  stochastic integral equation over 100,000 and 1,000,000 paths for 100 time
//'  intervals with skip=1, as calculated by R6+RccpParallel:
//'
//'     Unit: milliseconds              paths                paths
//'                        function   100,000            1,000,000
//'     ----------------------------------------------------------
//'                  StandardNormal   57.2060             588.7932
//'     ForwardPathIntegralEquation   18.9317  ________   187.9792  _________
//'                        subtotal             76.1377              776.7524
//'                   ForwardCountY   70.0579  ________   906.2226  _________
//'                           total            146.1956             1682.9750
//'
//' Times go up approximately linearly with the number of paths.  Looking at the
//'  subtotals for simulating from a standing start, 100,000 paths will take
//'  0.076 seconds and 1,000,000 paths will take 0.777 seconds.  Looking at the
//'  totals for counting from a standing start, 100,000 paths will take 0.146
//'  seconds and 1,000,000 paths will take 1.683 seconds.  About a third of that
//'  time is spent generating the standard normal variables, which can be reused.
//'  Subsequent simulations will only take 0.019 and 0.188 seconds.  Forward Paths
//'  are not reused and subsequent simulations plus counting will take 0.089 and
//'  1.094 seconds.
//'
//' The function ForwardCountY bins and counts means, variances, transition densities,
//'  transition probabilities and double integrals.  So five sets of plots can be
//'  drawn from one set of calculations. For comparison, a 3D plot by Plotly will
//'  take up to a second on an RTX 2070 GPU.  So calculations are only a part of
//'  the job.
//'
//' Rcpp versions of the functions were coded first.  All but one function were
//'  translated into RcppParallel versions.  RccpParallel uses Intel's Threading
//'  Building Blocks (TBB) on the CPU.  Unlike parallel processing on a GPU or
//'  accelerator, memory isn't copied and there is less overhead.  On trivially
//'  small problems, sequential versions calculate faster.  On large problems,
//'  parallel versions calculate much faster.
//'
//' RcppParallel and sitmo are optional packages.  If installed, they will be
//'  used.  Function RcppParallelInstalled() will enquire whether calculations
//'  will use RcppParallel or fall back to Rcpp only.  Function RcppsitmoInstalled()
//'  will enquire whether random numbers will be generated by RcppPrallel with
//'  sitmo() or fall back to Rcpp with rnorm().
//'
//' @details # From the Console
//' These functions are available in R, the RStudio console and RShiny apps.
//'  For example, a simulation of 1,000,000 forward paths over 100 time intervals
//'  with skip=1 would be:
//'
//'      stdnorm <- RcppOUPMCStandardNormal(101,1,1000000,9999)
//'      fwd <- RcppOUPMCForwardPathIntegralEquation(stdnorm,15,100,1,0.1,0.5,-15,15)
//'
//' A microbenchmark comparison of indirectly calling the RcppParallel functions
//'  from R6 with directly calling them from the console is:
//'
//'     Unit: milliseconds                     R6+        Console
//'                         function   RcppParallel  RcppParallel
//'     -------------------------------------------------------------------
//'      ForwardPathIntegralEquation       807.2069      837.4392
//'     BackwardPathIntegralEquation       807.4201      837.8897
//'      BoundedPathIntegralEquation      1200.3780      834.9122
//'
//' These timings are from a standing start, generating the random variables before
//'  simulating the paths.  For Forward and Backward Paths, the R6 object is faster,
//'  but for Bounded Paths, it is much slower.  Bounded Paths hit the boundary and have
//'  NA values thereafter.  We might speculate that R6 is slow with NA values.
//'
//' The R6 object has advantages over the console.  All inputs are optional and are
//'  coordinated across functions.  Enter an input once and calculate several outputs.
//'  The R6 object is reactive.  In other words, it stores the inputs and outputs and
//'  maps inputs to outputs.  If an input changes, dependent outputs are nullified and
//'  will be recalculated, as requested, but nothing is calculated twice.  The console
//'  stores outputs in the global environment, but there is no map to inputs and they
//'  can be stale.  Another advantage of the  R6 object is predefined plots with Plotly.
//'  The same simulation can plotted different ways without recalculation.
//'
//' The parallel code has more overhead and is slower on small problems. Here
//'  are microbenchmark median times for simulating a small number of Forward Paths
//'  by calling the Rcpp and RcppParallel functions from the console:
//'
//'     Unit: milliseconds  Console       Console
//'                  paths     Rcpp   RcppParallel
//'     ------------------------------------------
//'                    100  0.06970        0.07600
//'                  1,000  0.55425        0.23540
//'
//' For 100 paths, the Rcpp sequential function takes less time, but for 1,000 paths
//'  it takes over twice as long.  Most simulations will have more than 1,000 paths.
//'  So users get no choice.  By default, RcppParallel functions are compiled if
//'  RcppParallel is installed. Otherwise compilation falls back to Rcpp.
//'  
//' Monte Carlo simulations are memory and CPU intensive.  Here are microbenchmark
//'  median times by number of threads for generating standard normal variables
//'  and simulating 1,000,000 Forward Paths over 100 time intervals using the
//'  stochastic integral equation with parameter skip=1:
//'
//'     Unit: milliseconds
//'                threads    stdnorm  ForwardPath      total
//'     -----------------------------------------------------
//'                      1  5461.3963     667.0016  6129.3979
//'                      2  2843.4390     438.4120  3281.8510
//'                      3  2877.2132     437.0354  3314.2486
//'                      4  1525.4355     326.4385  1851.8740
//'                      5  1281.7121     301.5675  1583.2796
//'                      6  1280.7483     301.8346  1582.5829
//'                      7  1121.1044     293.7598  1414.8642
//'                      8  1038.6441     306.8126  1345.4567
//'                      9  1070.7125     259.2108  1329.9233 
//'                     10   847.3448     283.3598  1130.7046
//'                     11   850.3528     282.3086  1132.6614
//'                     12   805.8151     284.3971  1090.2122
//'                      
//' These times are longer than previous times.  My computer seems tired today.
//'  But the changes in times by thread number are instructive.  Generating the
//'  standard normal variables would benefit from more threads.  Simulating the
//'  Forward Paths only needs six or eight threads.  Once you have generated the
//'  random variables, you could get by with fewer threads using the RcppParallel
//'  commands:
//'  
//'      defaultNumThreads()
//'      setThreadOptions(numThreads=8)
//'
//' Potentially, the Rcpp functions could be imported into other packages.
//'
//' @name MonteCarlo_Rcpp

// Helpers (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

double OUPVisitingTimeProbabilityInf(double x, double k, double rho, double mu, double sigma)
{
  double pinf;
  if(k == R_PosInf || k == R_NegInf) { pinf = 0; }
  else
  {
    if(k == mu) { pinf = 0.5; }
    else if(pow(sigma,2) < 0.0000000001)
    {
      if(x == k) { pinf = 1.0; }
      else if(x > k)
      {
        if(k > mu) { pinf = 1.0; }
        else { pinf = 0.0; }
      }
      else
      {
        if(k < mu) { pinf = 1.0; }
        else { pinf = 0.0; }
      }
    }
    else
    {
      double v2 = rho*pow(((k-mu)/sigma),2);
      if(x == k) { pinf = (1.77245385090552+GammaSmallOneHalf(v2))/(2*1.77245385090552); }
      else if(x > k)
      {
        if(k > mu) { pinf = (1.77245385090552+GammaSmallOneHalf(v2))/(2*1.77245385090552); }
        else { pinf = GammaBigOneHalf(v2)/(2*1.77245385090552); }
      }
      else
      {
        if(k < mu) { pinf = (1.77245385090552+GammaSmallOneHalf(v2))/(2*1.77245385090552); }
        else { pinf = GammaBigOneHalf(v2)/(2*1.77245385090552); }
      }
    }
  }
  return pinf;
}

// Exports (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

#ifdef USE_PARALLEL
struct ROMCPMM : public Worker
{
  const RVector<double> matPaths;
  double min;
  double max;

  ROMCPMM(const NumericVector& matPaths)
    : matPaths(matPaths), min(std::numeric_limits<double>::infinity()), max(-std::numeric_limits<double>::infinity()) {}
  ROMCPMM(const ROMCPMM& banana, Split)
    : matPaths(banana.matPaths), min(std::numeric_limits<double>::infinity()), max(-std::numeric_limits<double>::infinity()) {}

  void operator()(std::size_t begin, std::size_t end) {
    double localMin = min;
    double localMax = max;
    for(std::size_t i = begin; i < end; i++)
    {
      if(matPaths[i] < localMin) { localMin = matPaths[i]; }
      if(matPaths[i] > localMax) { localMax = matPaths[i]; }
    }
    min = localMin;
    max = localMax;
  }
  void join(const ROMCPMM& rhs) {
    if(rhs.min < min) { min = rhs.min; }
    if(rhs.max > max) { max = rhs.max; }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCMinMax(matPaths)
//' @param  matPaths matrix of paths
//' @return minmax(2) <- RcppOUPMCMinMax()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPMCMinMax(NumericVector matPaths)
{
  std::size_t paths = matPaths.size();
  NumericVector minmax(2);
#ifdef USE_PARALLEL
  ROMCPMM worker(matPaths);
  parallelReduce(0, paths, worker);
  minmax[0] = worker.min;
  minmax[1] = worker.max;
#else
  minmax[0] = std::numeric_limits<double>::infinity();
  minmax[1] = -std::numeric_limits<double>::infinity();
  for(std::size_t i = 0; i < paths; i++)
  {
    if(matPaths[i] < minmax[0]) { minmax[0] = matPaths[i]; }
    if(matPaths[i] > minmax[1]) { minmax[1] = matPaths[i]; }
  }
#endif
  return minmax;
}

#if defined(USE_PARALLEL) && defined(USE_SITMO)
struct ROMCPSN : public Worker
{
  RVector<double> stdnorm;
  uint64_t seed;

  ROMCPSN(NumericVector& stdnorm, uint64_t seed)
    : stdnorm(stdnorm), seed(seed) {}

  void operator()(std::size_t begin, std::size_t end) {
    sitmo::prng_engine engine(seed+begin);
    std::uniform_real_distribution<double> unif(0.0,1.0);
    for(std::size_t j = begin; j < end; j+=2)
    {
      double u1 = unif(engine);
      double u2 = unif(engine);
      if(u1 <= 0.0) u1 = std::numeric_limits<double>::min();
      double r = std::sqrt(-2.0*std::log(u1));
      double theta = 2*M_PI*u2;
      stdnorm[j] = r*std::cos(theta);
      if(j+1 < end) { stdnorm[j+1] = r*std::sin(theta); }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCStandardNormal(m,skip,paths,seed)
//' @param  m    number of rows for states over time
//' @param  skip subdivide time interval but report every ds or dt 0<skip<20
//' @param  paths number of columns for paths
//' @param  seed seed for reproducibility
//' @return stdnorm((m-1)*skip,paths) <- RcppOUPMCStandardNormal()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCStandardNormal(int64_t m, int64_t skip, int64_t paths, uint64_t seed)
{
#if defined(USE_PARALLEL) && defined(USE_SITMO)
  if((m-1) % 2 == 0)
  {
    NumericVector stdnorm((m-1)*skip*paths);
    ROMCPSN worker(stdnorm, seed);
    parallelFor(0, (m-1)*skip*paths, worker, 2);
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
  else
  {
    NumericVector stdnorm(m*skip*paths);
    ROMCPSN worker(stdnorm, seed);
    parallelFor(0, m*skip*paths, worker, 2);
    stdnorm.attr("dim") = Dimension(m*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#else
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);
  if((m-1) % 2 == 0)
  {
    NumericVector stdnorm = Rcpp::rnorm((m-1)*skip*paths);
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
  else
  {
    NumericVector stdnorm = Rcpp::rnorm(m*skip*paths);
    stdnorm.attr("dim") = Dimension(m*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#endif
}

#ifdef USE_PARALLEL
struct ROMCPFwRK : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> forward;
  double x;
  std::size_t m;
  std::size_t skip;
  double dtau;
  double rho;
  double mu;
  double H;

  ROMCPFwRK(const NumericMatrix& stdnorm, NumericMatrix& forward, double x, std::size_t m, std::size_t skip, double dtau, double rho, double mu, double H)
    : stdnorm(stdnorm), forward(forward), x(x), m(m), skip(skip), dtau(dtau), rho(rho), mu(mu), H(H) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double y = x;
      forward(0,j) = y;
      for(std::size_t i = 1; i < m; i++)
      {
        for(std::size_t p = 0; p < skip; p++)
        {
          std::size_t q = (i-1)*skip+p;
          double Heps = H*stdnorm(q,j);
          double G0 = -rho*(y-mu);
          double x1 = y+0.5*G0*dtau+0.5*Heps;
          double G1 = -rho*(x1-mu);
          double x2 = y+0.5*G1*dtau+0.5*Heps;
          double G2 = -rho*(x2-mu);
          double x3 = y+G2*dtau+Heps;
          double G3 = -rho*(x3-mu);
          y = y+(G0+2*G1+2*G2+G3)*dtau/6+Heps;
        }
        forward(i,j) = y;
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  x       initial state or vector of backward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<50
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return forward(m,paths) <- RcppOUPMCForwardPathRungeKutta()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCForwardPathRungeKutta(NumericMatrix stdnorm, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix forward(m,paths);
  double dtau = dt/skip;
  double H = sigma*sqrt(dtau);
#ifdef USE_PARALLEL
  ROMCPFwRK worker(stdnorm, forward, x, m, skip, dtau, rho, mu, H);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double y = x;
    forward(0,j) = y;
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t p = 0; p < skip; p++)
      {
        std::size_t q = (i-1)*skip+p;
        double Heps = H*stdnorm(q,j);
        double G0 = -rho*(y-mu);
        double x1 = y+0.5*G0*dtau+0.5*Heps;
        double G1 = -rho*(x1-mu);
        double x2 = y+0.5*G1*dtau+0.5*Heps;
        double G2 = -rho*(x2-mu);
        double x3 = y+G2*dtau+Heps;
        double G3 = -rho*(x3-mu);
        y = y+(G0+2*G1+2*G2+G3)*dtau/6+Heps;
      }
      forward(i,j) = y;
    }
  }
#endif
  return forward;
}

#ifdef USE_PARALLEL
struct ROMCPBkRK : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> backward;
  double y;
  std::size_t m;
  std::size_t skip;
  double dtau;
  double rho;
  double mu;
  double H;

  ROMCPBkRK(const NumericMatrix& stdnorm, NumericMatrix& backward, double y, std::size_t m, std::size_t skip, double dtau, double rho, double mu, double H)
    : stdnorm(stdnorm), backward(backward), y(y), m(m), skip(skip), dtau(dtau), rho(rho), mu(mu), H(H) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double x = y;
      std::size_t i = 0;
      backward(i,j) = x;
      for(i = 1; i < m; i++)
      {
        for(std::size_t p = 0; p < skip; p++)
        {
          std::size_t q = (i-1)*skip+p;
          double Heps = H*stdnorm(q,j);
          double G0 = rho*(x-mu);
          double y1 = x+0.5*G0*dtau-0.5*Heps;
          double G1 = rho*(y1-mu);
          double y2 = x+0.5*G1*dtau-0.5*Heps;
          double G2 = rho*(y2-mu);
          double y3 = x+G2*dtau-Heps;
          double G3 = rho*(y3-mu);
          x = x+(G0+2*G1+2*G2+G3)*dtau/6-Heps;
        }
        backward(i,j) = x;
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBackwardPathRungeKutta(stdnorm,y,m,skip,ds,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  y       terminal state or vector of forward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<50
//' @param  ds      time interval for terminal value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return backward(m,paths) <- RcppOUPMCBackwardPathRungeKutta()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBackwardPathRungeKutta(NumericMatrix stdnorm, double y, std::size_t m, std::size_t skip, double ds, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix backward(m,paths);
  double dtau = ds/skip;
  double H = sigma*sqrt(dtau);
#ifdef USE_PARALLEL
  ROMCPBkRK worker(stdnorm, backward, y, m, skip, dtau, rho, mu, H);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double x = y;
    backward(0,j) = x;
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t p = 0; p < skip; p++)
      {
        std::size_t q = (i-1)*skip+p;
        double Heps = H*stdnorm(q,j);
        double G0 = rho*(x-mu);
        double y1 = x+0.5*G0*dtau-0.5*Heps;
        double G1 = rho*(y1-mu);
        double y2 = x+0.5*G1*dtau-0.5*Heps;
        double G2 = rho*(y2-mu);
        double y3 = x+G2*dtau-Heps;
        double G3 = rho*(y3-mu);
        x = x+(G0+2*G1+2*G2+G3)*dtau/6-Heps;
      }
      backward(i,j) = x;
    }
  }
#endif
  return backward;
}

#ifdef USE_PARALLEL
struct ROMCPBdRK : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> bndfpt;
  double k;
  double x;
  std::size_t m;
  std::size_t skip;
  double dtau;
  double rho;
  double mu;
  double H;

  ROMCPBdRK(const NumericMatrix& stdnorm, NumericMatrix& bndfpt, double k, double x, std::size_t m, std::size_t skip, double dtau, double rho, double mu, double H)
    : stdnorm(stdnorm), bndfpt(bndfpt), k(k), x(x), m(m), skip(skip), dtau(dtau), rho(rho), mu(mu), H(H) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double y = x;
      bndfpt(0,j) = y;
      bool hit = false;
      std::size_t i = 1;
      while(i < m && hit == false)
      {
        std::size_t p = 0;
        while(p < skip && hit == false)
        {
          std::size_t q = (i-1)*skip+p;
          double G0 = -rho*(y-mu);
          double x1 = y+0.5*G0*dtau+0.5*H*stdnorm(q,j);
          double G1 = -rho*(x1-mu);
          double x2 = y+0.5*G1*dtau+0.5*H*stdnorm(q,j);
          double G2 = -rho*(x2-mu);
          double x3 = y+G2*dtau+H*stdnorm(q,j);
          double G3 = -rho*(x3-mu);
          double newy = y+(G0+2*G1+2*G2+G3)*dtau/6+H*stdnorm(q,j);
          if(((x >= k) && (k >= newy)) || ((x <= k) && (k <= newy)))
          {
            hit = true;
            bndfpt(m,j) = (q+(k-y)/(newy-y))*dtau;
          }
          y = newy;
          p += 1;
        }
        bndfpt(i,j) = y;
        i += 1;
      }
      while(i < m)
      {
        bndfpt(i,j) = NA_REAL;
        i += 1;
      }
      if(!hit) { bndfpt(m,j) = NA_REAL; }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  k       threshold -inf<k<inf
//' @param  x       initial state or vector of backward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<20
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return bndfpt(m+1,paths) <- RcppOUPMCBoundedPathRungeKutta()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBoundedPathRungeKutta(NumericMatrix stdnorm, double k, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix bndfpt(m+1,paths);
  double dtau = dt/skip;
  double H = sigma*sqrt(dtau);
#ifdef USE_PARALLEL
  ROMCPBdRK worker(stdnorm, bndfpt, k, x, m, skip, dtau, rho, mu, H);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double y = x;
    bndfpt(0,j) = y;
    bool hit = false;
    std::size_t i = 1;
    while(i < m && hit == false)
    {
      std::size_t p = 0;
      while(p < skip && hit == false)
      {
        std::size_t q = (i-1)*skip+p;
        double G0 = -rho*(y-mu);
        double x1 = y+0.5*G0*dtau+0.5*H*stdnorm(q,j);
        double G1 = -rho*(x1-mu);
        double x2 = y+0.5*G1*dtau+0.5*H*stdnorm(q,j);
        double G2 = -rho*(x2-mu);
        double x3 = y+G2*dtau+H*stdnorm(q,j);
        double G3 = -rho*(x3-mu);
        double newy = y+(G0+2*G1+2*G2+G3)*dtau/6+H*stdnorm(q,j);
        if(((x >= k) && (k >= newy)) || ((x <= k) && (k <= newy)))
        {
          hit = true;
          bndfpt(m,j) = (q+(k-y)/(newy-y))*dtau;
        }
        y = newy;
        p += 1;
      }
      bndfpt(i,j) = y;
      i += 1;
    }
    while(i < m)
    {
      bndfpt(i,j) = NA_REAL;
      i += 1;
    }
    if(!hit) { bndfpt(m,j) = NA_REAL; }
  }
#endif
  return bndfpt;
}

#ifdef USE_PARALLEL
struct ROMCPFwIE : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> forward;
  double x;
  std::size_t m;
  std::size_t skip;
  double rho;
  double mu;
  double H;
  double exprhodt;

  ROMCPFwIE(const NumericMatrix& stdnorm, NumericMatrix& forward, double x, std::size_t m, std::size_t skip, double rho, double mu, double H, double exprhodt)
    : stdnorm(stdnorm), forward(forward), x(x), m(m), skip(skip), rho(rho), mu(mu), H(H), exprhodt(exprhodt) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double y = x;
      forward(0,j) = y;
      for(std::size_t i = 1; i < m; i++)
      {
        for(std::size_t p = 0; p < skip; p++)
        {
          std::size_t q = (i-1)*skip+p;
          y = mu+(y-mu)*exprhodt+H*stdnorm(q,j);
        }
        forward(i,j) = y;
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  x       initial state or vector of backward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<20
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return forward(m,paths) <- RcppOUPMCForwardPathIntegralEquation()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCForwardPathIntegralEquation(NumericMatrix stdnorm, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix forward(m,paths);
  double dtau = dt/skip;
  double H = sigma*sqrt(dtau);
  if(rho > 0) { H = sqrt(sigma*sigma/(2*rho)*(1-exp(-2*rho*dtau))); }
  double exprhodt = exp(-rho*dtau);
#ifdef USE_PARALLEL
  ROMCPFwIE worker(stdnorm, forward, x, m, skip, rho, mu, H, exprhodt);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double y = x;
    forward(0,j) = y;
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t p = 0; p < skip; p++)
      {
        std::size_t q = (i-1)*skip+p;
        y = mu+(y-mu)*exprhodt+H*stdnorm(q,j);
      }
      forward(i,j) = y;
    }
  }
#endif
  return forward;
}

#ifdef USE_PARALLEL
struct ROMCPBkIE : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> backward;
  double y;
  std::size_t m;
  std::size_t skip;
  double rho;
  double mu;
  double H;
  double exprhods;

  ROMCPBkIE(const NumericMatrix& stdnorm, NumericMatrix& backward, double y, std::size_t m, std::size_t skip, double rho, double mu, double H, double exprhods)
    : stdnorm(stdnorm), backward(backward), y(y), m(m), skip(skip), rho(rho), mu(mu), H(H), exprhods(exprhods) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double x = y;
      backward(0,j) = x;
      for(std::size_t i = 1; i < m; i++)
      {
        for(std::size_t p = 0; p < skip; p++)
        {
          std::size_t q = (i-1)*skip+p;
          x = mu+(x-mu)*exprhods-H*stdnorm(q,j);
        }
        backward(i,j) = x;
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBackwardPathIntegralEquation(stdnorm,y,m,skip,ds,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  y       terminal state or vector of forward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<20
//' @param  ds      time interval for terminal value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return backward(m,paths) <- RcppOUPMCBackwardPathIntegralEquation()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBackwardPathIntegralEquation(NumericMatrix stdnorm, double y, std::size_t m, std::size_t skip, double ds, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix backward(m,paths);
  double dtau = ds/skip;
  double H = sigma*sqrt(dtau);
  if(rho > 0) { H = sqrt(sigma*sigma/(2*rho)*(exp(2*rho*dtau)-1)); }
  double exprhods = exp(rho*dtau);
#ifdef USE_PARALLEL
  ROMCPBkIE worker(stdnorm, backward, y, m, skip, rho, mu, H, exprhods);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double x = y;
    backward(0,j) = x;
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t p = 0; p < skip; p++)
      {
        std::size_t q = (i-1)*skip+p;
        x = mu+(x-mu)*exprhods-H*stdnorm(q,j);
      }
      backward(i,j) = x;
    }
  }
#endif
  return backward;
}

#ifdef USE_PARALLEL
struct ROMCPBdIE : public Worker
{
  const RMatrix<double> stdnorm;
  RMatrix<double> bndfpt;
  double k;
  double x;
  std::size_t m;
  std::size_t skip;
  double dtau;
  double rho;
  double mu;
  double H;
  double exprhodt;

  ROMCPBdIE(const NumericMatrix& stdnorm, NumericMatrix& bndfpt, double k, double x, std::size_t m, std::size_t skip, double dtau, double rho, double mu, double H, double exprhodt)
    : stdnorm(stdnorm), bndfpt(bndfpt), k(k), x(x), m(m), skip(skip), dtau(dtau), rho(rho), mu(mu), H(H), exprhodt(exprhodt) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++)
    {
      double y = x;
      bndfpt(0,j) = y;
      bool hit = false;
      std::size_t i = 1;
      while(i < m && hit == false)
      {
        std::size_t p = 0;
        while(p < skip && hit == false)
        {
          std::size_t q = (i-1)*skip+p;
          double newy = mu+(y-mu)*exprhodt+H*stdnorm(q,j);
          if(((x >= k) && (k >= newy)) || ((x <= k) && (k <= newy)))
          {
            hit = true;
            bndfpt(m,j) = (q+(k-y)/(newy-y))*dtau;
          }
          y = newy;
          p += 1;
        }
        bndfpt(i,j) = y;
        i += 1;
      }
      while(i < m)
      {
        bndfpt(i,j) = NA_REAL;
        i += 1;
      }
      if(!hit) { bndfpt(m,j) = NA_REAL; }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma)
//' @param  stdnorm matrix of standard normal shocks
//' @param  x       initial state or vector of backward states
//' @param  k       threshold -inf<k<inf
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<20
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @return bndfpt(m+1,paths) <- RcppOUPMCBoundedPathIntegralEquation()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBoundedPathIntegralEquation(NumericMatrix stdnorm, double k, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma)
{
  std::size_t paths = stdnorm.ncol();
  NumericMatrix bndfpt(m+1,paths);
  double dtau = dt/skip;
  double H = sigma*sqrt(dtau);
  if(rho > 0) { H = sqrt(sigma*sigma/(2*rho)*(1-exp(-2*rho*dtau))); }
  double exprhodt = exp(-rho*dtau);
#ifdef USE_PARALLEL
  ROMCPBdIE worker(stdnorm, bndfpt, k, x, m, skip, dtau, rho, mu, H, exprhodt);
  parallelFor(0, paths, worker);
#else
  for(std::size_t j = 0; j < paths; j++)
  {
    double y = x;
    bndfpt(0,j) = y;
    bool hit = false;
    std::size_t i = 1;
    while(i < m && hit == false)
    {
      std::size_t p = 0;
      while(p < skip && hit == false)
      {
        std::size_t q = (i-1)*skip+p;
        double newy = mu+(y-mu)*exprhodt+H*stdnorm(q,j);
        if(((x >= k) && (k >= newy)) || ((x <= k) && (k <= newy)))
        {
          hit = true;
          bndfpt(m,j) = (q+(k-y)/(newy-y))*dtau;
        }
        y = newy;
        p += 1;
      }
      bndfpt(i,j) = y;
      i += 1;
    }
    while(i < m)
    {
      bndfpt(i,j) = NA_REAL;
      i += 1;
    }
    if(!hit) { bndfpt(m,j) = NA_REAL; }
  }
#endif
  return bndfpt;
}

#ifdef USE_PARALLEL
struct ROMCPmvpPPP : public Worker
{
  const RMatrix<double> forward;
  RMatrix<double> mvdpd;
  RMatrix<double> dens;
  RMatrix<double> prob;
  RMatrix<double> doub;
  std::size_t paths;
  std::size_t n;
  std::size_t n0;
  std::size_t n1;
  std::size_t n2;
  std::size_t nn;
  std::size_t offset;
  double width;
  double ymin;
  double psi;

  ROMCPmvpPPP(const NumericMatrix& forward, NumericMatrix& mvdpd, NumericMatrix& dens, NumericMatrix& prob, NumericMatrix& doub, std::size_t paths, std::size_t n, std::size_t n0, std::size_t n1, std::size_t n2, std::size_t nn, std::size_t offset, double width, double ymin, double psi)
    : forward(forward), mvdpd(mvdpd), dens(dens), prob(prob), doub(doub), paths(paths), n(n), n0(n0), n1(n1), n2(n2), nn(nn), offset(offset), width(width), ymin(ymin), psi(psi) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++)
    {
      for(std::size_t j = 0; j < paths; j++)
      {
        mvdpd(i,0) += forward(i,j);
        mvdpd(i,1) += pow(forward(i,j),2);
        std::size_t bin = static_cast<std::size_t>((forward(i,j)-ymin)/width);
        dens(i,bin) += 1;
      }
      mvdpd(i,0) /= paths;
      mvdpd(i,1) = mvdpd(i,1)/paths-pow(mvdpd(i,0),2);
      if(psi > 0)
      {
        dens(i,nn-1) = dens(i,nn-1)/paths;
        prob(i,nn-1) = dens(i,nn-1);
        doub(i,nn-1) = 0;
        dens(i,nn-1) = dens(i,nn-1)/width;
        for(std::size_t j = nn-1; j > 0+offset; j--)
        {
          dens(i,j-1) = dens(i,j-1)/paths;
          prob(i,j-1) = prob(i,j)+dens(i,j-1);
          doub(i,j-1) = doub(i,j)+prob(i,j)*width;
          dens(i,j-1) = dens(i,j-1)/width;
        }
      }
      else
      {
        dens(i,0) = dens(i,0)/paths;
        prob(i,0) = dens(i,0);
        doub(i,0) = 0;
        dens(i,0) = dens(i,0)/width;
        for(std::size_t j = 1; j < n+offset; j++)
        {
          dens(i,j) = dens(i,j)/paths;
          prob(i,j) = prob(i,j-1)+dens(i,j);
          doub(i,j) = doub(i,j-1)+prob(i,j-1)*width;
          dens(i,j) = dens(i,j)/width;
        }
      }
      for(std::size_t j = 0; j < n; j++)
      {
        mvdpd(i,j+n0) = dens(i,j+offset);
        mvdpd(i,j+n1) = prob(i,j+offset);
        mvdpd(i,j+n2) = doub(i,j+offset);
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCForwardCountY(forward,y,psi)
//' @param  forward matrix of forward paths
//' @param  y       terminal state or vector of forward states
//' @param  psi     <=0 for integral -inf to y, >0 for integral y to inf
//' @return mvdpd(m,3*n+2) <- RcppOUPMCForwardCountY()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCForwardCountY(NumericMatrix forward, NumericVector y, double psi)
{
  std::size_t paths = forward.ncol();
  std::size_t m = forward.nrow();
  std::size_t n = y.size();
  std::size_t n0 = 2;
  std::size_t n1 = n+2;
  std::size_t n2 = n*2+2;
  double x = forward(0,0);
  double width;
  double ymin;
  double ymax;
  std::size_t nn;
  std::size_t offset;
  NumericVector minmax(2);
  if(n > 1) { width = (y[n-1]-y[0])/(n-1); }
  else { width = 1; }
  minmax = RcppOUPMCMinMax(forward);
  if(minmax[0] > y[0]-0.5*width) { ymin = y[0]-0.5*width; }
  else { ymin = y[0]-0.5*width-static_cast<int>((y[0]-minmax[0])/width+0.5)*width; }
  if(minmax[1] < y[0]-0.5*width+n*width) { ymax = y[0]-0.5*width+n*width; }
  else { ymax = y[0]-0.5*width+n*width+static_cast<int>((minmax[1]-y[0])/width-n+1.5)*width; }
  nn = static_cast<int>((ymax-ymin)/width+0.5);
  offset = static_cast<int>((y[0]-ymin)/width);
  NumericMatrix dens(m,nn);
  NumericMatrix prob(m,nn);
  NumericMatrix doub(m,nn);
  NumericMatrix mvdpd(m,3*n+2); // means,variances,densities,probabilities,doubleintegrals
  mvdpd(0,0) = x;
  mvdpd(0,1) = 0;
  if(psi > 0)
  {
    for(std::size_t j = n; j > 0; j--)
    {
      if(y[j-1] > x+0.5*width || y[j-1] <= x-0.5*width) { mvdpd(0,j-1+n0) = 0; }
      else { mvdpd(0,j-1+n0) = 1/width; }
      if(y[j-1] > x)
      {
        mvdpd(0,j-1+n1) = 0;
        mvdpd(0,j-1+n2) = 0;
      }
      else if(y[j-1] < x)
      {
        mvdpd(0,j-1+n1) = 1;
        mvdpd(0,j-1+n2) = x-y[j-1];
      }
      else
      {
        mvdpd(0,j-1+n1) = 0.5;
        mvdpd(0,j-1+n2) = 0;
      }
    }
  }
  else
  {
    for(std::size_t j = 0; j < n; j++)
    {
      if(y[j] < x-0.5*width || y[j] >= x+0.5*width) { mvdpd(0,j+n0) = 0; }
      else { mvdpd(0,j+n0) = 1/width; }
      if(y[j] < x)
      {
        mvdpd(0,j+n1) = 0;
        mvdpd(0,j+n2) = 0;
      }
      else if(y[j] > x)
      {
        mvdpd(0,j+n1) = 1;
        mvdpd(0,j+n2) = y[j]-x;
      }
      else
      {
        mvdpd(0,j+n1) = 0.5;
        mvdpd(0,j+n2) = 0;
      }
    }
  }
#ifdef USE_PARALLEL
  ROMCPmvpPPP worker(forward, mvdpd, dens, prob, doub, paths, n, n0, n1, n2, nn, offset, width, ymin, psi);
  parallelFor(1, m, worker);
#else
  for(std::size_t i = 1; i < m; i++)
  {
    for(std::size_t j = 0; j < paths; j++)
    {
      mvdpd(i,0) += forward(i,j);
      mvdpd(i,1) += pow(forward(i,j),2);
      std::size_t bin = static_cast<int>((forward(i,j)-ymin)/width);
      dens(i,bin) += 1;
    }
    mvdpd(i,0) /= paths;
    mvdpd(i,1) = mvdpd(i,1)/paths-pow(mvdpd(i,0),2);
    if(psi > 0)
    {
      dens(i,nn-1) = dens(i,nn-1)/paths;
      prob(i,nn-1) = dens(i,nn-1);
      doub(i,nn-1) = 0;
      dens(i,nn-1) = dens(i,nn-1)/width;
      for(std::size_t j = nn-1; j > 0+offset; j--)
      {
        dens(i,j-1) = dens(i,j-1)/paths;
        prob(i,j-1) = prob(i,j)+dens(i,j-1);
        doub(i,j-1) = doub(i,j)+prob(i,j)*width;
        dens(i,j-1) = dens(i,j-1)/width;
      }
    }
    else
    {
      dens(i,0) = dens(i,0)/paths;
      prob(i,0) = dens(i,0);
      doub(i,0) = 0;
      dens(i,0) = dens(i,0)/width;
      for(std::size_t j = 1; j < n+offset; j++)
      {
        dens(i,j) = dens(i,j)/paths;
        prob(i,j) = prob(i,j-1)+dens(i,j);
        doub(i,j) = doub(i,j-1)+prob(i,j-1)*width;
        dens(i,j) = dens(i,j)/width;
      }
    }
    for(std::size_t j = 0; j < n; j++)
    {
      mvdpd(i,j+n0) = dens(i,j+offset);
      mvdpd(i,j+n1) = prob(i,j+offset);
      mvdpd(i,j+n2) = doub(i,j+offset);
    }
  }
#endif
  return mvdpd;
}

#ifdef USE_PARALLEL
struct ROMCPoOOO : public Worker
{
  const RMatrix<double> backward;
  RMatrix<double> dpo;
  RMatrix<double> dens;
  RMatrix<double> prob;
  RMatrix<double> optn;
  std::size_t paths;
  std::size_t n;
  std::size_t n0;
  std::size_t n1;
  std::size_t n2;
  std::size_t nn;
  std::size_t offset;
  double width;
  double xmin;
  double phi;
  double rho;
  double r;
  double ds;

  ROMCPoOOO(const NumericMatrix& backward, NumericMatrix& dpo, NumericMatrix& dens, NumericMatrix& prob, NumericMatrix& optn, std::size_t paths, std::size_t n, std::size_t n0, std::size_t n1, std::size_t n2, std::size_t nn, std::size_t offset, double width, double xmin, double phi, double rho, double r, double ds)
    : backward(backward), dpo(dpo), dens(dens), prob(prob), optn(optn), paths(paths), n(n), n0(n0), n1(n1), n2(n2), nn(nn), offset(offset), width(width), xmin(xmin), phi(phi), rho(rho), r(r), ds(ds) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++)
    {
      for(std::size_t j = 0; j < paths; j++)
      {
          std::size_t bin = static_cast<std::size_t>((backward(i,j)-xmin)/width);
          dens(i,bin) += 1;
      }
      if(phi > 0)
      {
        dens(i,0) = dens(i,0)/paths;
        prob(i,0) = dens(i,0);
        optn(i,0) = 0;
        dens(i,0) = dens(i,0)/width;
        for(std::size_t j = 1; j < n+offset; j++)
        {
          dens(i,j) = dens(i,j)/paths;
          prob(i,j) = prob(i,j-1)+dens(i,j);
          optn(i,j) = optn(i,j-1)+prob(i,j-1)*width;
          dens(i,j) = dens(i,j)/width;
        }
      }
      else
      {
        dens(i,nn-1) = dens(i,nn-1)/paths;
        prob(i,nn-1) = dens(i,nn-1);
        optn(i,nn-1) = 0;
        dens(i,nn-1) = dens(i,nn-1)/width;
        for(std::size_t j = nn-1; j > 0+offset; j--)
        {
          dens(i,j-1) = dens(i,j-1)/paths;
          prob(i,j-1) = prob(i,j)+dens(i,j-1);
          optn(i,j-1) = optn(i,j)+prob(i,j)*width;
          dens(i,j-1) = dens(i,j-1)/width;
        }
      }
      double exprhods = exp(rho*i*ds);
      double exprhords= exp(-(rho+r)*i*ds);
      for(std::size_t j = 0; j < n; j++)
      {
        dpo(i,j+n0) = exprhods*dens(i,j+offset);
        dpo(i,j+n1) = prob(i,j+offset);
        dpo(i,j+n2) = exprhords*optn(i,j+offset);
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBackwardCountX(backward,x,phi,rho,r,ds)
//' @param  backward matrix of backward paths
//' @param  x       initial state or vector of backward states
//' @param  phi     <=0 for integral -inf to x, >0 for integral x to inf
//' @param  rho      rate parameter 0<=rho<inf
//' @param  r        discount rate 0<r
//' @param  ds       time interval for terminal value problems
//' @return dpo(m,3*n) <- RcppOUPMCBackwardCountX()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBackwardCountX(NumericMatrix backward, NumericVector x, double phi, double rho, double r, double ds)
{
  std::size_t paths = backward.ncol();
  std::size_t m = backward.nrow();
  std::size_t n = x.size();
  std::size_t n0 = 0;
  std::size_t n1 = n;
  std::size_t n2 = n*2;
  double y = backward(0,0);
  double width;
  double xmin;
  double xmax;
  std::size_t nn;
  std::size_t offset;
  NumericVector minmax(2);
  if(n > 1) { width = (x[n-1]-x[0])/(n-1); }
  else { width = 1; }
  minmax = RcppOUPMCMinMax(backward);
  if(minmax[0] > x[0]-0.5*width) { xmin = x[0]-0.5*width; }
  else { xmin = x[0]-0.5*width-static_cast<int>((x[0]-minmax[0])/width+0.5)*width; }
  if(minmax[1] < x[0]-0.5*width+n*width) { xmax = x[0]-0.5*width+n*width; }
  else { xmax = x[0]-0.5*width+n*width+static_cast<int>((minmax[1]-x[0])/width-n+1.5)*width; }
  nn = static_cast<int>((xmax-xmin)/width+0.5);
  offset = static_cast<int>((x[0]-xmin)/width);
  NumericMatrix dens(m,nn);
  NumericMatrix prob(m,nn);
  NumericMatrix optn(m,nn);
  NumericMatrix dpo(m,3*n); // densities[,(0):(n-1)], probabilities[,(n):(2n-1)], doubleintegrals[,(2n):(3n-1)]
  if(phi > 0)
  {
    for(std::size_t j = 0; j < n; j++)
    {
      if(x[j] < y-0.5*width || x[j] >= y+0.5*width) { dpo(0,j+n0) = 0; }
      else { dpo(0,j+n0) = 1/width; }
      if(x[j] < y)
      {
        dpo(0,j+n1) = 0;
        dpo(0,j+n2) = 0;
      }
      else if(x[j] > y)
      {
        dpo(0,j+n1) = 1;
        dpo(0,j+n2) = x[j]-y;
      }
      else
      {
        dpo(0,j+n1) = 0.5;
        dpo(0,j+n2) = 0;
      }
    }
  }
  else
  {
    for(std::size_t j = n; j > 0; j--)
    {
      if(x[j-1] > y+0.5*width || x[j-1] <= y-0.5*width) { dpo(0,j-1+n0) = 0; }
      else { dpo(0,j-1+n0) = 1/width; }
      if(x[j-1] > y)
      {
        dpo(0,j-1+n1) = 0;
        dpo(0,j-1+n2) = 0;
      }
      else if(x[j-1] < y)
      {
        dpo(0,j-1+n1) = 1;
        dpo(0,j-1+n2) = y-x[j-1];
      }
      else
      {
        dpo(0,j-1+n1) = 0.5;
        dpo(0,j-1+n2) = 0;
      }
    }
  }
#ifdef USE_PARALLEL
  ROMCPoOOO worker(backward, dpo, dens, prob, optn, paths, n, n0, n1, n2, nn, offset, width, xmin, phi, rho, r, ds);
  parallelFor(1, m, worker);
#else
  for(std::size_t i = 1; i < m; i++)
  {
    for(std::size_t j = 0; j < paths; j++)
    {
      std::size_t bin = static_cast<int>((backward(i,j)-xmin)/width);
      dens(i,bin) += 1;
    }
    if(phi > 0)
    {
      dens(i,0) = dens(i,0)/paths;
      prob(i,0) = dens(i,0);
      optn(i,0) = 0;
      dens(i,0) = dens(i,0)/width;
      for(std::size_t j = 1; j < n+offset; j++)
      {
        dens(i,j) = dens(i,j)/paths;
        prob(i,j) = prob(i,j-1)+dens(i,j);
        optn(i,j) = optn(i,j-1)+prob(i,j-1)*width;
        dens(i,j) = dens(i,j)/width;
      }
    }
    else
    {
      dens(i,nn-1) = dens(i,nn-1)/paths;
      prob(i,nn-1) = dens(i,nn-1);
      optn(i,nn-1) = 0;
      dens(i,nn-1) = dens(i,nn-1)/width;
      for(std::size_t j = nn-1; j > 0+offset; j--)
      {
        dens(i,j-1) = dens(i,j-1)/paths;
        prob(i,j-1) = prob(i,j)+dens(i,j-1);
        optn(i,j-1) = optn(i,j)+prob(i,j)*width;
        dens(i,j-1) = dens(i,j-1)/width;
      }
    }
    double exprhods = exp(rho*i*ds);
    double exprhords= exp(-(rho+r)*i*ds);
    for(std::size_t j = 0; j < n; j++)
    {
      dpo(i,j+n0) = exprhods*dens(i,j+offset);
      dpo(i,j+n1) = prob(i,j+offset);
      dpo(i,j+n2) = exprhords*optn(i,j+offset);
    }
  }
#endif
  return dpo;
}

#ifdef USE_PARALLEL
struct ROMCPpctpvPv : public Worker
{
  const RMatrix<double> forward;
  RMatrix<double> pctdp;
  std::size_t paths;
  double k;
  double x;
  double mu;

  ROMCPpctpvPv(const NumericMatrix& forward, NumericMatrix& pctdp, std::size_t paths, double k, double x, double mu)
    : forward(forward), pctdp(pctdp), paths(paths), k(k), x(x), mu(mu) {}

  void operator()(std::size_t begin, std::size_t end) {
    if((x > k) || (x == k && k >= mu))
    {
      for(std::size_t i = begin; i < end; i++)
      {
        for(std::size_t j = 0; j < paths; j++)
        {
          if(k >= forward(i,j)) { pctdp(i,4) += 1; }
        }
        pctdp(i,4) /= paths;
      }
    }
    else if((x < k) || (x == k && k < mu))
    {
      for(std::size_t i = begin; i < end; i++)
      {
        for(std::size_t j = 0; j < paths; j++)
        {
          if(k <= forward(i,j)) { pctdp(i,4) += 1; }
        }
        pctdp(i,4) /= paths;
      }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
//' @param  forward matrix of forward paths
//' @param  k       threshold -inf<k<inf
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @param  Ppct    probability for a percentile 0.01<pct<0.99
//' @return pctdp(m,5) <- RcppOUPMCForwardCountT()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCForwardCountT(NumericMatrix forward, double k, double dt, double rho, double mu, double sigma, double Ppct)
{
  std::size_t paths = forward.ncol();
  std::size_t m = forward.nrow();
  std::size_t m5 = m;
  if(m5 < 5) { m5 = 5; }
  NumericMatrix pctdp(m5,5); // tpercentile,pvpercentile,Pvpercentile,densities,probabilities
  double pctlow;
  double pcthigh;
  double pinf;
  double x = forward(0,0);
  if(x == k)
  {
    pctdp(0,4) = 0.5;
    pctdp(0,3) = 9;
  }
  else
  {
    pctdp(0,4) = 0;
    pctdp(0,3) = 0;
  }
#ifdef USE_PARALLEL
  ROMCPpctpvPv worker(forward, pctdp, paths, k, x, mu);
  parallelFor(1, m, worker);
#else
  if((x > k) || (x == k && k >= mu))
  {
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t j = 0; j < paths; j++)
      {
        if(k >= forward(i,j)) { pctdp(i,4) += 1; }
      }
      pctdp(i,4) /= paths;
    }
  }
  else if((x < k) || (x == k && k < mu))
  {
    for(std::size_t i = 1; i < m; i++)
    {
      for(std::size_t j = 0; j < paths; j++)
      {
        if(k <= forward(i,j)) { pctdp(i,4) += 1; }
      }
      pctdp(i,4) /= paths;
    }
  }
#endif
  for(std::size_t i = 1; i < m-1; i++) { pctdp(i,3) = (pctdp(i+1,4)-pctdp(i-1,4))/(2*dt); }
  if(m > 1)
  {
    pctdp(m-1,3) = (pctdp(m-1,4)-pctdp(m-2,4))/dt;
    pinf = OUPVisitingTimeProbabilityInf(x,k,rho,mu,sigma);
    for(std::size_t i = 1; i < m; i++)
    {
      if(pctdp(0,1) < pctdp(i,3))
      {
        pctdp(0,0) = i*dt;
        pctdp(0,1) = pctdp(i,3);
        pctdp(0,2) = pctdp(i,4);
      }
      pctdp(2,0) += pctdp(i,4);
    }
    pctdp(2,0) = m-1-(pctdp(2,0)-0.5*(pctdp(0,4)+pctdp(m-1,4)))/pinf;
    std::size_t i = static_cast<int>(pctdp(2,0));
    if(i < 0 || i > m-1)
    {
      pctdp(2,1) = NA_REAL;
      pctdp(2,2) = NA_REAL;
    }
    else if(i == m-1)
    {
      pctdp(2,1) = pctdp(m-1,3);
      pctdp(2,2) = pctdp(m-1,4);
    }
    else
    {
      pctdp(2,1) = pctdp(i,3)*(i+1-pctdp(2,0))+pctdp(i+1,3)*(pctdp(2,0)-i);
      pctdp(2,2) = pctdp(i,4)*(i+1-pctdp(2,0))+pctdp(i+1,4)*(pctdp(2,0)-i);
    }
    pctdp(2,0) *= dt;
    if(Ppct < 0.5)
    {
      pctlow = Ppct;
      pcthigh = 1-Ppct;
    }
    else
    {
      pctlow = 1-Ppct;
      pcthigh = Ppct;
    }
    i = 1;
    while(pctdp(i,4) < pctlow*pinf && i < m-1) { i += 1;}
    if(pctdp(i,4) >= pctlow*pinf)
    {
      pctdp(3,0) = (i-(pctdp(i,4)-pctlow*pinf)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(3,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(3,0)/dt);
      pctdp(3,2) = pctlow*pinf;
    }
    else
    {
      pctdp(3,0) = NA_REAL;
      pctdp(3,1) = NA_REAL;
      pctdp(3,2) = NA_REAL;
    }
    while(pctdp(i,4) < 0.5*pinf && i < m-1) { i += 1;}
    if(pctdp(i,4) >= 0.5*pinf)
    {
      pctdp(1,0) = (i-(pctdp(i,4)-0.5*pinf)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(1,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(1,0)/dt);
      pctdp(1,2) = 0.5*pinf;
    }
    else
    {
      pctdp(1,0) = NA_REAL;
      pctdp(1,1) = NA_REAL;
      pctdp(1,2) = NA_REAL;
    }
    while(pctdp(i,4) < pcthigh*pinf && i < m-1) { i += 1;}
    if(pctdp(i,4) >= pcthigh*pinf)
    {
      pctdp(4,0) = (i-(pctdp(i,4)-pcthigh*pinf)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(4,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(4,0)/dt);
      pctdp(4,2) = pcthigh*pinf;
    }
    else
    {
      pctdp(4,0) = NA_REAL;
      pctdp(4,1) = NA_REAL;
      pctdp(4,2) = NA_REAL;
    }
  }
  return pctdp;
}

#ifdef USE_PARALLEL
struct ROMCPpctdp : public Worker
{
  const RVector<double> fpt;
  RMatrix<double> pctdp;
  double dt;

  ROMCPpctdp(const NumericVector& fpt, NumericMatrix& pctdp, double dt)
    : fpt(fpt), pctdp(pctdp), dt(dt) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t cnt = 0;
    for(std::size_t j = begin; j < end; j++)
    {
      if(!Rcpp::traits::is_na<REALSXP>(fpt[j]))
      {
        std::size_t bin = static_cast<int>(fpt[j]/dt+0.5);
        pctdp(bin,3) += 1;
        pctdp(2,0) += fpt[j];
        cnt += 1;
      }
    }
    pctdp(2,0) /= cnt;
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)
//' @param  fpt     vector of first passage times
//' @param  m       number of rows for states over time
//' @param  dt      time interval for initial value problems
//' @param  Ppct    probability for a percentile 0.01<pct<0.99
//' @return pctdp(m,5) <- RcppOUPMCBoundedCountT()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBoundedCountT(NumericVector fpt, std::size_t m, double dt, double Ppct)
{
  std::size_t paths = fpt.size();
  std::size_t m5 = m;
  if(m5 < 5) { m5 = 5; }
  NumericMatrix pctdp(m5,5); // tpercentile,pfpercentile,Pfpercentile,densities,probabilities
  double pctlow;
  double pcthigh;
#ifdef USE_PARALLEL
  ROMCPpctdp worker(fpt, pctdp, dt);
  parallelFor(1, paths, worker);
#else
  std::size_t cnt = 0;
  for(std::size_t j = 0; j < paths; j++)
  {
    if(!Rcpp::traits::is_na<REALSXP>(fpt[j]))
    {
      std::size_t bin = static_cast<int>(fpt[j]/dt+0.5);
      pctdp(bin,3) += 1;
      pctdp(2,0) += fpt[j];
      cnt += 1;
    }
  }
  pctdp(2,0) /= cnt;
#endif
  pctdp(0,3) /= paths;
  pctdp(0,4) = pctdp(0,3);
  pctdp(0,3) /= dt;
  for(std::size_t i = 1; i < m; i++)
  {
    pctdp(i,3) /= paths;
    pctdp(i,4) = pctdp(i-1,4) + pctdp(i,3);
    pctdp(i,3) /= dt;
  }
  pctdp(m-1,4) += pctdp(m-1,3)*dt;
  pctdp(m-1,3) *= 2;
  if(m > 1)
  {
    for(std::size_t i = 1; i < m; i++)
    {
      if(pctdp(0,1) < pctdp(i,3))
      {
        pctdp(0,0) = i*dt;
        pctdp(0,1) = pctdp(i,3);
        pctdp(0,2) = pctdp(i,4);
      }
    }
    std::size_t i = static_cast<int>(pctdp(2,0)/dt);
    if(i < 0 || i > m-1)
      {
      pctdp(2,1) = NA_REAL;
      pctdp(2,2) = NA_REAL;
      }
    else if(i == m-1)
    {
      pctdp(2,1) = pctdp(m-1,3);
      pctdp(2,2) = pctdp(m-1,4);
    }
    else
    {
      pctdp(2,1) = pctdp(i,3)*(i+1-pctdp(2,0)/dt)+pctdp(i+1,3)*(pctdp(2,0)/dt-i);
      pctdp(2,2) = pctdp(i,4)*(i+1-pctdp(2,0)/dt)+pctdp(i+1,4)*(pctdp(2,0)/dt-i);
    }
    if(Ppct < 0.5)
    {
      pctlow = Ppct;
      pcthigh = 1-Ppct;
    }
    else
    {
      pctlow = 1-Ppct;
      pcthigh = Ppct;
    }
    i = 1;
    while(pctdp(i,4) < pctlow && i < m-1) { i += 1;}
    if(pctdp(i,4) >= pctlow)
    {
      pctdp(3,0) = (i-(pctdp(i,4)-pctlow)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(3,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(3,0)/dt);
      pctdp(3,2) = pctlow;
    }
    else
    {
      pctdp(3,0) = NA_REAL;
      pctdp(3,1) = NA_REAL;
      pctdp(3,2) = NA_REAL;
    }
    while(pctdp(i,4) < 0.5 && i < m-1) { i += 1;}
    if(pctdp(i,4) >= 0.5)
    {
      pctdp(1,0) = (i-(pctdp(i,4)-0.5)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(1,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(1,0)/dt);
      pctdp(1,2) = 0.5;
    }
    else
    {
      pctdp(1,0) = NA_REAL;
      pctdp(1,1) = NA_REAL;
      pctdp(1,2) = NA_REAL;
    }
    while(pctdp(i,4) < pcthigh && i < m-1) { i += 1;}
    if(pctdp(i,4) >= pcthigh)
    {
      pctdp(4,0) = (i-(pctdp(i,4)-pcthigh)/(pctdp(i,4)-pctdp(i-1,4)))*dt;
      pctdp(4,1) = pctdp(i,3)-(pctdp(i,3)-pctdp(i-1,3))*(i-pctdp(4,0)/dt);
      pctdp(4,2) = pcthigh;
    }
    else
    {
      pctdp(4,0) = NA_REAL;
      pctdp(4,1) = NA_REAL;
      pctdp(4,2) = NA_REAL;
    }
  }
  return pctdp;
}

#ifdef USE_PARALLEL
struct ROMCPheat : public Worker
{
  const RMatrix<double> matPaths;
  RMatrix<double> heat;
  RMatrix<double> dens;
  std::size_t paths;
  std::size_t n;
  std::size_t offset;
  double width;
  double zmin;

  ROMCPheat(const NumericMatrix& matPaths, NumericMatrix& heat, NumericMatrix& dens, std::size_t paths, std::size_t n, std::size_t offset, double width, double zmin)
    : matPaths(matPaths), heat(heat), dens(dens), paths(paths), n(n), offset(offset), width(width), zmin(zmin) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++)
    {
      for(std::size_t j = 0; j < paths; j++)
      {
        if(!NumericVector::is_na(matPaths(i,j)))
        {
          std::size_t bin = static_cast<int>((matPaths(i,j)-zmin)/width);
          dens(i,bin) += 1;
        }
      }
      for(std::size_t j = 0; j < n; j++) { heat(i,j) = dens(i,j+offset)/(paths*width); }
    }
  }
};
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCHeatCountZ(matPaths,z)
//' @param  matPaths matrix of paths
//' @param  z        vector of states
//' @return heat(m,n) <- RcppOUPMCHeatCountZ()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCHeatCountZ(NumericMatrix matPaths, NumericVector z)
{
  std::size_t paths = matPaths.ncol();
  std::size_t m = matPaths.nrow();
  std::size_t n = z.size();
  double x = matPaths(0,0);
  double width;
  double zmin;
  double zmax;
  NumericVector minmax(2);
  if(n > 1) { width = (z[n-1]-z[0])/(n-1); }
  else { width = 1; }
  minmax = RcppOUPMCMinMax(matPaths);
  if(minmax[0] > z[0]-0.5*width) { zmin = z[0]-0.5*width; }
  else { zmin = z[0]-0.5*width-static_cast<int>((z[0]-minmax[0])/width+0.5)*width; }
  if(minmax[1] < z[0]-0.5*width+n*width) { zmax = z[0]-0.5*width+n*width; }
  else { zmax = z[0]-0.5*width+n*width+static_cast<int>((minmax[1]-z[0])/width-n+1.5)*width; }
  std::size_t nn = static_cast<int>((zmax-zmin)/width+0.5);
  std::size_t offset = static_cast<int>((z[0]-zmin)/width);
  NumericMatrix dens(m,nn);
  NumericMatrix heat(m,n);
  for(std::size_t j = 0; j < n; j++)
  {
    if(z[j] < x-0.5*width || z[j] >= x+0.5*width) { heat(0,j) = 0; }
    else { heat(0,j) = 1/width; }
  }
#ifdef USE_PARALLEL
  ROMCPheat worker(matPaths, heat, dens, paths, n, offset, width, zmin);
  parallelFor(1, m, worker);
#else
  for(std::size_t i = 1; i < m; i++)
  {
    for(std::size_t j = 0; j < paths; j++)
    {
      if(!NumericVector::is_na(matPaths(i,j)))
      {
        std::size_t bin = static_cast<int>((matPaths(i,j)-zmin)/width);
        dens(i,bin) += 1;
      }
    }
    for(std::size_t j = 0; j < n; j++) { heat(i,j) = dens(i,j+offset)/(paths*width); }
  }
#endif
  return heat;
}
