#include <Rcpp.h>
using namespace Rcpp;
#ifdef USE_DQRNG
#include <dqrng.h>
#include <dqrng_distribution.h>
#endif
#if defined(USE_SITMO)
#include <sitmo.h>
#endif
#ifdef USE_PARALLEL
#include <RcppParallel.h>
using namespace RcppParallel;
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
//'     stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed,engine)
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
//'     bndfpt <- RcppOUPMCBoundedPaths(stdnorm,k,x,m,skip,dt,rho,mu,sigma,method)
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
//'  which locks up the computer.  Rccp can be hundreds of times faster and makes
//'  Monte-Carlo simulations practical for interactive applications such as RStudio
//'  and RShiny. RcppParallel speeds the calculations another five to eight times.
//'
//' Monte Carlo simulations are possible paths taken by the stochastic integral
//'  equation.  Paths can go forward or backward and may be bounded.  Paths are
//'  binned and counted to approximate several solutions. The approximations
//'  converge to analytical solutions as the number of paths increases.  Binning
//'  and counting 1,000,000 paths will be accurate to 3 or 4 significant digits.
//'
//' Here are microbenchmark median times for 100,000 and 1,000,000 paths over
//'  100 time intervals, as calculated by  R6+RccpParallel on an i7 CPU with
//'  12 threads running at a maximum of 4.5 GHz:
//'
//'     Unit: milliseconds     paths                paths
//'               function   100,000            1,000,000
//'     ------------------------------------------------------------
//'         StandardNormal   21.6702             218.0919
//'           ForwardPaths   19.6883  ________   189.7758  _________
//'               Subtotal             41.4642              407.8677
//'            Probability   54.8138  ________   965.1143  _________
//'                  Total             96.2780             1372.9820
//'
//' The R6 object is reactive and will call the StandardNormal function only once.
//'  After that the standard normal variables will be passed to ForwardPaths and
//'  the forward paths will be passed to Probability.  To calculate a median time,
//'  the R6 object is tricked into recalculating by changing an input, calculating,
//'  changing the input back to the original, and calculating the original again.
//'  Do this 11 times and record the sixth fastest time as the median.  Times for
//'  ForwardPaths, the Subtotal and the Total are calculating by trickery. Times
//'  for StandardNormal and Probability are inferred by subtraction.
//'
//' Times for StandardNormal and ForwardPaths go up approximately ten-fold with
//'  a ten-fold increase in paths.  Times for Probability go up almost 18-fold
//'  with a ten-fold increase in paths.
//'
//' The R6 function Probability calls the Rcpp function ForwardCountY which bins
//'  and counts means, variances, transition densities, transition probabilities
//'  and  double integrals.  So five sets of plots can be drawn from one simulation
//'  followed by a count.
//'
//' The R6 object manages inputs and outputs and draws plots.  All calculations
//'  are in Rcpp and RcppParallel functions.  RccpParallel uses Intel's Threading
//'  Building Blocks (TBB) on the CPU.  Unlike parallel processing on a GPU or
//'  accelerator, memory isn't copied and there is less overhead.  On trivially
//'  small problems, sequential versions calculate faster.  On large problems,
//'  parallel versions calculate much faster.
//'
//' RcppParallel is an optional package.  If it is installed, it will be used.
//'  Function RcppParallelInstalled() will enquire whether code is compiled with
//'  RcppParallel or has fallen back to Rcpp.  Function RcppParallelThreads will
//'  return the number of threads. Optional packages for random number generation
//'  are dqrng and sitmo.  The functions RcppdqrngInstalled() and
//'  RcppsitmoInstalled() will enquire whether they are installed.
//'
//' @details # From the Console
//' Rcpp and RcppParallel functions are available in R, the RStudio console and
//'  RShiny apps.  From the console, a simulation of 1,000,000 forward paths
//'  over 100 time intervals would be:
//'
//'      stdnorm <- RcppOUPMCStandardNormal(101,1,1000000,9999,1)
//'      fwd <- RcppOUPMCForwardPaths(stdnorm,15,101,1,0.1,0.5,-15,15,5)
//'
//' The R6 object doesn't give users a choice, but from the console there are four
//'  random number generators available: dqrng, std::mt19937, sitmo, and rnorm.
//'  In the last argument of the function, these are requested as engine
//'  1, 2, 3 or 4, respectively.  Microbenchmark median times for 100,000,000
//'  standard normal variables are:
//'
//'     Unit: milliseconds
//'               language        engine    transform  StandardNormal
//'     -------------------------------------------------------------
//'                   Rcpp         rnorm    inversion       4145.0140
//'                   Rcpp   sitmo::prng   Box-Muller       3937.8950
//'                   Rcpp  std::mt19937   Box-Muller       3231.4860
//'                   Rcpp  dqrng::pcg64     Ziggurat        727.3865
//'           RcppParallel   sitmo::prng   Box-Muller        585.1609
//'           RcppParallel  std::mt19937   Box-Muller        547.9659
//'           RcppParallel  dqrng::pcg64     Ziggurat        247.1802
//'
//' The random number generators generate uniform random variables.  More time is
//'  spent transforming uniform to normal random variables.  The method of transform
//'  is also listed.  These include inverting the cumulative normal, the Polar
//'  Box-Muller transform and the Ziggurat transform.  The Ziggurat transform is
//'  a sophisticated lookup table and much faster.  Results on your computer may
//'  vary. On Unix-alike operating systems, the Polar Box-Muller transform has been
//'  replaced with the Ziggurat transform.
//'
//' Both dqrng and sitmo have other random number generators, but dqrng::pcg64 and
//'  sitmo::prng are the defaults.  Both packages are optional.  If dqrng is not
//'  installed, the fall back is std::mt19937, which is always available.  Therefore,
//'  sitmo::prng and rnorm will only used if requested as engines 3 and 4.  If
//'  sitmo::prng is requested but not installed the fallback is rnorm.
//'
//' Another choice available from the console is the method of simulation, either
//'  a 4th-order Runge-Kutta numerical integration or the stochastic integral
//'  equation, itself.  The last argument of the function is the method, with 4
//'  for 4th-order Runge-Kutta and 5 for stochastic integral equation.  Arguments
//'  1, 2 and 3 are reserved for possible future implementations of 1st-order
//'  Euler and 2nd and 3rd order Maryuma or Runge-Kutta methods.  But this could
//'  be dangerous. Users might use them.  The purpose would be to demonstrate that
//'  low-order numerical methods only converge with short time intervals.
//'
//' In the function above, the skip argument divides the time intervals.  For example,
//'  if the number of times is 101, there are 100 time intervals.  Argument
//'  skip=10 subdivides 100 into 1000 time intervals for the calculations and
//'  reports results at the 101 times.
//'
//' Here are microbenchmark times for 1,000,000 paths over 100 time intervals with
//'  increasing skips.  Also shown are the maximum and minimum differences.
//'
//'     Unit: milliseconds   Standard   Integral     Runge-
//'                   skip     Normal   Equation      Kutta  max dif  min dif
//'     ---------------------------------------------------------------------
//'                      1   208.6749   190.2418   360.3870  8.4e-03  -8.9e-03
//'                      2   431.8739   214.1760   608.2229  2.1e-03  -2.1e-03
//'                      4   978.7199   295.2814  1120.7290  5.3e-04  -5.3e-04
//'                      8  2132.0220   424.4842  2086.3920  1.3e-04  -1.3e-04
//'
//' Even larger skips will calculate, but microbenchmark becomes pac man and
//'  starts chomping memory.  For skip=8, the paths are the same to within four
//'  significant digits.  But the Runge-Kutta method is much slower.  The time
//'  for the standard normal variables plus the time fore the Runge-Kutta
//'  simulation is 4.2 seconds.  The integral equation is not improved by larger
//'  skips.  For skip=1, the integral equation does the job in 0.4 seconds.
//'
//' A microbenchmark comparison of indirectly calling RcppParallel functions
//'  from R6 with directly calling them from the console is:
//'
//'     Unit: milliseconds            R6+       Console
//'               function   RcppParallel  RcppParallel
//'     -----------------------------------------------
//'           ForwardPaths       540.7106      587.6297
//'          BackwardPaths       543.9790      584.8225
//'           BoundedPaths       845.6546      579.5391
//'
//' These timings are from a standing start, generating the standard normal variables
//'  before simulating the paths.  For Forward and Backward Paths, the R6 object is faster,
//'  but for Bounded Paths, it is much slower.  Bounded Paths hit the boundary and have
//'  NA values thereafter.  We might speculate that R6 is slow with NA values.
//'
//' The R6 object has advantages over the console.  All inputs are optional and are
//'  coordinated across functions.  Enter an input once and calculate several outputs.
//'  The R6 object is reactive.  In other words, it stores the inputs and outputs and
//'  maps inputs to outputs.  If an input changes, dependent outputs are nullified and
//'  will be recalculated, as requested, but nothing is calculated twice.  The console
//'  stores outputs in the global environment, but there is no map of inputs to outputs
//'  and outputs can be stale.  Another advantage of the R6 object are pre-programmed plots
//'  with Plotly.  The same simulation can be plotted different ways without recalculation.
//'
//'
//' Parallel processing has more overhead and is slower on small problems. Here
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
//'  So users get no choice.  RcppParallel functions are compiled if RcppParallel
//'  is installed. Otherwise compilation falls back to Rcpp.
//'
//' More threads may be faster but also have more overhead.  Here are microbenchmark
//'  median times by number of threads for generating standard normal variables and
//'  simulating 1,000,000 Forward Paths over 100 time intervals:
//'
//'     Unit: milliseconds
//'                threads   stdnorm  ForwardPaths      total
//'     -----------------------------------------------------
//'                      1	742.8735	    497.8163	1240.6898
//'                      2	422.2265	    307.3954	 729.6219
//'                      3	396.3184	    307.5473	 703.8657
//'                      4	257.9281	    214.4442	 472.3723
//'                      5	260.6473	    197.7865	 458.4338
//'                      6	248.4076	    197.5342	 445.9418
//'                      7	246.5812	    191.5214	 438.1026
//'                      8	231.4124	    191.8745	 423.2869
//'                      9	232.0408	    194.9968	 427.0376
//'                     10	235.8238	    191.1864	 427.0102
//'                     11	216.9428	    192.1811	 409.1239
//'                     12	211.6724	    191.6413	 403.3137
//'
//' More is better, but a few is pretty good.  You could set fewer threads using
//'  the RcppParallel commands:
//'
//'      library(RcppParallel)
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
    else if(sigma*sigma < 0.0000000001)
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
      double v2 = rho*((k-mu)/sigma)*((k-mu)/sigma);
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

#ifdef USE_PARALLEL
#ifdef USE_DQRNG
struct ROMCPSNdqrng : public RcppParallel::Worker {
  RVector<double> stdnorm;
  uint64_t seed;

  ROMCPSNdqrng(Rcpp::NumericVector& stdnorm, uint64_t seed)
    : stdnorm(stdnorm), seed(seed) {}

  void operator()(std::size_t begin, std::size_t end) {
    auto rng = dqrng::generator<>(seed+begin);
    dqrng::normal_distribution dist(0.0,1.0);
    for (std::size_t i = begin; i < end; i++) { stdnorm[i] = dist(*rng); }
  }
};
#endif

struct ROMCPSNcpp : public Worker
{
  RVector<double> stdnorm;
  uint64_t seed;

  ROMCPSNcpp(NumericVector& stdnorm, uint64_t seed)
    : stdnorm(stdnorm), seed(seed) {}

  void operator()(std::size_t begin, std::size_t end) {
    thread_local std::mt19937 rng(seed+begin);
    std::normal_distribution<double> dist(0.0,1.0);
    for(std::size_t i = begin; i < end; i++) { stdnorm[i] = dist(rng); }
  }
};

#ifdef USE_SITMO
struct ROMCPSNsitmo : public Worker
{
  RVector<double> stdnorm;
  uint64_t seed;

  ROMCPSNsitmo(NumericVector& stdnorm, uint64_t seed)
    : stdnorm(stdnorm), seed(seed) {}

  void operator()(std::size_t begin, std::size_t end) {
    sitmo::prng_engine rng(seed+begin);
    std::normal_distribution<double> dist(0.0,1.0);
    for(std::size_t i = begin; i < end; i++) { stdnorm[i] = dist(rng); }
  }
};
#endif
#endif

//' @rdname MonteCarlo_Rcpp
//' @usage  RcppOUPMCStandardNormal(m,skip,paths,seed,engine)
//' @param  m      number of rows for states over time
//' @param  skip   subdivide time interval but report every ds or dt 0<skip<20
//' @param  paths  number of columns for paths
//' @param  seed   seed for reproducibility
//' @param  engine random number generator, 1 dqrng, 2 mt19937, 3 sitmo, 4 rnorm
//' @return stdnorm((m-1)*skip,paths) <- RcppOUPMCStandardNormal()
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppOUPMCStandardNormal(int64_t m, int64_t skip, int64_t paths, uint64_t seed, uint64_t engine)
{
#ifdef USE_PARALLEL
#ifdef USE_DQRNG
  if(engine < 2)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    ROMCPSNdqrng worker(stdnorm,seed);
    parallelFor(0,n,worker);
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#endif
  if(engine < 3)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    ROMCPSNcpp worker(stdnorm,seed);
    parallelFor(0,n,worker);
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#ifdef USE_SITMO
  if(engine < 4)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    ROMCPSNsitmo worker(stdnorm,seed);
    parallelFor(0,n,worker);
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#endif
#else
#ifdef USE_DQRNG
  if(engine < 2)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    auto rng = dqrng::generator<>(seed);
    dqrng::normal_distribution dist(0.0,1.0);
    for(std::size_t i = 0; i < n; i++) { stdnorm[i] = dist(*rng); }
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#endif
  if(engine < 3)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(0.0,1.0);
    for(std::size_t i = 0; i < n; i++) { stdnorm[i] = dist(rng); }
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#ifdef USE_SITMO
  if(engine < 4)
  {
    std::size_t n = (m-1)*skip*paths;
    NumericVector stdnorm(n);
    sitmo::prng_engine rng(seed);
    std::normal_distribution<double> dist(0.0,1.0);
    for(std::size_t i = 0; i < n; i++) { stdnorm[i] = dist(rng); }
    stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
    return as<NumericMatrix>(stdnorm);
  }
#endif
#endif
  std::size_t n = (m-1)*skip*paths;
  NumericVector stdnorm(n);
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed);
  stdnorm = rnorm(n);
  stdnorm.attr("dim") = Dimension((m-1)*skip,paths);
  return as<NumericMatrix>(stdnorm);
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
//' @usage  RcppOUPMCForwardPaths(stdnorm,x,m,skip,dt,rho,mu,sigma,method)
//' @param  stdnorm matrix of standard normal shocks
//' @param  x       initial state or vector of backward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<50
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @param  method  4 for 4th order Runge-Kutta, 5 for integral equation
//' @return forward(m,paths) <- RcppOUPMCForwardPaths()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCForwardPaths(NumericMatrix stdnorm, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma, std::size_t method)
{
  std::size_t nrows = stdnorm.nrow();
  std::size_t paths = stdnorm.ncol();
  NumericMatrix forward(m,paths);
  if(nrows == (m-1)*skip)
  {
    if(method == 4)
    {
      double dtau = dt/skip;
      double H = sigma*std::sqrt(dtau);
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
    }
    else
    {
      double dtau = dt/skip;
      double H = sigma*std::sqrt(dtau);
      if(rho > 0) { H = std::sqrt(sigma*sigma/(2*rho)*(1-std::exp(-2*rho*dtau))); }
      double exprhodt = std::exp(-rho*dtau);
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
    }
  }
  else { Rcout << "stdnorm requires " << (m-1)*skip << " rows but has " << nrows << " rows instead." << std::endl; }

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
//' @usage  RcppOUPMCBackwardPaths(stdnorm,y,m,skip,ds,rho,mu,sigma,method)
//' @param  stdnorm matrix of standard normal shocks
//' @param  y       terminal state or vector of forward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<50
//' @param  ds      time interval for terminal value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @param  method  4 for 4th order Runge-Kutta, 5 for integral equation
//' @return backward(m,paths) <- RcppOUPMCBackwardPaths()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBackwardPaths(NumericMatrix stdnorm, double y, std::size_t m, std::size_t skip, double ds, double rho, double mu, double sigma, std::size_t method)
{
  std::size_t nrows = stdnorm.nrow();
  std::size_t paths = stdnorm.ncol();
  NumericMatrix backward(m,paths);
  if(nrows == (m-1)*skip)
  {
    if(method == 4)
    {
      double dtau = ds/skip;
      double H = sigma*std::sqrt(dtau);
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
    }
    else
    {
      double dtau = ds/skip;
      double H = sigma*std::sqrt(dtau);
      if(rho > 0) { H = std::sqrt(sigma*sigma/(2*rho)*(std::exp(2*rho*dtau)-1)); }
      double exprhods = std::exp(rho*dtau);
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
    }
  }
  else { Rcout << "stdnorm requires " << (m-1)*skip << " rows but has " << nrows << " rows instead." << std::endl; }

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
//' @usage  RcppOUPMCBoundedPaths(stdnorm,k,x,m,skip,dt,rho,mu,sigma,method)
//' @param  stdnorm matrix of standard normal shocks
//' @param  k       threshold -inf<k<inf
//' @param  x       initial state or vector of backward states
//' @param  m       number of rows for states over time
//' @param  skip    subdivide time interval but report every ds or dt 0<skip<20
//' @param  dt      time interval for initial value problems
//' @param  rho     rate parameter 0<=rho<inf
//' @param  mu      location parameter -inf<mu<inf
//' @param  sigma   scale parameter -inf<sigma<inf
//' @param  method  4 for 4th order Runge-Kutta, 5 for integral equation
//' @return bndfpt(m+1,paths) <- RcppOUPMCBoundedPaths()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMCBoundedPaths(NumericMatrix stdnorm, double k, double x, std::size_t m, std::size_t skip, double dt, double rho, double mu, double sigma, std::size_t method)
{
  std::size_t nrows = stdnorm.nrow();
  std::size_t paths = stdnorm.ncol();
  NumericMatrix bndfpt(m+1,paths);
  if(nrows == (m-1)*skip)
  {
    if(method == 4)
    {
      double dtau = dt/skip;
      double H = sigma*std::sqrt(dtau);
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
    }
    else
    {
      double dtau = dt/skip;
      double H = sigma*std::sqrt(dtau);
      if(rho > 0) { H = std::sqrt(sigma*sigma/(2*rho)*(1-std::exp(-2*rho*dtau))); }
      double exprhodt = std::exp(-rho*dtau);
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
    }
  }
  else { Rcout << "stdnorm requires " << (m-1)*skip << " rows but has " << nrows << " rows instead." << std::endl; }

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
        mvdpd(i,1) += forward(i,j)*forward(i,j);
        std::size_t bin = static_cast<std::size_t>((forward(i,j)-ymin)/width);
        dens(i,bin) += 1;
      }
      mvdpd(i,0) /= paths;
      mvdpd(i,1) = mvdpd(i,1)/paths-mvdpd(i,0)*mvdpd(i,0);
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
      mvdpd(i,1) += forward(i,j)*forward(i,j);
      std::size_t bin = static_cast<int>((forward(i,j)-ymin)/width);
      dens(i,bin) += 1;
    }
    mvdpd(i,0) /= paths;
    mvdpd(i,1) = mvdpd(i,1)/paths-mvdpd(i,0)*mvdpd(i,0);
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
      double exprhods = std::exp(rho*i*ds);
      double exprhords= std::exp(-(rho+r)*i*ds);
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
    double exprhods = std::exp(rho*i*ds);
    double exprhords= std::exp(-(rho+r)*i*ds);
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
