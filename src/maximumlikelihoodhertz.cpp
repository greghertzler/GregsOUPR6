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

//' @title MaximumLikelihood_Rcpp functions to estimate an Ornstein-Uhlenbeck Process
//'
//' @description
//' Calculations for the R6 class 'MaximumLikelihood', with sequential and parallel
//'  processing.
//'
//' @details # Notes on Values
//' Return values are vectors and matrices allocated in Rcpp.  The dimensions are
//'  shown for information.  Of course, do not include them in R calls.  For example:
//'
//'     logL <- RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
//'
//' The return value:
//'
//'     logL(4)
//'
//'  is a vector containing the LogLikelihood, NA, alpha and m-1, subset in R as:
//'
//'     logL <- RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
//'     lnL <- logL[1]
//'     k <- logL[2]
//'     alpha <- logL[3]
//'     m1 <- logL[4]
//'
//' The return value:
//'
//'     nmstart(2,4)
//'
//'  is a matrix containing the starting theta, starting steps, length of tau
//'  and number of non-zero steps.  Each non-zero step identifies a free
//'  parameter in theta.  Subset in R as:
//'
//'     nmstart <- RcppOUPMLNMStart(tau,z,rhor,mur,sigmar,rhos,mus,sigmas)
//'     thetas <- nmstart[1,1:3]
//'     stepss <- nmstart[2,1:3]
//'     m <- nmstart[1,4]
//'     nk <- nmstart[2,4]
//'
//' The return value:
//'
//'     theta(3)
//'
//'  is a vector containing the estimated parameters, subset in R as:
//'
//'     theta <- RcppOUPMLNelderMead(tau,z,thetas,steps)
//'     rho <- theta[1]
//'     mu <- theta[2]
//'     sigma <- theta[3]
//'
//' The return value:
//'
//'     gof(2,6)
//'
//'  is a matrix containing two vectors for Invariant and Scaled Brownian Motion
//'  comparisons.  In R, the results are combined into four lists:
//'
//'     gof <- RcppOUPMLGoodnessOfFit(tau,z,lnL,alpha)
//'     theta_i <- list(rhor=gof[1,1],mu=gof[1,2],sigma=gof[1,3],lnLr=gof[1,4],k=gof[1,5],alphar=gof[1,6],m1=gof[1,7])
//'     Inv <- list(R2=gof[1,8],PVal=gof[1,9])
//'     theta_s <- list(rhor=gof[2,1],mur=gof[2,2],sigma=gof[2,3],lnLr=gof[2,4],k=gof[2,5],alphar=gof[2,6],m1=gof[2,7])
//'     SBM <- list(R2=gof[2,8],Pval=gof[2,9])
//'
//' The return value:
//'
//'     lrt(2)
//'
//'  is a vector containing R2 and Pval, subset in R as:
//'
//'     lrt <- RcppOUPMLLikelihoodRatioTest(lnL,alpha,m,lnLr)
//'     R2 <- lrt[1]
//'     PVal <- lrt[2]
//'
//' @details # Discussion
//' First, maximum likelihood estimation was implemented in R6 as a
//'  single-threaded application.  Then the R6 code was translated into
//'  Rcpp sequential code and, where possible, into RcppParallel code.  The
//'  Nelder-Mead algorithm and the Likelihood Ratio Test are hopelessly sequential.
//'  The Log Likelihood and the Goodness of Fit test are amenable to parallel
//'  processing.  Therefore, functions for maximum likelihood estimation use
//'  both sequential and parallel processing.  The R6 single-threaded code was
//'  archived.
//'
//' Below are microbenchmark median times for estimating the Ornstein-Uhlenbeck
//'  Process with 19,312 observations on an i7 CPU with 12 threads running at a
//'  maximum speed of 4.5 GHz.
//'
//'     Unit: milliseconds             R6      R6+      R6+Rcpp+
//'               function  single-thread     Rcpp  RcppParallel
//'     --------------------------------------------------------
//'               Estimate        3951.54  2386.59        144.35
//'
//' R6 single-thread takes almost 4 seconds.  R6+Rcpp, with sequential processing
//'  for both the Nelder-Mead algorithm and the Log Likelihood function, calculates
//'  1.6 times faster.  R6+Rcpp+RcppParallel, with parallel processing of the
//'  Log Likelihood function, calculates calculates 16.5 times faster than R6+Rcpp
//'  and 27.4 times faster than R6 single-thread.  Times increase linearly with
//'  sample size.  The median time for estimating with 59,256 observations is
//'  0.45198 seconds.
//'
//' RccpParallel uses Intel's Threading Building Blocks (TBB) on the CPU.  Unlike
//'  parallel processing on a GPU or accelerator, memory isn't copied and there
//'  is very little overhead.  On small problems, both the parallel and sequential
//'  code calculate within a clock tick.
//'
//' RcppParallel is an optional package.  If it is installed, it will be used.
//'  Function RcppParallelInstalled() will enquire whether code is compiled with
//'  RcppParallel or has fallen back to Rcpp.
//'
//' @details # From the Console
//' These functions are available in R, the RStudio console and RShiny apps.
//'  For example, an estimation of unrestricted parameters would be:
//'
//'     df <- OUPDataRead("OUP_Convergence")
//'     tau <- df[[1]]
//'     z <- df[[2]]
//'     nmstart <- RcppOUPMLNMStart(tau,z)
//'     thetas <- nmstart[1,1:3]
//'     steps <- nmstart[2,1:3]
//'     theta <- RcppOUPMLNelderMead(tau,z,thetas,steps)
//'
//' Good starting values are calculated from the data in RcppOUPMLNMStart(), or
//'  starting values can be chosen some other way. Both thetas and steps are called
//'  by reference and modified within RcppOUPMLNelderMead().  Calling RcppOUPMLNelderMead()
//'  again, will start from where it left off with thetas and steps.  This can
//'  be used to iterate over similar data sets without calling RcppOUPMLNMStart()
//'  each time.
//'
//' Calling functions directly from the console is slightly faster than calling
//'  them indirectly through R6 objects.  Here are microbenchmark median times for
//'  the data set with 19,312 observations:
//'
//'     Unit: milliseconds       R6+Rcpp+  Console Rcpp+
//'                function  RcppParallel   RcppParallel
//'     -----------------------------------------
//'           LogLikelihood       0.56320        0.48475
//'                 NMStart       0.40560        0.11865
//'              NelderMead     144.35000      143.42000
//'           GoodnessOfFit       1.45590        1.19890
//'     LikelihoodRatioTest       0.13540        0.00440
//'
//' The extra times taken by the R6 object are fractions of a millisecond.  For a
//'  small time penalty, the R6 object is more convenient.  All inputs are optional
//'  and are coordinated across functions.  Enter an input once and calculate several
//'  outputs.  The R6 object is reactive.  In other words, it stores the inputs and
//'  outputs and maps inputs to outputs.  If an input changes, dependent outputs are
//'  nullified and will be recalculated, as requested, but nothing is calculated twice.
//'  The console stores outputs in the global environment, but there is no map to
//'  inputs and outputs can be stale.  Another advantage of the R6 object is predefined
//'  plots with Plotly. The same simulation can plotted different ways without
//'  recalculation.
//'
//' The overhead of setting up RcppParallel means that small problems will
//'  calculate more slowly.  But small problems calculate in microseconds, anyway.
//'
//' Sequential calculations are reproducible, but parallel calculations are not.
//'  Two runs of the same problem will agree to about 12 significant digits, but
//'  disagree thereafter.  For exact arithmetic, order doesn't matter.  For
//'  floating-point arithmetic, it does.  For example, a Log Likelihood is rounded
//'  to 15 digits as Log Likelihoods for each observation are added.  Changing the
//'  order changes the rounding and may give slightly different answers.  The TBB
//'  scheduler determines the order.
//'
//' Potentially, the functions could be imported into other packages.
//'
//' @name MaximumLikelihood_Rcpp

// Helpers (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

#ifdef USE_PARALLEL
struct NMlnLrhozero : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double sgma;
  double logL;

  NMlnLrhozero(const NumericVector& tau, const NumericVector& z, double sgma)
    : tau(tau), z(z), sgma(sgma), logL(0.0) {}
  NMlnLrhozero(const NMlnLrhozero& banana, Split)
    : tau(banana.tau), z(banana.z), sgma(banana.sgma), logL(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_logL = logL;
    for(std::size_t i = begin; i < end; i++)
    {
      double mean = z[i];
      double variance = pow(sgma,2)*(tau[i+1]-tau[i]);
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      local_logL += lnu-v;
    }
    logL = local_logL;
  }
  void join(const NMlnLrhozero& rhs) {
    logL += rhs.logL;
  }
};

struct NMlnL : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double rho;
  double mu;
  double sgma;
  double logL;

  NMlnL(const NumericVector& tau, const NumericVector& z, double rho, double mu, double sgma)
    : tau(tau), z(z), rho(rho), mu(mu), sgma(sgma), logL(0.0) {}
  NMlnL(const NMlnL& banana, Split)
    : tau(banana.tau), z(banana.z), rho(banana.rho), mu(banana.mu), sgma(banana.sgma), logL(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_logL = logL;
    for(std::size_t i = begin; i < end; i++)
    {
      double mean = mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]));
      double variance = pow(sgma,2)/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])));
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      local_logL += lnu-v;
    }
    logL = local_logL;
  }
  void join(const NMlnL& rhs) {
    logL += rhs.logL;
  }
};
#endif

double NMLogLikelihood(NumericVector tau, NumericVector z, double rho, double mu, double sigma)
{
  std::size_t m = tau.size();
  double sgma = sigma;
  if(sgma < 0.000001) { sgma = 0.000001; }
  double logL = 0.0;
  if(abs(rho) < 0.0000000001)
  {
#ifdef USE_PARALLEL
    NMlnLrhozero worker(tau,z,sgma);
    parallelReduce(0,m-1,worker);
    logL = worker.logL;
#else
    for(std::size_t i = 0; i < m-1; i++)
    {
      double mean = z[i];
      double variance = pow(sgma,2)*(tau[i+1]-tau[i]);
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      logL += lnu-v;
    }
#endif
  }
  else
  {
#ifdef USE_PARALLEL
    NMlnL worker(tau,z,rho,mu,sgma);
    parallelReduce(0,m-1,worker);
    logL = worker.logL;
#else
    for(std::size_t i = 0; i < m-1; i++)
    {
      double mean = mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]));
      double variance = pow(sgma,2)/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])));
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      logL += lnu-v;
    }
#endif
  }
  return logL;
}

// Exports (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

#ifdef USE_PARALLEL
struct ROMLPlnLrhozero : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double sgma;
  double logL;

  ROMLPlnLrhozero(const NumericVector& tau, const NumericVector& z, double sgma)
    : tau(tau), z(z), sgma(sgma), logL(0.0) {}
  ROMLPlnLrhozero(const ROMLPlnLrhozero& banana, Split)
    : tau(banana.tau), z(banana.z), sgma(banana.sgma), logL(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_logL = logL;
    for(std::size_t i = begin; i < end; i++)
    {
      double mean = z[i];
      double variance = pow(sgma,2)*(tau[i+1]-tau[i]);
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      local_logL += lnu-v;
    }
    logL = local_logL;
  }
  void join(const ROMLPlnLrhozero& rhs) {
    logL += rhs.logL;
  }
};

struct ROMLPlnL : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double rho;
  double mu;
  double sgma;
  double logL;
  double alpha;

  ROMLPlnL(const NumericVector& tau, const NumericVector& z, double rho, double mu, double sgma)
    : tau(tau), z(z), rho(rho), mu(mu), sgma(sgma), logL(0.0), alpha(0.0) {}
  ROMLPlnL(const ROMLPlnL& banana, Split)
    : tau(banana.tau), z(banana.z), rho(banana.rho), mu(banana.mu), sgma(banana.sgma), logL(0.0), alpha(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_logL = logL;
    double local_alpha = alpha;
    for(std::size_t i = begin; i < end; i++)
    {
      double mean = mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]));
      double variance = pow(sgma,2)/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])));
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      local_logL += lnu-v;
      local_alpha += 1+exp(-2*rho*(tau[i+1]-tau[i]));
    }
    logL = local_logL;
    alpha = local_alpha;
  }
  void join(const ROMLPlnL& rhs) {
    logL += rhs.logL;
    alpha += rhs.alpha;
  }
};
#endif

//' @rdname MaximumLikelihood_Rcpp
//' @usage  RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
//' @param  tau   vector of times
//' @param  z     vector of states
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @return logL(4) <- RcppOUPMLLogLikelihood()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPMLLogLikelihood(NumericVector tau, NumericVector z, double rho, double mu, double sigma)
{
  std::size_t m = tau.size();
  double sgma = sigma;
  if(sgma < 0.000001) { sgma = 0.000001; }
  NumericVector logL(4);
  logL[1] = NA_REAL;
  logL[3] = m-1;
  if(abs(rho) < 0.0000000001)
  {
#ifdef USE_PARALLEL
    ROMLPlnLrhozero worker(tau,z,sgma);
    parallelReduce(0,m-1,worker);
    logL[0] = worker.logL;
    logL[2] = 1.0;
#else
    for(std::size_t i = 0; i < m-1; i++)
    {
      double mean = z[i];
      double variance = pow(sgma,2)*(tau[i+1]-tau[i]);
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      logL[0] += lnu-v;
    }
    logL[2] = 0.5;
#endif
  }
  else
  {
#ifdef USE_PARALLEL
    ROMLPlnL worker(tau,z,rho,mu,sgma);
    parallelReduce(0,m-1,worker);
    logL[0] = worker.logL;
    logL[2] = 0.5*worker.alpha/(m-1);
#else
    for(std::size_t i = 0; i < m-1; i++)
    {
      double mean = mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]));
      double variance = pow(sgma,2)/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])));
      double lnu = -0.5*log(2*3.14159265358979*variance);
      double v = 0.5*pow(z[i+1]-mean,2)/variance;
      logL[0] += lnu-v;
      logL[2] += 1+exp(-2*rho*(tau[i+1]-tau[i]));
    }
    logL[2] *= 0.5/(m-1);
#endif
  }
  return logL;
}

#ifdef USE_PARALLEL
struct ROMLPNMS : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double sum;
  double invsq;
  double sbmsq;

  ROMLPNMS(const NumericVector& tau, const NumericVector& z)
    : tau(tau), z(z), sum(0.0), invsq(0.0), sbmsq(0.0) {}
  ROMLPNMS(const ROMLPNMS& banana, Split)
    : tau(banana.tau), z(banana.z), sum(0.0), invsq(0.0), sbmsq(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_sum = sum;
    double local_invsq = invsq;
    double local_sbmsq = sbmsq;
    for(std::size_t i = begin; i < end; i++)
    {
      local_sum += z[i+1];
      local_invsq += z[i+1]*z[i+1];
      local_sbmsq += (z[i+1]-z[i])*(z[i+1]-z[i])/(tau[i+1]-tau[i]);
    }
    sum = local_sum;
    invsq = local_invsq;
    sbmsq = local_sbmsq;
  }
  void join(const ROMLPNMS& rhs) {
    sum += rhs.sum;
    invsq += rhs.invsq;
    sbmsq += rhs.sbmsq;
  }
};
#endif

//' @rdname MaximumLikelihood_Rcpp
//' @usage  RcppOUPMLNMStart(tau,z,rhor,mur,sigmar,rhos,mus,sigmas)
//' @param  tau    vector of times
//' @param  z      vector of states
//' @param  rhor   optional constant to fix the rate parameter 0<=rhor<inf
//' @param  mur    optional constant to fix the location parameter -inf<mur<inf
//' @param  sigmar optional constant to fix scale parameter -inf<sigmar<inf
//' @param  rhos   optional starting value for the rate parameter 0<=rhos<inf
//' @param  mus    optional starting value for the location parameter -inf<mus<inf
//' @param  sigmas optional starting value the scale parameter -inf<sigmas<inf
//' @return start(2,4) <- RcppOUPMLNMStart()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMLNMStart(NumericVector tau, NumericVector z, Nullable<double> rhor=R_NilValue, Nullable<double> mur=R_NilValue, Nullable<double> sigmar=R_NilValue, Nullable<double> rhos=R_NilValue, Nullable<double> mus=R_NilValue, Nullable<double> sigmas=R_NilValue)
{
  std::size_t m = tau.size();
  NumericMatrix thetasteps(2,4);
#ifdef USE_PARALLEL
  ROMLPNMS worker(tau,z);
  parallelReduce(0,m-1,worker);
  double sum = worker.sum/(m-1);
  double invsq = worker.invsq/(m-1)-sum*sum;
  double sbmsq = worker.sbmsq/(m-1);
#else
  double sum = 0.0;
  double invsq = 0.0;
  double sbmsq = 0.0;
  for(std::size_t i = 0; i < m-1; i++)
  {
    sum += z[i+1];
    invsq += z[i+1]*z[i+1];
    sbmsq += (z[i+1]-z[i])*(z[i+1]-z[i])/(tau[i+1]-tau[i]);
  }
  sum /= (m-1);
  invsq = invsq/(m-1)-sum*sum;
  sbmsq /= (m-1);
#endif
  if(sigmar.isNotNull())
  {
    thetasteps(0,2) = as<double>(sigmar);
    thetasteps(1,2) = 0;
  }
  else
  {
    if(sigmas.isNotNull()) { thetasteps(0,2) = as<double>(sigmas); }
    else { thetasteps(0,2) = pow(sbmsq,0.5); }
    thetasteps(1,2) = 1;
  }
  if(mur.isNotNull())
  {
    thetasteps(0,1) = as<double>(mur);
    thetasteps(1,1) = 0;
  }
  else
  {
    if(mus.isNotNull()) { thetasteps(0,1) = as<double>(mus); }
    else { thetasteps(0,1) = sum; }
    thetasteps(1,1) = 1;
  }
  if(rhor.isNotNull())
  {
    thetasteps(0,0) = as<double>(rhor);
    if(thetasteps(0,0) <= 0.0)
    {
      thetasteps(0,0) = 0.0;
      thetasteps(0,1) = 0.0;
      thetasteps(1,1) = 0;
    }
    thetasteps(1,0) = 0;
  }
  else
  {
    if(rhos.isNotNull())
    {
      thetasteps(0,0) = as<double>(rhos);
      if(thetasteps(0,0) <= 0) { thetasteps(0,0) = 0.0; }
    }
    else { thetasteps(0,0) = 0.5*sbmsq/invsq; }
    thetasteps(1,0) = 1;
  }
  thetasteps(0,3) = m;
  thetasteps(1,3) = thetasteps(1,0)+thetasteps(1,1)+thetasteps(1,2);
  return thetasteps;
}

//' @rdname MaximumLikelihood_Rcpp
//' @usage  RcppOUPMLNelderMead(tau,z,theta,steps)
//' @param  tau   vector of times
//' @param  z     vector of states
//' @param  theta vector containing rho, mu and sigma
//' @param  steps vector of 0's or 1's, 0 for fixed, 1 for variable
//' @return theta(3) <- RcppOUPMLNelderMead()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPMLNelderMead(NumericVector tau, NumericVector z, NumericVector& theta, NumericVector& steps)
{
  // algorithm parameters
  double rho = 1.5;       // reflect
  double epsilon = 2.0;   // expand
  double mu = 1.1;        // move (not in the usual Nelder-Mead algorithm)
  double chi = 0.5;       // contract
  double sigma = 0.5;     // shrink
  double kappa = 0.1;     // minimum step size (not in the usual Nelder-Mead algorithm)
  double iota = 0.3;      // steps increment (not in the usual Nelder-Mead algorithm)
  // iteration parameters
  int sigdig = 12;
  int cmax = 9999;
  double tensig = pow(0.1,sigdig);
  // cement constant thetas and create index of thetas in the simplex
  int k = theta.size();
  IntegerVector Ix(k);
  NumericVector tj(k);
  NumericVector tbar(k);
  NumericVector tr(k);
  NumericVector te(k);
  NumericVector tm(k);
  NumericVector tc(k);
  NumericVector ts(k);
  int nk = 0;
  for(int i = 0; i < k; i++)
  {
    if(steps[i] == 0)
    {
      tj[i] = theta[i];
      tbar[i] = theta[i];
      tr[i] = theta[i];
      te[i] = theta[i];
      tm[i] = theta[i];
      tc[i] = theta[i];
      ts[i] = theta[i];
    }
    else
    {
      steps[i] = theta[i]/iota;
      if(steps[i] < 0.1) { steps[i] = 0.1; }
      Ix[nk] = i;
      nk += 1;
    }
  }
// Rcout << "steps" << std::endl;
// Rcout << steps << std::endl;
// Rcout << "Ix" << std::endl;
// Rcout << Ix << std::endl;
  // outer loop
  double LnL = -NMLogLikelihood(tau,z,theta[0],theta[1],theta[2]);
  if(nk > 0)
  {
    NumericMatrix tplex(nk,nk+1);
    NumericVector Lplex(nk+1);
    double Lprev = std::numeric_limits<double>::max();
    int sign = -1;
    int cnt = 0;
    int starts = 0;
    while(abs(Lprev-LnL) > abs(LnL*tensig) && cnt < cmax)
    {
      starts += 1;
      sign *= -1;
      // initial simplex
      for(int i = 0; i < nk; i++)
      {
        tplex(i,0) = theta[Ix[i]];
      }
      Lplex[0] = LnL;
      for(int j = 1; j < nk+1; j++)
      {
        for(int i = 0; i < j-1; i++)
        {
          tplex(i,j) = theta[Ix[i]];
          tj[Ix[i]] = tplex(i,j);
        }
        tplex(j-1,j) = theta[Ix[j-1]]+sign*steps[Ix[j-1]];
        tj[Ix[j-1]] = tplex(j-1,j);
        for(int i = j; i < nk; i++)
        {
          tplex(i,j) = theta[Ix[i]];
          tj[Ix[i]] = tplex(i,j);
        }
        Lplex[j] = -NMLogLikelihood(tau,z,tj[0],tj[1],tj[2]);
      }
      // inner loop
      int jmin;
      int jmax;
      double tdev = std::numeric_limits<double>::max();
      double Ldev = std::numeric_limits<double>::max();
      while(tdev > tensig && Ldev > tensig && cnt < cmax)
      {
// Rcout << "\n" << "tplex" << "\n" << tplex << "Lplex" << "\n" << Lplex << std::endl;
        //   minimum and maximum function values
        double Lmin = Lplex[0];
        jmin = 0;
        double Lmax = Lplex[0];
        jmax = 0;
        for(int j = 1; j < nk+1; j++)
        {
          if(Lplex[j] < Lmin)
          {
            Lmin = Lplex[j];
            jmin = j;
          }
          else if(Lplex[j] > Lmax)
          {
            Lmax = Lplex[j];
            jmax = j;
          }
        }
// Rcout << "Lmin,jmin,Lmax,jmax\n" << Lmin << "," << jmin << "," << Lmax << "," << jmax << std::endl;
        //  penultimate function value
        double Lpenult = Lmin;
        for(int j = 0; j < nk+1; j++)
        {
          if(Lpenult < Lplex[j] && Lplex[j] < Lmax)
          {
            Lpenult = Lplex[j];
          }
        }
// Rcout << "Lpenult\n" << Lpenult << std::endl;
        //   calculate centroid of all but Lmax
        for(int i = 0; i < nk; i++)
        {
          tbar[Ix[i]] = 0;
          for(int j = 0; j < jmax; j++)
          {
            tbar[Ix[i]] = tbar[Ix[i]]+tplex(i,j);
          }
          for(int j = jmax+1; j < nk+1; j++)
          {
            tbar[Ix[i]] = tbar[Ix[i]]+tplex(i,j);
          }
          tbar[Ix[i]] = tbar[Ix[i]]/nk;
        }
        double Lbar = -NMLogLikelihood(tau,z,tbar[0],tbar[1],tbar[2]);
// Rcout << "tbar\n" << tbar << "\nLbar\n" << Lbar << std::endl;
        //    calculate reflection of Lmax
        for(int i = 0; i < nk; i++)
        {
          tr[Ix[i]] = tbar[Ix[i]]+rho*(tbar[Ix[i]]-tplex(i,jmax));
        }
        double Lr = -NMLogLikelihood(tau,z,tr[0],tr[1],tr[2]);
// Rcout << "Lr\n" << Lr << std::endl;
        //    expansion
        if(Lr < Lmin)
        {
          for(int i = 0; i < nk; i++)
          {
            te[Ix[i]] = tbar[Ix[i]]+epsilon*(tr[Ix[i]]-tbar[Ix[i]]);
          }
          double Le = -NMLogLikelihood(tau,z,te[0],te[1],te[2]);
          if(Le < Lr)
          {
            for(int i = 0; i < nk; i++)
            {
              tplex(i,jmax) = te[Ix[i]];
            }
            Lplex[jmax] = Le;
          }
          else
          {
            for(int i = 0; i < nk; i++)
            {
              tplex(i,jmax) = tr[Ix[i]];
            }
            Lplex[jmax] = Lr;
          }
// Rcout << "expansion, Lr, Lmin, Le " << Lr << "," << Lmin << "," << Le << std::endl;
        }
        //    reflection
        else if(Lmin <= Lr && Lr < Lpenult)
        {
          for(int i = 0; i < nk; i++)
          {
            tplex(i,jmax) = tr[Ix[i]];
            tm[Ix[i]] = tbar[Ix[i]]+mu*(tplex(i,jmin)-tbar[Ix[i]]);
          }
          Lplex[jmax] = Lr;
// Rcout << "reflection, Lr, Lmin, Lpenult " << Lr << "," << Lmin << "," << Lpenult << std::endl;
          //    move minimum
          double Lm = -NMLogLikelihood(tau,z,tm[0],tm[1],tm[2]);
          if(Lm < Lplex[jmin])
          {
            for(int i = 0; i < nk; i++)
            {
              tplex(i,jmin) = tm[Ix[i]];
            }
            Lplex[jmin] = Lm;
// Rcout << "move, Lr, Lmin, Lm " << Lr << "," << Lmin << "," << Lm << std::endl;
          }
        }
        //    outside contraction
        else if(Lr < Lmax)
        {
          for(int i = 0; i < nk; i++)
          {
            tc[Ix[i]] = tbar[Ix[i]]+chi*(tr[Ix[i]]-tbar[Ix[i]]);
          }
          double Lc = -NMLogLikelihood(tau,z,tc[0],tc[1],tc[2]);
          if(Lc < Lr)
          {
            for(int i =0; i < nk; i++)
            {
              tplex(i,jmax) = tc[Ix[i]];
            }
            Lplex[jmax] = Lc;
// Rcout << "outside contraction, Lr, Lmax, Lc " << Lr << "," << Lmax << "," << Lc << std::endl;
          }
          //   shrink
          else
          {
            for(int j = 0; j < jmin; j++)
            {
              for(int i = 0; i < nk; i++)
              {
                ts[Ix[i]] = tplex(i,jmin)+sigma*(tplex(i,j)-tplex(i,jmin));
                tplex(i,j) = ts[Ix[i]];
              }
              Lplex[j] = -NMLogLikelihood(tau,z,ts[0],ts[1],ts[2]);
            }
            for(int j = jmin+1; j < nk+1; j++)
            {
              for(int i = 0; i < nk; i++)
              {
                ts[Ix[i]] = tplex(i,jmin)+sigma*(tplex(i,j)-tplex(i,jmin));
                tplex(i,j) = ts[Ix[i]];
              }
              Lplex[j] = -NMLogLikelihood(tau,z,ts[0],ts[1],ts[2]);
            }
// Rcout << "outside contraction shrink, Lr, Lmax, Lc " << Lr << "," << Lmax << "," << Lc << std::endl;
          }
        }
        //    inside contraction
        else
        {
          for(int i = 0; i < nk; i++)
          {
            tc[Ix[i]] = tbar[Ix[i]]+chi*(tplex(i,jmax)-tbar[Ix[i]]);
          }
          double Lc = -NMLogLikelihood(tau,z,tc[0],tc[1],tc[2]);
          if(Lc < Lmax)
          {
            for(int i = 0; i < nk; i++)
            {
              tplex(i,jmax) = tc[Ix[i]];
            }
            Lplex[jmax] = Lc;
// Rcout << "inside contraction, Lr, Lmax, Lc " << Lr << "," << Lmax << "," << Lc << std::endl;
          }
          //   shrink
          else
          {
            for(int j = 0; j < jmin; j++)
            {
              for(int i = 0; i < nk; i++)
              {
                ts[Ix[i]] = tplex(i,jmin)+sigma*(tplex(i,j)-tplex(i,jmin));
                tplex(i,j) = ts[Ix[i]];
              }
              Lplex[j] = -NMLogLikelihood(tau,z,ts[0],ts[1],ts[2]);
            }
            for(int j = jmin+1; j < nk+1; j++)
            {
              for(int i = 0; i < nk; i++)
              {
                ts[Ix[i]] = tplex(i,jmin)+sigma*(tplex(i,j)-tplex(i,jmin));
                tplex(i,j) = ts[Ix[i]];
              }
              Lplex[j] = -NMLogLikelihood(tau,z,ts[0],ts[1],ts[2]);
            }
// Rcout << "inside contraction shrink, Lr, Lmax, Lc " << Lr << "," << Lmax << "," << Lc << std::endl;
          }
        }
        //   deviations
        tdev = 0;
        Ldev = 0;
        for(int j = 0; j < nk+1; j++)
        {
          for(int i = 0; i < nk; i++)
          {
            if(abs((tplex(i,j)-tbar[Ix[i]])/steps[Ix[i]]) > tdev)
            {
              tdev = abs((tplex(i,j)-tbar[Ix[i]])/steps[Ix[i]]);
            }
          }
          if(abs(Lplex[j]/Lbar-1) > Ldev)
          {
            Ldev = abs(Lplex[j]/Lbar-1);
          }
        }
        cnt += +1;
// Rcout << "tdev, Ldev, starts, cnt " << tdev << "," << Ldev << "," << starts << "," << cnt << std::endl;
      }
      // new minimum
      Lprev = LnL;
      for(int i = 0; i < nk; i++)
      {
        theta[Ix[i]] = tplex(i,jmin);
        steps[Ix[i]] = abs(theta[Ix[i]])*iota;
        if(steps[i] < kappa) { steps[i] = kappa; }
      }
      LnL = Lplex[jmin];
    }
    Rcout << "Nelder-Mead " << starts << ":" << cnt << std::endl;
  }
  return(theta);
}

#ifdef USE_PARALLEL
struct ROMLPGoF : public Worker
{
  const RVector<double> tau;
  const RVector<double> z;
  double sum;
  double invsq;
  double sbmsq;
  double sumlntau;

  ROMLPGoF(const NumericVector& tau, const NumericVector& z)
    : tau(tau), z(z), sum(0.0), invsq(0.0), sbmsq(0.0), sumlntau(0.0) {}
  ROMLPGoF(const ROMLPGoF& banana, Split)
    : tau(banana.tau), z(banana.z), sum(0.0), invsq(0.0), sbmsq(0.0), sumlntau(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    double local_sum = sum;
    double local_invsq = invsq;
    double local_sbmsq = sbmsq;
    double local_sumlntau = sumlntau;
    for(std::size_t i = begin; i < end; i++)
    {
      local_sum += z[i+1];
      local_invsq += z[i+1]*z[i+1];
      local_sbmsq += (z[i+1]-z[i])*(z[i+1]-z[i])/(tau[i+1]-tau[i]);
      local_sumlntau += log(tau[i+1]-tau[i]);
    }
    sum = local_sum;
    invsq = local_invsq;
    sbmsq = local_sbmsq;
    sumlntau = local_sumlntau;
  }
  void join(const ROMLPGoF& rhs) {
    sum += rhs.sum;
    invsq += rhs.invsq;
    sbmsq += rhs.sbmsq;
    sumlntau += rhs.sumlntau;
  }
};
#endif

//' @rdname MaximumLikelihood_Rcpp
//' @usage  RcppOUPMLGoodnessOfFit(tau,z,lnL,alpha)
//' @param  tau   vector of times
//' @param  z     vector of states
//' @param  lnL   log likelihood of estimates lnL<0
//' @param  alpha identifies the distribution of lnL 0.5<alpha<1
//' @return gof(2,6) <- RcppOUPMLGoodnessOfFit()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPMLGoodnessOfFit(NumericVector tau, NumericVector z, double lnL, double alpha)
{
  std::size_t m = tau.size();
  NumericMatrix gof(2,9);
#ifdef USE_PARALLEL
  ROMLPGoF worker(tau,z);
  parallelReduce(0,m-1,worker);
  double sum = worker.sum/(m-1);
  double invsq = worker.invsq/(m-1)-sum*sum;
  double sbmsq = worker.sbmsq/(m-1);
  double sumlntau = worker.sumlntau;
#else
  double sum = 0.0;
  double invsq = 0.0;
  double sbmsq = 0.0;
  double sumlntau = 0.0;
  for(std::size_t i = 0; i < m-1; i++)
  {
    sum += z[i+1];
    invsq += z[i+1]*z[i+1];
    sbmsq += (z[i+1]-z[i])*(z[i+1]-z[i])/(tau[i+1]-tau[i]);
    sumlntau += log(tau[i+1]-tau[i]);
  }
  sum /= (m-1);
  invsq = invsq/(m-1)-sum*sum;
  sbmsq /= (m-1);
#endif
  double lnLinv = -0.5*(m-1)*(log(2*3.14159265358979*invsq)+1);
  double upsinv = 2*(lnL-lnLinv);
  double lnLsbm = -0.5*(m-1)*(log(2*3.14159265358979*sbmsq)+1)-0.5*sumlntau;
  double upssbm = 2*(lnL-lnLsbm);
  gof(0,0) = NA_REAL;
  gof(0,1) = sum;
  gof(0,2) = pow(invsq,0.5);
  gof(0,3) = lnLinv;
  gof(0,4) = 2;
  gof(0,5) = 0.5;
  gof(0,6) = m-1;
  if(upsinv < 0)
  {
    gof(0,7) = NA_REAL;
    gof(0,8) = NA_REAL;
  }
  else
  {
    gof(0,7) = 1-exp(-log(2)*0.5*upsinv/(alpha*(m-1)));
    gof(0,8) = GammaBigRatio(alpha*(m-1),0.5*upsinv);
  }
  gof(1,0) = 0.0;
  gof(1,1) = NA_REAL;
  gof(1,2) = pow(sbmsq,0.5);
  gof(1,3) = lnLsbm;
  gof(1,4) = 1;
  gof(1,5) = 1.0;
  gof(1,6) = m-1;
  if(upssbm < 0)
  {
    gof(1,7) = NA_REAL;
    gof(1,8) = NA_REAL;
  }
  else
  {
    gof(1,7) = 1-exp(-log(2)*0.5*upssbm/(alpha*(m-1)));
    gof(1,8) = GammaBigRatio(alpha*(m-1),0.5*upssbm);
  }
  return gof;
}

//' @rdname MaximumLikelihood_Rcpp
//' @usage  RcppOUPMLLikelihoodRatioTest(lnL,alpha,m,lnLr)
//' @param  lnL   log likelihood of estimates lnL<0
//' @param  alpha identifies the distribution of lnL 0.5<alpha<1
//' @param  m     number of times t in the data set
//' @param  lnLr  the log likelihood for restricted parameters lnLr<=lnL
//' @return lrt(2) <- RcppOUPMLLikelihoodRatioTest()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPMLLikelihoodRatioTest(double lnL, double alpha, double m, double lnLr)
{
  NumericVector lrt(2);
  double upsilon = 2*(lnL-lnLr);
  if(upsilon < 0)
  {
    lrt[0] = NA_REAL;
    lrt[1] = NA_REAL;
  }
  else
  {
    lrt[0] = 1-exp(-log(2)*0.5*upsilon/(alpha*(m-1)));
    lrt[1] = GammaBigRatio(alpha*(m-1),0.5*upsilon);
  }
  return lrt;
}
