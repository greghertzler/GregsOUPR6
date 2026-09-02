#include <Rcpp.h>
using namespace Rcpp;

// roxygen (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @title Optional packages for parallel processing
//'
//' @description
//' Queries whether functions are compiled with RcppParallel or fall back
//'  to Rcpp only.  Also queries whether random numbers are generated
//'  by RcppParallel using sitmo() or fall back to Rcpp using rnorm().
//'
//' @details # Discussion
//' Rcpp calculates hundreds of times faster than R6 objects.  RcppParallel
//'  calculates five to eight times faster than Rcpp on a typical laptop with
//'  12 threads, and thousands of times faster than R6 objects.  For Monte Carlo
//'  simulations, RcppParallel using sitmo() generates random numbers six times
//'  faster Rcpp using rnorm() and nine times faster R6 using rnorm().  Both
//'  RcppParallel and sitmo are optional but recommended:
//'
//'      install.packages("RcppParallel", "sitmo")
//'
//' @name OptionalPackages

// Exports (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @rdname OptionalPackages
//' @usage  RcppParallelInstalled()
//' @return bool <- RcppParallelInstalled()
//' @export
// [[Rcpp::export]]
bool RcppParallelInstalled()
{
#ifdef USE_PARALLEL
  return true;
#else
  return false;
#endif
}

//' @rdname OptionalPackages
//' @usage RcppsitmoInstalled()
//' @return bool <- RcppsitmoInstalled()
//' @export
// [[Rcpp::export]]
bool RcppsitmoInstalled()
{
#ifdef USE_SITMO
  return true;
#else
  return false;
#endif
}
