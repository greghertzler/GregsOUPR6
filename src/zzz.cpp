#include <Rcpp.h>
using namespace Rcpp;

// roxygen (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @title Optional packages
//'
//' @description
//' Queries whether functions are compiled with RcppParallel or fall back
//'  to Rcpp only.  Also queries whether random number packages dqrng and sitmo
//'  are installed.
//'
//' @details # Discussion
//' Rcpp calculates hundreds of times faster than R6 objects.  RcppParallel
//'  calculates five to eight times faster than Rcpp on a typical laptop and
//'  thousands of times faster than R6 objects.  Random number generation with
//'  the R function rnorm() is slow.  The packages dqrng and sitmo are
//'  alternatives:
//'
//'      install.packages("RcppParallel", "dqrng", "sitmo")
//'
//' If RcppParallel is installed it will be used for almost every calculation.
//'  If dqrng is installed, it will be the default for random number generation.
//'  Otherwise, the default is std::mt19937.  If sitmo is installed, it can be
//'  selected as an option in the function RcppOUPStandardNormal().
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
//' @usage RcppdqrngInstalled()
//' @return bool <- RcppdqrngInstalled()
//' @export
// [[Rcpp::export]]
bool RcppdqrngInstalled()
{
#ifdef USE_DQRNG
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
