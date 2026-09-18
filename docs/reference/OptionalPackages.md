# Optional packages

Queries whether functions are compiled with RcppParallel or fall back to
Rcpp only. Also queries whether random number packages dqrng and sitmo
are installed.

## Usage

``` r
RcppParallelInstalled()

RcppdqrngInstalled()

RcppsitmoInstalled()
```

## Value

bool \<- RcppParallelInstalled()

bool \<- RcppdqrngInstalled()

bool \<- RcppsitmoInstalled()

## Discussion

Rcpp calculates hundreds of times faster than R6 objects. RcppParallel
calculates five to eight times faster than Rcpp on a typical laptop and
thousands of times faster than R6 objects. Random number generation with
the R function rnorm() is slow. The packages dqrng and sitmo are
alternatives:

     install.packages("RcppParallel", "dqrng", "sitmo")

If RcppParallel is installed it will be used for almost every
calculation. If dqrng is installed, it will be the default for random
number generation. Otherwise, the default is std::mt19937. If sitmo is
installed, it can be selected as an option in the function
RcppOUPStandardNormal().
