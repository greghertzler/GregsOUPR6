#include "gammahertz.h"
#include <Rcpp.h>
#include <cmath>
#include <limits>

// Gamma functions (((((((((((((((((((((((((((((())))))))))))))))))))))))))))))

double GammaSmallOneHalf(double x)
{
  double gamma;
  double dgm;
  double gm;
  int sigdig = 13;
  int cmax = 2000;
  if(x <= 0) { gamma = 0; }
  else if(x >= 33.8551300352457) { gamma = 1.77245385090552; }
  else
  {
    dgm = 1;
    gm = dgm;
    int cnt = 0;
    while(std::pow(10,sigdig)*dgm >= gm && cnt < cmax)
    {
      dgm *= x/(1.5+cnt);
      gm += dgm;
      cnt += 1;
    }
    if(cnt < cmax) { gamma = 2*std::pow(x,0.5)*std::exp(-x)*gm; }
    else { gamma = NA_REAL; }
  }
  return gamma;
}

double GammaBigOneHalf(double x)
{
  double gamma;
  double dgm;
  double gm;
  int sigdig = 13;
  int cmax = 2000;
  if(x <= 0) { gamma = 1.77245385090552; }
  else if(x > 708.396418532264) { gamma = 0.0; }
  else if(x < 27)
  {
    dgm = 1;
    gm = dgm;
    int cnt = 0;
    while(std::pow(10,sigdig)*dgm >= gm && cnt < cmax)
    {
      dgm *= x/(1.5+cnt);
      gm += dgm;
      cnt += 1;
    }
    if(cnt < cmax) { gamma = 1.77245385090552-2*std::pow(x,0.5)*std::exp(-x)*gm; }
    else { gamma = NA_REAL; }
  }
  else
  {
    int casym = static_cast<int>(std::round(x));
    dgm = 1;
    gm = 1;
    int cnt = 0;
    while(cnt < casym)
    {
      cnt += 1;
      dgm *= (0.5-cnt)/x;
      gm += dgm;
    }
    gamma = std::pow(x,-0.5)*std::exp(-x)*gm;
  }
  return gamma;
}

double PoincaireExpansion(double a)
{
  Rcpp::NumericVector B2k(16);
  // absolute value of Bernoulli numbers (they alternate in sign)
  B2k[0] = 0.166666666666667;
  B2k[1] = 3.33333333333333E-02;
  B2k[2] = 2.38095238095238E-02;
  B2k[3] = 3.33333333333333E-02;
  B2k[4] = 7.57575757575758E-02;
  B2k[5] = 0.253113553113553;
  B2k[6] = 1.16666666666667;
  B2k[7] = 7.0921568627451;
  B2k[8] = 54.9711779448622;
  B2k[9] = 529.124242424242;
  B2k[10] = 6192.1231884058;
  B2k[11] = 86580.2531135531;
  B2k[12] = 1425517.16666667;
  B2k[13] = 27298231.0678161;
  B2k[14] = 601580873.900642;
  B2k[15] = 15116315767.0922;
  double gammaln = (a-0.5)*std::log(a)-a+0.918938533204673;
  int k = 0;
  while(k < 15)
  {
    k += 1;
    double dgm = std::log(B2k[k-1])-std::log(2*k*(2*k-1))-(2*k-1)*std::log(a);
    gammaln += std::exp(dgm);
    k += 1;
    dgm = std::log(B2k[k-1])-std::log(2*k*(2*k-1))-(2*k-1)*std::log(a);
    gammaln -= std::exp(dgm);
  }
  return(gammaln);
}

double GammaComplete(double a)
{
  double gamma;
  // small to medium a
  if(a < 50)
  {
    // split a into integer portion, p, and remainder, r
    double p = std::floor(a);
    double r = a-p;
    // check for a equal to zero or negative integer
    if(r == 0 && p < 1) { gamma = NA_REAL; }
    else
    {
      // special cases
      if(2*a-std::floor(2*a) == 0)
      {
        if(r == 0)
        {
          gamma = 1.0;
          double j = 0.0;
          while(p-1-j > 1e-15)
          {
            j += 1;
            gamma *= j;
          }
        }
        else
        {
          gamma = 1.77245385090552;
          if(a > 0)
          {
            double j = -0.5;
            while(p-0.5-j > 1e-15)
            {
              j += 1;
              gamma *= j;
            }
          }
          else
          {
            double j = 0.5;
            while(j-(p+0.5) > 1e-15)
            {
              j -= 1;
              gamma /= j;
            }
          }
        }
      }
      // other cases
      else
      {
        int sigdig = 15;
        int cmax = 2000;
        // confluent hypergeometric in well-behaved region 1 < r < 2
        double dgm = 1.0;
        double gm = dgm;
        int cnt = 0;
        while(cnt < cmax && std::pow(10,sigdig)*dgm >= gm)
        {
          cnt += 1;
          dgm *= 199/(r+1+cnt);
          gm += dgm;
        }
        gm *= std::pow(199,r+1)/(r+1)*std::exp(-199);
        // factorial calculations (gm is gamma of r+1)
        double j = r;
        while(p-1+r-j > 1e-15)
        {
          j += 1;
          gm *= j;
        }
        j = 1+r;
        while(j-(p+r) > 1e-15)
        {
          j -= 1;
          gm /= j;
        }
        gamma = gm;
      }
    }
  }
  // large a
  else { gamma = std::exp(PoincaireExpansion(a)); }

  return(gamma);
}

double GammaLn(double a)
{
  double gammaln;
  // positive a
  if(a <= 0) { gammaln = NA_REAL; }
  // ln of GammaComplete for small to medium a
  else if(a < 50) { gammaln = std::log(GammaComplete(a)); }
  // large a
  else { gammaln = PoincaireExpansion(a); }

  return(gammaln);
}

double GammaSmallRatio(double a, double x)
{
  double gammaratio;
  // zero or negative a
  if(a <= 0) { gammaratio = NA_REAL; }
  // degenerate solution
  else if(x <= 0) { gammaratio = 0.0; }
  // series expansion using logs
  else
  {
    int sigdig = 15;
    int cmax = 30000;
    double lnx = std::log(x);
    double sumlnxlna = 0.0;
    double dgm = 1.0;
    double gm = dgm;
    int cnt = 1;
    while(std::pow(10,sigdig)*dgm >= gm && cnt <= cmax)
    {
      sumlnxlna += lnx-std::log(cnt+a);
      dgm = std::exp(sumlnxlna);
      gm += dgm;
      cnt += 1;
    }
    double lngm = -std::log(a)+a*lnx-x+std::log(gm)-GammaLn(a);
    if(lngm <= -27.6310211159285) { gammaratio = 0.0; }
    else if(lngm >= -1e-12) { gammaratio = 1.0; }
    else { gammaratio = std::exp(lngm); }
  }
  return(gammaratio);
}

double GammaBigRatio(double a, double x)
{
  double gammaratio;
  // zero or negative a
  if(a <= 0) { gammaratio = NA_REAL; }
  // degenerate solution
  else if(x <= 0) { gammaratio = 1.0; }
  // by subtraction
  else { gammaratio = 1-GammaSmallRatio(a, x); }

  return(gammaratio);
}
