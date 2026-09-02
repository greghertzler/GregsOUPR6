#ifndef GAMMA_H
#define GAMMA_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

double GammaSmallOneHalf(double x);

double GammaBigOneHalf(double x);

double PoincaireExpansion(double a);

double GammaComplete(double a);

double GammaSmallRatio(double a, double x);

double GammaBigRatio(double a, double x);

#endif // GAMMA_H
