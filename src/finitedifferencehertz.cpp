#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <limits>

// roxygen (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @title FiniteDifference_Rcpp functions for numerical solutions of an
//'  Ornstein-Uhlenbeck Process
//'
//' @description
//' Calculations for the R6 class 'FiniteDifference', with sequential processing.
//'
//' @details # Notes on Values
//' Return values are vectors and matrices allocated in Rcpp.  The dimensions are
//'  shown for information.  Of course, do not include them in R calls.  For example:
//'
//'     g <- RcppOUPFDDrift(x,rho,mu)
//'
//' The return values:
//'
//'     g(n)
//'     h2(n)
//'     V(n)
//'
//'  are vectors of the same dimension as x.
//'
//' The return value:
//'
//'     c(m,n)
//'
//'  is a matrix with row dimension from s and column dimension from x.
//'
//' The return value:
//'
//'     env(2,n)
//'
//'  is a matrix with two row vectors for the option prices along the envelope
//'  and the corresponding times.  It is subset in R as:
//'
//'     OOenv <- env[1,]
//'     tsenv <- env[2,]
//'
//' The return value:
//'
//'     dec(2)
//'
//'  is a vector with two elements for the state at the decision threshold and
//'  the option price at that state.  It is subset in R as:
//'
//'     k <- dec[1]
//'     OOhat <- dec[2]
//'
//' @details # Discussion
//' First, the finite difference method was implemented in R6 as a
//'  single-threaded application.  Then the R6 code was translated into Rcpp
//'  sequential code.  The finite difference method steps backward in time.  For
//'  each time, it uses LU decomposition, followed by forward and backward
//'  substitution to solve matrices.  The finite difference method is inherently
//'  sequential.  Global searches for envelopes of option prices and decision
//'  thresholds are also sequential.  There is no parallel processing in this
//'  module.  The R6 single-threaded code has been archived, leaving the Rcpp code.
//'
//' Below are microbenchmark median times to calculate 400,000 Option Prices.
//'  Parameter skip=10 and only 40,000 Option Prices are reported.  Other rows
//'  show median times for the Option Envelope and the Decision Threshold.
//'  Calculations are on an i7 CPU with a maximum speed of 4.5 GHz.
//'
//'     Unit: milliseconds             R6     R6+
//'              function   single-thread    Rcpp
//'     -----------------------------------------
//'                Option        342.1143  7.7743
//'        OptionEnvelope        608.9655  4.6913
//'     DecisionThreshold          6.0607  0.1089
//'
//' The Option Envelope is an input into the calculation of the Decision Threshold.
//'  If the Decision Threshold is calculated from a standing start, add the two times.
//'  R6+Rcpp calculates from 44.0 to 129.8 times faster than R6 single-thread.
//'  The algorithms are not identical.  R6+Rcpp takes advantage of void subroutines
//'  and call-by-reference to manipulate vectors and matrices without copying.
//'
//' @details # From the Console
//' These functions are available in R, the RStudio console and RShiny apps.  For
//'  example, a  calculation of option prices would be:
//'
//'     s <- seq(from=30,to=10,by=-0.1)
//'     x <- seq(from=-40,to=60,by=0.5)
//'     n <- length(x)
//'     V <- vector("double",n)
//'     for(j in 1:n) { V[j] <- max(0,x[j]) }
//'     OO <- RcppOUPFDOption(s,x,V,0.05,0.5,10,0.5,15,15)
//'
//' Calling functions from the console is slightly faster than calling them
//'  through R6 objects.  The table below compares the R6+Rcpp median times
//'  above with the median times from the console:
//'
//'     Unit: milliseconds     R6+  Console
//'              function     Rcpp     Rcpp
//'     -----------------------------------
//'                Option   7.7743   7.3196
//'        OptionEnvelope   4.6913   4.5448
//'     DecisionThreshold   0.1089   0.0174
//'
//' There is less than half a millisecond time penalty for using the R6 object
//'  and it is much more convenient.  It is reactive.  In other words, it
//'  stores inputs and outputs and maps inputs to outputs.  Changing an input will
//'  nullify dependent outputs, eliminating any danger of reporting a stale output.
//'  Outputs are calculated only as needed and only once. Then they are reused.
//'  In plots using Plotly, for example.
//'
//' Potentially, the Rcpp functions could be imported into other packages.
//'
//' @name FiniteDifference_Rcpp

// Helpers (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

void OptionA(int n, NumericMatrix& A, NumericVector g, NumericVector h2, double r, double theta, double ds, double dx)
{
  // Rcout << "OptionA" << std::endl;
  // Rcout << n << ", " << r << ", " << theta << ", " << ds << ", " << dx << std::endl;
  double dx2 = dx*dx;
  // first row has 3 entries
  A(0,0) = 1/ds+theta*r+0.5*theta*(3*g[0]/dx-h2[0]/dx2);
  A(0,1) = -theta*(2*g[0]/dx-h2[0]/dx2);
  A(0,2) = 0.5*theta*(g[0]/dx-h2[0]/dx2);
  // Rcout << 0 << ": " << A(0,0) << ", " << A(0,1) << ", " << A(0,2) << std::endl;
  // middle rows are tridiagonal
  for(int j = 1; j < n-1; j++)
  {
    A(j,j-1) = 0.5*theta*(g[j]/dx-h2[j]/dx2);
    A(j,j) = 1/ds+theta*r+theta*h2[j]/dx2;
    A(j,j+1) = -0.5*theta*(g[j]/dx+h2[j]/dx2);
    // Rcout << j << ": " << A(j,j-1) << ", " << A(j,j) << ", " << A(j,j+1) << std::endl;
  }
  // last row has 3 entries
  A(n-1,n-3) = -0.5*theta*(g[n-1]/dx+h2[n-1]/dx2);
  A(n-1,n-2) = theta*(2*g[n-1]/dx+h2[n-1]/dx2);
  A(n-1,n-1) = 1/ds+theta*r-0.5*theta*(3*g[n-1]/dx+h2[n-1]/dx2);
  // Rcout << n-1 << ": " << A(n-1,n-3) << ", " << A(n-1,n-2) << ", " << A(n-1,n-1) << std::endl;
}

void OptionLU(int n, NumericMatrix& A)
{
  // LU is in A with implicit 1 on diagonal of L
  // begin in second row and stop before last row is hit
  for(int j = 0; j < n-3; j++)
  {
    int i = j+1;
    A(i,j) /= A(j,j);
    for(int k = j+1; k < j+3; k++) { A(i,k) -= A(i,j)*A(j,k); }
  }
  // last 2 rows done this way because last row has 3 entries
  for(int j = n-3; j < n-1; j++)
  {
    for(int i = j+1; i < n; i++)
    {
      A(i,j) /= A(j,j);
      for(int k = j+1; k < n; k++) { A(i,k) -= A(i,j)*A(j,k); }
    }
  }
}

void Optionb(int i, int n, NumericVector& b, NumericMatrix c, NumericVector g, NumericVector h2, double r, double theta, double ds, double dx)
{
  // new b for each iteration
  double dx2 = dx*dx;
  b[0] = (1/ds-(1-theta)*r)*c(i-1,0)+0.5*(1-theta)*(g[0]*(-3*c(i-1,0)+4*c(i-1,1)-c(i-1,2))/dx+h2[0]*(c(i-1,0)-2*c(i-1,1)+c(i-1,2))/dx2);
  for(int j = 1; j < n-1; j++) { b[j] = (1/ds-(1-theta)*r)*c(i-1,j)+0.5*(1-theta)*(g[j]*(-c(i-1,j-1)+c(i-1,j+1))/dx+h2[j]*(c(i-1,j-1)-2*c(i-1,j)+c(i-1,j+1))/dx2); }
  b[n-1] = (1/ds-(1-theta)*r)*c(i-1,n-1)+0.5*(1-theta)*(g[n-1]*(c(i-1,n-3)-4*c(i-1,n-2)+3*c(i-1,n-1))/dx+h2[n-1]*(c(i-1,n-3)-2*c(i-1,n-2)+c(i-1,n-1))/dx2);
}

void OptionSolve(int i, int n, NumericMatrix A, NumericVector& b, NumericMatrix& c)
{
  // forward substitution to find q as solution of Lq=b
  for(int j = 1; j < n-1; j++) { b[j] -= A(j,j-1)*b[j-1]; }
  b[n-1] = b[n-1]-A(n-1,n-3)*b[n-3]-A(n-1,n-2)*b[n-2];
  // backward substitution to find c as solution of Uc=q
  c(i,n-1) = b[n-1]/A(n-1,n-1);
  for(int j = n-2; j > 0; j--) { c(i,j) = (b[j]-A(j,j+1)*c(i,j+1))/A(j,j); }
  c(i,0) = (b[0]-A(0,1)*c(i,1)-A(0,2)*c(i,2))/A(0,0);
}

double Lagrange3Point(double OO0, double OO1, double OO2, double x0, double x1, double x2, double xx)
{
  double P0 = ((xx-x1)*(xx-x2))/((x0-x1)*(x0-x2));
  double P1 = ((xx-x0)*(xx-x2))/((x1-x0)*(x1-x2));
  double P2 = ((xx-x0)*(xx-x1))/((x2-x0)*(x2-x1));

  return OO0*P0+OO1*P1+OO2*P2;
}

double Lagrange7Point(NumericVector OO, NumericVector x, double xx)
{
  NumericVector P(7);
  P[0] = ((xx-x[1])*(xx-x[2])*(xx-x[3])*(xx-x[4])*(xx-x[5])*(xx-x[6]))/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3])*(x[0]-x[4])*(x[0]-x[5])*(x[0]-x[6]));
  P[1] = ((xx-x[0])*(xx-x[2])*(xx-x[3])*(xx-x[4])*(xx-x[5])*(xx-x[6]))/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3])*(x[1]-x[4])*(x[1]-x[5])*(x[1]-x[6]));
  P[2] = ((xx-x[0])*(xx-x[1])*(xx-x[3])*(xx-x[4])*(xx-x[5])*(xx-x[6]))/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3])*(x[2]-x[4])*(x[2]-x[5])*(x[2]-x[6]));
  P[3] = ((xx-x[0])*(xx-x[1])*(xx-x[2])*(xx-x[4])*(xx-x[5])*(xx-x[6]))/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2])*(x[3]-x[4])*(x[3]-x[5])*(x[3]-x[6]));
  P[4] = ((xx-x[0])*(xx-x[1])*(xx-x[2])*(xx-x[3])*(xx-x[5])*(xx-x[6]))/((x[4]-x[0])*(x[4]-x[1])*(x[4]-x[2])*(x[4]-x[3])*(x[4]-x[5])*(x[4]-x[6]));
  P[5] = ((xx-x[0])*(xx-x[1])*(xx-x[2])*(xx-x[3])*(xx-x[4])*(xx-x[6]))/((x[5]-x[0])*(x[5]-x[1])*(x[5]-x[2])*(x[5]-x[3])*(x[5]-x[4])*(x[5]-x[6]));
  P[6] = ((xx-x[0])*(xx-x[1])*(xx-x[2])*(xx-x[3])*(xx-x[4])*(xx-x[5]))/((x[6]-x[0])*(x[6]-x[1])*(x[6]-x[2])*(x[6]-x[3])*(x[6]-x[4])*(x[6]-x[5]));

  return OO[0]*P[0]+OO[1]*P[1]+OO[2]*P[2]+OO[3]*P[3]+OO[4]*P[4]+OO[5]*P[5]+OO[6]*P[6];
}

// Exports (((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDDrift(x,rho,mu)
//' @param  x   vector of states
//' @param  rho rate parameter 0<=rho<inf
//' @param  mu  location parameter -inf<mu<inf
//' @return g(n) <- RcppOUPFDDrift()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDDrift(NumericVector x, double rho, double mu)
{
  int n = x.size();
  NumericVector g(n);
  for(int j = 0; j < n; j++) { g[j] = -rho*(x[j]-mu); }

  return g;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDDiffusion(x,sigma)
//' @param  x     vector of states
//' @param  sigma scale parameter -inf<sigma<inf
//' @return h2(n) <- RcppOUPFDDiffusion()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDDiffusion(NumericVector x, double sigma)
{
  int n = x.size();
  NumericVector h2(n);
  for(int j = 0; j < n; j++) { h2[j] = sigma*sigma; }

  return h2;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Linear(x,xo,vs)
//' @param  x  vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  vs slope
//' @return V(n) <- RcppOUPFDTerminalValue_Linear()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Linear(NumericVector x, double xo, double vs)
{
  int n = x.size();
  NumericVector V(n);
  for(int j = 0; j < n; j++) { V[j] = vs*(x[j]-xo); }

  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Degenerate(x,xo,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Degenerate()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Degenerate(NumericVector x, double xo, double Vmax, double Vmin)
{
  int n = x.size();
  double dx = (x[n-1]-x[0])/(n-1);
  NumericVector V(n);
  for(int j = 0; j < n; j++)
  {
    if(x[j] >= xo-0.5*dx && x[j] < xo+0.5*dx) { V[j] = Vmax; }
    else { V[j] = Vmin; }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Stepped(x,xo,vs,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  vs   slope
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Stepped()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Stepped(NumericVector x, double xo, double vs, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  if(vs > 0)
  {
    for(int j = 0; j < n; j++)
    {
      if(x[j] < xo) { V[j] = Vmin; }
      else if(x[j] == xo) { V[j] = 0.5*(Vmax+Vmin); }
      else { V[j] = Vmax; }
    }
  }
  else
  {
    for(int j = 0; j < n; j++)
    {
      if(x[j] < xo) { V[j] = Vmax; }
      else if(x[j] == xo) { V[j] = 0.5*(Vmax+Vmin); }
      else { V[j] = Vmin; }
    }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Kinked(x,xo,vs,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  vs   slope
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Kinked()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Kinked(NumericVector x, double xo, double vs, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  for(int j = 0; j < n; j++)
  {
    V[j] = vs*(x[j]-xo);
    if(V[j] > Vmax) { V[j] = Vmax; }
    else if(V[j] < Vmin) { V[j] = Vmin; }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Butterfly(x,xo,xm,vs,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  xm   state at the maximum or right wing
//' @param  vs   slope
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Butterfly()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Butterfly(NumericVector x, double xo, double xm, double vs, double Vmax, double Vmin)
{
  int n = x.size();
  double slope = abs(vs);
  NumericVector V(n);
  for(int j = 0; j < n; j++)
  {
    double left = -slope*(x[j]-xo);
    double right = slope*(x[j]-xm);
    if(left > right) { V[j] = left; }
    else { V[j] = right; }
    if(V[j] > Vmax) { V[j] = Vmax; }
    else if(V[j] < Vmin) { V[j] = Vmin; }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Mitscherlich(x,xo,vr,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  vr   rate of change
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Mitscherlich()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Mitscherlich(NumericVector x, double xo, double vr, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  for(int j = 0; j < n; j++)
  {
    if(vr*(x[j]-xo) > 0) { V[j] = Vmin+(Vmax-Vmin)*(1-exp(-vr*(x[j]-xo))); }
    else { V[j] = Vmin; }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Gompertz(x,xi,vr,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xi   state at the inflection point
//' @param  vr   rate of change
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Gompertz()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Gompertz(NumericVector x, double xi, double vr, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  if(vr > 0)
  {
    for(int j = 0; j < n; j++) { V[j] = Vmin+(Vmax-Vmin)*exp(-exp(-vr*(x[j]-xi))); }
  }
  else
  {
    for(int j = 0; j < n; j++) { V[j] = Vmax-(Vmax-Vmin)*exp(-exp(vr*(x[j]-xi))); }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Logistic(x,xi,vr,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xi   state at the inflection point
//' @param  vr   rate of change
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Logistic()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Logistic(NumericVector x, double xi, double vr, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  if(vr > 0)
  {
    for(int j = 0; j < n; j++) { V[j] = Vmin+(Vmax-Vmin)/(1+exp(-vr*(x[j]-xi))); }
  }
  else
  {
    for(int j = 0; j < n; j++) { V[j] = Vmax-(Vmax-Vmin)/(1+exp(vr*(x[j]-xi))); }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_Transcendental(x,xo,xi,xm,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  xi   state at the inflection point
//' @param  xm   state at the maximum or right wing
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_Transcendental()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_Transcendental(NumericVector x, double xo, double xi, double xm, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  if(abs(xm-xi) < 0.0000000001 || Vmax-Vmin < 0.0000000001 || (xm-xi)*(xi-xo) < 0)
  {
    for(int j = 0; j < n; j++) { V[j] = Vmin; }
  }
  else
  {
    double b = pow(((xm-xo)/(xm-xi)),2);
    double c = (xm-xo)/pow((xm-xi),2);
    for(int j = 0; j < n; j++)
    {
      if((x[j]-xo)*(xm-xo) > 0)
      {
        double lnv = log(Vmax-Vmin)+b*log((x[j]-xo)/(xm-xo))-c*(x[j]-xm);
        V[j] = Vmin+exp(lnv);
      }
      else { V[j] = Vmin; }
    }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDTerminalValue_YieldIndex(x,xo,xi,xm,Vmax,Vmin)
//' @param  x    vector of states
//' @param  xo   state at the intercept, step, kink or left wing
//' @param  xi   state at the inflection point
//' @param  xm   state at the maximum or right wing
//' @param  Vmax maximum terminal value
//' @param  Vmin minimum terminal value
//' @return V(n) <- RcppOUPFDTerminalValue_YieldIndex()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDTerminalValue_YieldIndex(NumericVector x, double xo, double xi, double xm, double Vmax, double Vmin)
{
  int n = x.size();
  NumericVector V(n);
  if(abs(xm-xi) < 0.0000000001 || Vmax-Vmin < 0.0000000001 || (xm-xi)*(xi-xo) < 0)
  {
    for(int j = 0; j < n; j++) { V[j] = Vmax; }
  }
  else
  {
    double b = pow(((xm-xo)/(xm-xi)),2);
    double c = (xm-xo)/pow((xm-xi),2);
    for(int j = 0; j < n; j++)
    {
      if((x[j]-xo)*(xm-xo) > 0)
      {
        double lnv = log(Vmax-Vmin)+b*log((x[j]-xo)/(xm-xo))-c*(x[j]-xm);
        V[j] = Vmax-exp(lnv);
        if(V[j] < 0) { V[j] = 0; }
      }
      else { V[j] = Vmax; }
    }
  }
  return V;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDOption(s,x,V,r,theta,skip,rho,mu,sigma)
//' @param  s     vector of times
//' @param  x     vector of states
//' @param  V     vector of terminal values
//' @param  r     discount rate 0<=r<inf
//' @param  theta weight of current time in time stepping 0.5<=theta<=1
//' @param  skip  subdivide time intervals but report at times s 1<=skip<=1000
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @return c(m,n) <- RcppOUPFDOption()
//' @export
// [[Rcpp::export]]
NumericMatrix RcppOUPFDOption(NumericVector s, NumericVector x, NumericVector V, double r, double theta, int skip, double rho, double mu, double sigma)
{
  Rcout << "Option" << std::endl;
  int m = s.size();
  int n = x.size();
  Rcout << m << ":" << n << ", " << s[0] << ", " << s[m-1] << ", " << x[0] << ", " << x[n-1] << std::endl;
  NumericMatrix A(n,n);
  NumericMatrix c(m,n);
  NumericMatrix cskip(skip,n);
  NumericVector b(n);
  NumericVector g(n);
  NumericVector h2(n);
  for(int j = 0; j < n; j++)
  {
    c(0,j) = V[j];
    g[j] = -rho*(x[j]-mu);
    h2[j] = sigma*sigma;
  }
  double ds = abs(s[0]-s[m-1])/(m-1);
  double dsskip = ds/skip;
  double dx = abs(x[n-1]-x[0])/(n-1);
  Rcout << ds << ", " << dsskip << ", " << dx << std::endl;
  OptionA(n,A,g,h2,r,theta,dsskip,dx);
  OptionLU(n,A);
  for(int i = 1; i < m; i++)
  {
    Optionb(i,n,b,c,g,h2,r,theta,dsskip,dx);
    for(int iskip = 1; iskip < skip; iskip++)
    {
      OptionSolve(iskip,n,A,b,cskip);
      Optionb(iskip+1,n,b,cskip,g,h2,r,theta,dsskip,dx);
    }
    OptionSolve(i,n,A,b,c);
  }
  Rcout << "Option " << (m-1)*skip+1 << ":" << n << std::endl;

  return c;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDOptionEnvelope(s,x,V,r,theta,skip,rho,mu,sigma)
//' @param  s     vector of times
//' @param  x     vector of states
//' @param  V     vector of terminal values
//' @param  r     discount rate 0<=r<inf
//' @param  theta weight of current time in time stepping 0.5<=theta<=1
//' @param  skip  subdivide time intervals but report at times s 1<=skip<=1000
//' @param  rho   rate parameter 0<=rho<inf
//' @param  mu    location parameter -inf<mu<inf
//' @param  sigma scale parameter -inf<sigma<inf
//' @return env(2,n) <- RcppOUPFDOptionEnvelope()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDOptionEnvelope(NumericVector s, NumericVector x, NumericVector V, double r, double theta, int skip, double rho, double mu, double sigma)
{
  int m = s.size();
  int n = x.size();
  int mskip = 500*skip+1;
  NumericMatrix A(n,n);
  NumericMatrix cskip(mskip,n);
  NumericVector b(n);
  NumericVector g(n);
  NumericVector h2(n);
  NumericMatrix env(2,n);
  LogicalVector up(n);
  bool hit = true;
  for(int j = 0; j < n; j++)
  {
    cskip(0,j) = V[j];
    env(0,j) = V[j];
    env(1,j) = 0;
    g[j] = -rho*(x[j]-mu);
    h2[j] = sigma*sigma;
    up[j] = false;
  }
  double ds = abs(s[0]-s[m-1])/(m-1);
  double dsskip = ds/skip;
  double dx = abs(x[n-1]-x[0])/(n-1);
  OptionA(n,A,g,h2,r,theta,dsskip,dx);
  OptionLU(n,A);
  int i;
  for(i = 1; i < mskip && hit; i++)
  {
    Optionb(i,n,b,cskip,g,h2,r,theta,dsskip,dx);
    OptionSolve(i,n,A,b,cskip);
    hit = false;
    for(int j = 0; j < n; j++)
    {
      if(cskip(i,j) > cskip(i-1,j))
      {
        hit = true;
        up[j] = true;
        if(cskip(i,j) > env(0,j))
        {
          env(0,j) = cskip(i,j);
          env(1,j) = i*dsskip;
        }
      }
      else if(cskip(i,j) <= cskip(i-1,j) && up[j] && i > 1)
      {
        up[j] = false;
        env(1,j) = (i-1+0.5*(cskip(i,j)-cskip(i-2,j))/(2*cskip(i-1,j)-cskip(i-2,j)-cskip(i,j)))*dsskip;
        env(0,j) = Lagrange3Point(cskip(i-2,j),cskip(i-1,j),cskip(i,j),(i-2)*dsskip,(i-1)*dsskip,i*dsskip,env(1,j));
        if(env(0,j) < V[j])
        {
          env(0,j) = V[j];
          env(1,j) = 0;
        }
      }
    }
  }
  Rcout << "Option Envelope " << i << ":" << n << std::endl;

  return env;
}

//' @rdname FiniteDifference_Rcpp
//' @usage  RcppOUPFDDecisionThreshold(x,V,OOenv,phi)
//' @param  x      vector of states
//' @param  V      vector of terminal values
//' @param  OOenv  envelope of option values
//' @param  phi    search direction for exit or entry options
//' @return dec(2) <- RcppOUPFDDecisionThreshold()
//' @export
// [[Rcpp::export]]
NumericVector RcppOUPFDDecisionThreshold(NumericVector x, NumericVector V, NumericVector OOenv, double phi)
{
  int n = x.size();
  NumericVector OO(7);
  NumericVector dec(2);
  double dx = (x[n-1]-x[0])/(n-1);
  double dxx = dx/100;
  // search for enter to right
  if(phi > 0 || (phi == 0 && V[n-1] > V[0]))
  {
    bool hit = false;
    dec[0] = x[n-1];
    dec[1] = OOenv[n-1];
    // moving right, find last point where envelope exceeds terminal value
    for(int j = 6; j < n-1 && hit == false; j++)
    {
      if(OOenv[j] > V[j] && OOenv[j+1] <= V[j+1])
      {
        hit = true;
        // extrapolate OOenv using 7 point Lagrange polynomial
        bool doublehit = false;
        double xx = x[j];
        double thisdiff = OOenv[j]-V[j];
        for(int k = 0; k < 200 && doublehit == false; k++)
        {
          xx += dxx;
          double OO = Lagrange7Point(OOenv[seq(j-6,j)],x[seq(j-6,j)],xx);
          // compare with interpolated V
          double VV = V[j]+(V[j+1]-V[j])/(x[j+1]-x[j])*(xx-x[j]);
          double prevdiff = thisdiff;
          thisdiff = OO-VV;
          if(thisdiff < 0 || thisdiff > prevdiff)
          {
            doublehit = true;
            dec[0] = xx;
            dec[1] = OO;
          }
        }
      }
    }
  }
  // search for exit to left
  else
  {
    bool hit = false;
    dec[0] = x[0];
    dec[1] = OOenv[0];
    // moving left, find last point where envelope exceeds terminal value
    for(int j = n-7; j > 0 && hit == false; j--)
    {
      if(OOenv[j] > V[j] && OOenv[j-1] <= V[j-1])
      {
        hit = true;
        // extrapolate OOenv using 7 point Lagrange polynomial
        bool doublehit = false;
        double xx = x[j];
        double thisdiff = OOenv[j]-V[j];
        for(int k = 0; k < 200 && doublehit == false; k++)
        {
          xx -= dxx;
          double OO = Lagrange7Point(OOenv[seq(j,j+6)],x[seq(j,j+6)],xx);
          // compare with interpolated V
          double VV = V[j]+(V[j-1]-V[j])/(x[j]-x[j-1])*(x[j]-xx);
          double prevdiff = thisdiff;
          thisdiff = OO-VV;
          if(thisdiff < 0 || thisdiff > prevdiff)
          {
            doublehit = true;
            dec[0] = xx;
            dec[1] = OO;
          }
        }
      }
    }
  }
  return dec;
}
