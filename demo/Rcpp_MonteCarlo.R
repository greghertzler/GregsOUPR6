# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Rcpp_MonteCarlo.R",echo=TRUE)', or
# Type 'demo(Rcpp_MonteCarlo)'
# Rcpp functions are exported to console
# m=100,skip=1,paths=1000000,seed=99
stdnorm <- RcppOUPMCStandardNormal(100,1,1000000,99)
# x=20,m=100,skip=1,dt=0.1,rho=0.5,mu=-15,sigma=15
forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,20,100,1,0.1,0.5,-15,15)
y <- seq(from=-30,to=30,by=0.6)
# bin into y and sum from left to right with psi=-1
mvdpd <- RcppOUPMCForwardCountY(forward,y,-1)
# subset
means <- mvdpd[,1,drop=FALSE]
variances <- mvdpd[,2,drop=FALSE]
densities <- mvdpd[,3:(100+2),drop=FALSE]
probabilities <- mvdpd[,(100+3):(2*100+2),drop=FALSE]
doubleintegrals <- mvdpd[,(2*100+3):(3*100+2),drop=FALSE]
# print means
means
