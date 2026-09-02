# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Rcpp_MaximumLikelihood.R",echo=TRUE)', or
# Type 'demo(Rcpp_MaximumLikelihood)'
# read data set
df <- OUPDataRead("OUP_Convergence.csv")
# allocate matrix for results
n <- ncol(df)
estimates <- matrix("double",n-1,3)
# find starting values for first estimation
tau <- df[[1]]
z <- df[[2]]
# Rcpp functions are exported to console
thetasteps <- RcppOUPMLNMStart(tau,z)
# first estimation, data must be 'clean' or this will fail
estimates[1,] <- RcppOUPMLNelderMead(tau,z,thetasteps[1,1:3],thetasteps[2,1:3])
# iterate through data set, previous estimates become starting values
i <- 2
while(i < n)
{
  z <- df[[i+1]]
  estimates[i,] <- RcppOUPMLNelderMead(tau,z,thetasteps[1,1:3],thetasteps[2,1:3])
  i <- i+1
}
estimates
