# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Rcpp_FiniteDifference.R",echo=TRUE)', or
# Type 'demo(Rcpp_FiniteDifference)'
# calculate option prices
s <- seq(from=30,to=10,by=-0.1)
x <- seq(from=-40,to=60,by=0.5)
n <- length(x)
V <- vector("double",n)
for(j in 1:n) { V[j] <- max(0,x[j]) }
# Rcpp functions are exported to console
OO <- RcppOUPFDOption(s,x,V,0.05,0.5,10,0.5,15,15)
# print option prices
OO
