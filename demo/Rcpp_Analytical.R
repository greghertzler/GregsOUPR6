# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Rcpp_Analytical.R",echo=TRUE)', or
# Type 'demo(Rcpp_Analytical)'
# calculate option prices
s <- seq(from=30,to=10,by=-0.1)
x <- seq(from=-40,to=60,by=0.5)
# Rcpp functions are exported to console
OO <- RcppOUPAOption(s,x,30,0,0.5,15,15,0.05,1,0,0)
# print option prices
OO
