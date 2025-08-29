# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_RateOfConvergence.R",echo=TRUE)', or
# Type 'demo(Estimation_RateOfConvergence)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
df<-read.csv("data/OUP_Convergence.csv")
# rho = 0.1
ML$Estimates()
ML$Estimates(df=df)
# rho = 0.5
ML$Estimates(df=df,taucol=1,zcol=3)
# rho = 2.5
ML$Estimates(df=df,taucol=1,zcol=4)
