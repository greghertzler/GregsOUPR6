# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_GoodnessOfFit.R",echo=TRUE)', or
# Type 'demo(ML_GoodnessOfFit)'
# R6 object
ML <- MaximumLikelihood$new()
# default data and true parameters
ML$GoodnessOfFit()
# bad parameters
ML$GoodnessOfFit(rho=0.4,mu=-10,sigma=25)
# unrestricted estimates
ML$Estimates(plotit=FALSE)
ML$GoodnessOfFit()
# invariant stochastic process
ML$Estimates(rhor=10000,plotit=FALSE)
ML$GoodnessOfFit()
# scaled brownian motion
ML$Estimates(rhor=0,plotit=FALSE)
ML$GoodnessOfFit()
# other data
df<-myDataRead("OUP_NotMissing")
ML$Estimates(df=df,plotit=FALSE)
ML$GoodnessOfFit()
# other columns in data
ML$Estimates(df=df,taucol=1,zcol=3,plotit=FALSE)
ML$GoodnessOfFit()
# with plots this time
ML$PlotTimeSeries(df=df)
ML$Estimates()
ML$GoodnessOfFit()
# customize plot
ML$PlotEstimates(title="My Goodness of Fit",tbeg=50,tend=100)
