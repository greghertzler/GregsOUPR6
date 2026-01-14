# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_LogLikelihood.R",echo=TRUE)', or
# Type 'demo(ML_LogLikelihood)'
# R6 object
ML <- MaximumLikelihood$new()
# default data and true parameters
ML$LogLikelihood()
# other parameters
ML$LogLikelihood(rho=0.4,mu=-10,sigma=25)
# other data
df<-OUPReadData("OUP_NotMissing")
ML$LogLikelihood(rho=0.5,mu=-15,sigma=15,df=df)
# other columns in data
ML$LogLikelihood(df=df,taucol=1,zcol=3)
# try again
ML$LogLikelihood(rho=10,sigma=67)
# with plots this time
ML$PlotTimeSeries(df=df)
ML$LogLikelihood(rho=0.5,mu=-15,sigma=15)
ML$PlotEstimates()
# customize plot
ML$PlotEstimates(title="My Log Likelihood",tbeg=100,tend=150)
