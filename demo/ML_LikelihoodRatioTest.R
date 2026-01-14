# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_LikelihoodRatioTest.R",echo=TRUE)', or
# Type 'demo(ML_LikelihoodRatioTest)'
# R6 object
ML <- MaximumLikelihood$new()
# default data comparing unrestricted and true parameters
ML$LikelihoodRatioTest()
# other data comparing unrestricted and true parameters
df<-OUPReadData("OUP_NotMissing")
ML$Estimates(df=df,plotit=FALSE)
ML$Estimates(rhor=0.5,mur=-15,sigmar=15,plotit=FALSE)
ML$LikelihoodRatioTest()
# or doing it the more transparent way
u <- ML$Estimates(df=df,plotit=FALSE)
r <- ML$Estimates(rhor=0.5,mur=-15,sigmar=15,plotit=FALSE)
ML$LikelihoodRatioTest(u$lnLu,r$lnLr,r$alphar,r$m1)
# other columns in data
ML$Estimates(df=df,taucol=1,zcol=3,plotit=FALSE)
ML$Estimates(rhor=10,mur=-15,sigmar=67.08204,plotit=FALSE)
ML$LikelihoodRatioTest()
# 95% lower bound on sigma
ML$Estimates(sigmar=35.00838,plotit=FALSE)
ML$set_timeseries_info(estimation="P(sigma>35.01)=95%")
ML$LikelihoodRatioTest()
# with plots this time
ML$PlotTimeSeries(df=df,taucol=1,zcol=3)
ML$Estimates()
ML$Estimates(sigmar=35.00838)
ML$set_timeseries_info(estimation="P(sigma>35.01)=95%")
ML$LikelihoodRatioTest()
# customize plot
ML$PlotEstimates(tbeg=4,tend=24)
