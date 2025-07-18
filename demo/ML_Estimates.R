# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_Estimates.R",echo=TRUE)', or
# Type 'demo(ML_Estimates)'
# R6 object
ML <- MaximumLikelihood$new()
# default data and unrestricted estimates
ML$Estimates(plotit=FALSE)
# restricted estimates
ML$Estimates(rhor=0.5,mur=-15,plotit=FALSE)
# other data
df<-read.csv("data/ML_OUP_UnequalIntervals.csv")
ML$Estimates(df=df,plotit=FALSE)
# other columns in data
ML$Estimates(df=df,taucol=3,zcol=4,plotit=FALSE)
# with plot this time
ML$Estimates(df=df,taucol=3,zcol=4)
# customize plot
ML$PlotEstimates(title="My Estimates",tbeg=4,tend=6)
