# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_Estimates.R",echo=TRUE)', or
# Type 'demo(ML_Estimates)'
# R6 object
ML <- MaximumLikelihood$new()
# default data and unrestricted estimates
ML$Estimates()
# restricted estimates
ML$Estimates(rhor=0.5,mur=-15)
# other data
df<-OUPDataRead("OUP_NotMissing")
ML$Estimates(df=df)
# other columns in data
ML$Estimates(df=df,taucol=1,zcol=3)
# with automatic plot this time
ML$set_flags(plotit=TRUE)
ML$Estimates(df=df,taucol=1,zcol=3)
# customize plot
ML$PlotEstimates(title="My Estimates",tbeg=200)
