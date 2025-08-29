# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_ObservationInterval.R",echo=TRUE)', or
# Type 'demo(Estimation_ObservationInterval)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
df<-read.csv("data/OUP_ObservationInterval.csv")
# plot by day
ML$PlotTimeSeries(df,taucol=1,zcol=2)
# plot by year
ML$PlotTimeSeries(df,taucol=3,zcol=2)
# estimate by day
ML$Estimates(df,taucol=1,zcol=2)
# estimate by year
ML$Estimates(df=df,taucol=3,zcol=2)
