# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_NoMissingObservations.R",echo=TRUE)', or
# Type 'demo(Estimation_NoMissingObservations)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
df<-read.csv("data/OUP_NotMissing.csv")
# estimate for equal observation intervals
ML$Estimates(df,taucol=1,zcol=2)
# estimate for unequal observation intervals
ML$Estimates(df=df,taucol=1,zcol=3)
# estimate with averages to make equal intervals
ML$Estimates(df,taucol=1,zcol=4)
# estimate with means to make equal intervals
ML$Estimates(df=df,taucol=1,zcol=5)
