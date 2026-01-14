# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_SmoothedData.R",echo=TRUE)', or
# Type 'demo(Estimation_SmoothedData)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
df<-OUPReadData("OUP_SmoothedData")
# estimate for raw data
ML$Estimates(df,taucol=1,zcol=2)
# estimate for one smoothing
ML$Estimates(df=df,taucol=1,zcol=3)
# estimate for two smoothings
ML$Estimates(df,taucol=1,zcol=4)
# estimate for three smoothings
ML$Estimates(df=df,taucol=1,zcol=5)
# estimate for four smoothings
ML$Estimates(df,taucol=1,zcol=6)
# estimate for five smoothings
ML$Estimates(df=df,taucol=1,zcol=7)
# estimate for six smoothings
ML$Estimates(df,taucol=1,zcol=8)
# estimate for seven smoothings
ML$Estimates(df=df,taucol=1,zcol=9)
# estimate for eight smoothings
ML$Estimates(df=df,taucol=1,zcol=10)
# estimate for nine smoothings
ML$Estimates(df,taucol=1,zcol=11)
