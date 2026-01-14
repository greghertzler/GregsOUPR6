# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_SampleSize.R",echo=TRUE)', or
# Type 'demo(Estimation_SampleSize)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
filePath <- paste0(myDataPath(),"OUP_SampleSize.csv")
df<-read.csv(filePath)
# estimate for small
ML$Estimates(df,taucol=1,zcol=2)
# estimate for medium
ML$Estimates(df=df,taucol=1,zcol=3)
# estimate for large
ML$Estimates(df,taucol=1,zcol=4)
