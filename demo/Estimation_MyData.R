# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_MyData.R",echo=TRUE)', or
# Type 'demo(Estimation_MyData)'
# R6 object
ML <- MaximumLikelihood$new()
# read data
datapath<-OUPDataPath()  # this is a path to your data directory
filepath<-paste0(path,"MyData.csv")  # this is your filepath
filepath
df<-read.csv(filepath)
# estimate
ML$Estimates(df,taucol=1,zcol=2)
