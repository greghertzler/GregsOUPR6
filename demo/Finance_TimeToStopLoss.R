# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_TimeToStopLoss.R",echo=TRUE)', or
# Type 'demo(Finance_TimeToStopLoss)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# Read data and estimate
filePath <- paste0(myDataPath(),"Finance_KansasCity_WheatFutures.csv")
df<-read.csv(filePath)
ML$Estimates(df=df,tau=1,z=5)
A$set_t_stoch_args(k=525,x=575)
A$axes_t_stoch()
# Low stop loss
A$PassageTimePercentiles()
# High stop loss
A$PassageTimePercentiles(k=625)
