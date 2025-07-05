# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Time_HittingFuturesStopLoss.R",echo=TRUE)', or
# Type 'demo(Time_HittingFuturesStopLoss)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# wheat futures
df <- read.csv("data/ML_WheatFutures_KansasCity.csv")
# estimate
ML$Estimates(df=df,taucol=1,zcol=5)
# goodness of fit
ML$GoodnessOfFit()
# set stop loss k and initial futures price x
A$set_t_stoch_args(k=600,x=648,t=seq(from=20,to=70,by=0.5),z=seq(from=350,to=750,by=40))
# range of initial futures prices
A$PlotPassageTimeModeMedianMean()
# plot mode, median and mean with density
A$PlotPassageTimeModeMedianMean(type=1,ptmax=0.03)
# plot mode, median and mean with probability
A$PlotPassageTimeModeMedianMean(type=2)
# calculate and print quartiles
A$PassageTimePercentiles(Ppct=0.25,plotit=FALSE)
# plot quartiles with probability
A$PlotPassageTimePercentiles(title="Passage Time Quartiles")
# plot default
A$PlotPassageTimePercentiles(title="Passage Time Quartiles",type=3)
