# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Resilience_Tipping.R",echo=TRUE)', or
# Type 'demo(Resilience_Tipping)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# read data
df<-OUPDataRead("Ecosys_SydneyWater")
# Warragamba Dam
ML$PlotTimeSeries(df,taucol=14,zcol=9)
ML$PlotTimeSeries(tbeg=2025)
# estimates for all observations, not affected by tbeg
ML$Estimates()
# first passage times for full dam
A$PassageTimePercentiles(k=894,x=2031,z=seq(from=0,to=2100,by=21),omega=1)
A$PlotPassageTimePercentiles(type=2)
# first passage times for dam at long-term mean
A$PassageTimePercentiles(x=1792.55,plotit=FALSE)
A$PlotPassageTimePercentiles(type=3)
A$PlotPassageTimePercentiles(type=2)
# first passage times for dam almost empty
A$PassageTimePercentiles(x=1000,plotit=FALSE)
A$PlotPassageTimePercentiles(type=3)
A$PlotPassageTimePercentiles(type=2)
