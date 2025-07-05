# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimeMode.R",echo=TRUE)', or
# Type 'demo(A_PassageTimeMode)'
# R6 object
A <- Analytical$new()
# type 'mode' to print numbers
mode <- A$PassageTimeMode(plotit=FALSE)
# print and plot
mode <- A$PassageTimeMode()
# new mode
mode <- A$PassageTimeMode(k=5,x=-20,rho=0.1,mu=15)
# z vector manually
mode <- A$PassageTimeMode(z=seq(from=-30,to=10,by=0.4))
# horizontal axis automatically
A$axes_t_stoch()
mode <- A$PassageTimeMode()
# new mode (median but no mean) using set and plot
A$set_oup_params(rho=0,mu=5)
A$set_t_stoch_args(x=-10)
A$PlotPassageTimeModeMedianMean(title="My Title")
# plot types
A$PlotPassageTimeModeMedianMean(title="type=4",type=4)
A$PlotPassageTimeModeMedianMean(title="type=5 (click on the legend)",type=5)
A$PlotPassageTimeModeMedianMean(title="type=6 (click on the legend)",type=6)
A$PlotPassageTimeModeMedianMean(title="type=1",type=1)
A$PlotPassageTimeModeMedianMean(title="type=2",type=2)
A$PlotPassageTimeModeMedianMean(title="type=3 (default)",type=3)
# finer resolution
A$set_t_stoch_args(z=seq(from=-51,to=29,by=0.8))
A$PlotPassageTimeModeMedianMean(title="type=5 (finer resolution)",type=5)
