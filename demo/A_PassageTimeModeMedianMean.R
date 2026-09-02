# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimeModeMedianMean.R",echo=TRUE)', or
# Type 'demo(A_PassageTimeModeMedianMean)'
# R6 object
A <- Analytical$new()
# type 'mmm' to print numbers
mmm <- A$PassageTimeModeMedianMean()
# automatic plots with calculations
A$set_flags(plotit=TRUE)
mmm <- A$PassageTimeModeMedianMean()
# z vector manually
mmm <- A$PassageTimeModeMedianMean(z=seq(from=-30,to=20,by=0.5))
# horizontal axis automatically
A$axes_t_stoch()
mmm <- A$PassageTimeModeMedianMean()
# new percentiles using set and plot
A$set_oup_params(mu=5,sigma=10)
A$set_t_stoch_args(x=-20)
A$axes_t_stoch()
A$PlotPassageTimeModeMedianMean(title="My Title")
# plot types
A$PlotPassageTimeModeMedianMean(title="type=1 (click on the legend)",type=1)
A$PlotPassageTimeModeMedianMean(title="type=2 (click on the legend)",type=2)
A$PlotPassageTimeModeMedianMean(title="type=-3",type=-3)
A$PlotPassageTimeModeMedianMean(title="type=-2",type=-2)
A$PlotPassageTimeModeMedianMean(title="type=-1",type=-1)
A$PlotPassageTimeModeMedianMean(title="type=0 (default)",type=0)
