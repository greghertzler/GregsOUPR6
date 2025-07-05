# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimeMean.R",echo=TRUE)', or
# Type 'demo(A_PassageTimeMean)'
# R6 object
A <- Analytical$new()
# type 'mean' to print numbers
mean <- A$PassageTimeMean(plotit=FALSE)
# print and plot
mean <- A$PassageTimeMean()
# new mean
mean <- A$PassageTimeMean(k=5,x=-20,rho=0.1,mu=15)
# z vector manually
mean <- A$PassageTimeMean(z=seq(from=-30,to=10,by=0.4))
# horizontal axis automatically
A$axes_t_stoch()
mean <- A$PassageTimeMean()
# mean does not exist using set and plot
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
# coarser and finer resolution
A$set_oup_params(rho=0.5)
A$PlotPassageTimeModeMedianMean(title="type=5 (coarser resolution)",type=5)
A$PlotPassageTimeModeMedianMean(title="type=6 (coarser resolution)",type=6)
A$set_t_stoch_args(z=seq(from=-51,to=29,by=0.8))
A$PlotPassageTimeModeMedianMean(title="type=5 (finer resolution)",type=5)
A$PlotPassageTimeModeMedianMean(title="type=6 (finer resolution)",type=6)
