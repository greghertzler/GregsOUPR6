# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotPassageTimeModeMedianMean.R",echo=TRUE)', or
# Type 'demo(A_PlotPassageTimeModeMedianMean)'
# R6 object
A <- Analytical$new()
# plot mode, median and mean
A$PlotPassageTimeModeMedianMean()
# with custom title
A$PlotPassageTimeModeMedianMean(title="My Title")
# new mode, median and mean using set and plot
A$set_oup_params(rho=0,mu=5)
A$set_t_stoch_args(x=-10)
A$PlotPassageTimeModeMedianMean()
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
A$set_t_stoch_args(z=seq(from=-14,to=6,by=0.2))
A$PlotPassageTimeModeMedianMean(title="type=5 (finer resolution)",type=5)
A$PlotPassageTimeModeMedianMean(title="type=6 (finer resolution)",type=6)
