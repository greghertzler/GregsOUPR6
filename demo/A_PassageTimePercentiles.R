# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimePercentiles.R",echo=TRUE)', or
# Type 'demo(A_PassageTimePercentiles)'
# R6 object
A <- Analytical$new()
# type 'percentiles' to print numbers
percentiles <- A$PassageTimePercentiles(plotit=FALSE)
# print and plot
percentiles <- A$PassageTimePercentiles()
# quartile
percentiles <- A$PassageTimePercentiles(Ppct=0.25)
# z vector manually
percentiles <- A$PassageTimePercentiles(z=seq(from=-30,to=10,by=0.4))
# horizontal axis automatically
A$axes_t_stoch()
percentiles <- A$PassageTimePercentiles()
# new percentile (upper and lower) using set and plot
A$set_oup_params(rho=0,mu=5)
A$set_t_stoch_args(x=-10)
A$axes_t_stoch()
A$PlotPassageTimePercentiles(title="My Title")
# plot types
A$PlotPassageTimePercentiles(title="type=4",type=4)
A$PlotPassageTimePercentiles(title="type=5 (click on the legend)",type=5)
A$PlotPassageTimePercentiles(title="type=6 (click on the legend)",type=6)
A$PlotPassageTimePercentiles(title="type=1",type=1)
A$PlotPassageTimePercentiles(title="type=2",type=2)
A$PlotPassageTimePercentiles(title="type=3 (default)",type=3)
# finer resolution
A$set_t_stoch_args(z=seq(from=-14,to=6,by=0.2))
A$PlotPassageTimePercentiles(title="type=5 (finer resolution)",type=5)
