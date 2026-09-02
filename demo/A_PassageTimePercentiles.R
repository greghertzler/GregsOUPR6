# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimePercentiles.R",echo=TRUE)', or
# Type 'demo(A_PassageTimePercentiles)'
# R6 object
A <- Analytical$new()
# type 'pct' to print numbers
pct <- A$PassageTimePercentiles()
# automatic plots with calculations
A$set_flags(plotit=TRUE)
pct <- A$PassageTimePercentiles()
# equivalent of one normal standard deviation
pct <- A$PassageTimePercentiles(Ppct=0.15)
# z vector manually, quartiles
pct <- A$PassageTimePercentiles(z=seq(from=-30,to=20,by=0.5),Ppct=0.25)
# horizontal axis automatically
A$axes_t_stoch()
pct <- A$PassageTimePercentiles()
# new percentiles using set and plot
A$set_oup_params(mu=5,sigma=10)
A$set_t_stoch_args(x=-20)
A$axes_t_stoch()
A$PlotPassageTimePercentiles(title="My Title")
# plot types
A$PlotPassageTimePercentiles(title="type=1 (click on the legend)",type=1)
A$PlotPassageTimePercentiles(title="type=2 (click on the legend)",type=2)
A$PlotPassageTimePercentiles(title="type=-3",type=-3)
A$PlotPassageTimePercentiles(title="type=-2",type=-2)
A$PlotPassageTimePercentiles(title="type=-1",type=-1)
A$PlotPassageTimePercentiles(title="type=0 (default)",type=0)
