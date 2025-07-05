# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotPassageTimePercentiles.R",echo=TRUE)', or
# Type 'demo(A_PlotPassageTimePercentiles)'
# R6 object
A <- Analytical$new()
# plot percentiles
A$PlotPassageTimePercentiles()
# with custom title
A$PlotPassageTimePercentiles(title="My Title")
# quartiles
A$PassageTimePercentiles(Ppct=0.25)
# plot types
A$PlotPassageTimePercentiles(title="type=4",type=4)
A$PlotPassageTimePercentiles(title="type=5 (click on the legend)",type=5)
A$PlotPassageTimePercentiles(title="type=6 (click on the legend)",type=6)
A$PlotPassageTimePercentiles(title="type=1",type=1)
A$PlotPassageTimePercentiles(title="type=2",type=2)
A$PlotPassageTimePercentiles(title="type=3 (default)",type=3)
# finer resolution
A$set_t_stoch_args(t=seq(from=20,to=25.35,by=0.0535),z=seq(from=-30,to=30,by=0.6))
A$PlotPassageTimePercentiles(title="type=5 (finer resolution)",type=5)
