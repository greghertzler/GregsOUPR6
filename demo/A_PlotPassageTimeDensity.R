# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotPassageTimeDensity.R",echo=TRUE)', or
# Type 'demo(A_PlotPassageTimeDensity)'
# R6 object
A <- Analytical$new()
# plot density
A$PlotPassageTimeDensity()
# with custom title
A$PlotPassageTimeDensity(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_t_stoch_args(t=seq(from=0,to=5,by=0.05),z=seq(from=-50,to=50,by=1),omega=0.5)
A$PlotPassageTimeDensity()
# vertical axis
A$PlotPassageTimeDensity(ptmax=0.8)
# plot types
A$PlotPassageTimeDensity(title="type=4 (click on the legend)",type=4)
A$PlotPassageTimeDensity(title="type=5 (click on the legend)",type=5)
A$PlotPassageTimeDensity(title="type=1 (click on the legend)",type=1)
A$PlotPassageTimeDensity(title="type=2",type=2)
A$PlotPassageTimeDensity(title="type=3 (default)",type=3)
