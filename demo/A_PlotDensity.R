# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotDensity.R",echo=TRUE)', or
# Type 'demo(A_PlotDensity)'
# R6 object
A <- Analytical$new()
# plot density
A$PlotDensity()
# with custom title
A$PlotDensity(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=-5,sigma=50)
A$set_y_stoch_args(t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=15)
A$PlotDensity()
# vertical axis
A$PlotDensity(pmax=0.03)
# plot types
A$PlotDensity(title="type=4",type=4)
A$PlotDensity(title="type=5",type=5)
A$PlotDensity(title="type=2",type=2)
A$PlotDensity(title="type=3 (default)",type=3)
