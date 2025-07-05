# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotMean.R",echo=TRUE)', or
# Type 'demo(A_PlotMean)'
# R6 object
A <- Analytical$new()
# plot mean
A$PlotMean()
# with custom title
A$PlotMean(title="My Title")
# using set and plot
A$set_oup_params(rho=0.3,mu=5)
A$set_y_stoch_args(x=-10)
A$PlotMean()
# plot types
A$PlotMean(title="type=4",type=4)
A$PlotMean(title="type=5",type=5)
A$PlotMean(title="type=3 (default)",type=3)
