# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotPassageTimeVariance.R",echo=TRUE)', or
# Type 'demo(A_PlotPassageTimeVariance)'
# R6 object
A <- Analytical$new()
# plot variance
A$PlotPassageTimeVariance()
# with custom title
A$PlotPassageTimeVariance(title="My Title")
# using set and plot
A$set_oup_params(mu=5)
A$set_t_stoch_args(x=-20)
A$PlotPassageTimeVariance()
# plot types
A$PlotPassageTimeVariance(title="type=4",type=4)
A$PlotPassageTimeVariance(title="type=3 (default)",type=3)
