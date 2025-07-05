# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotVariance.R",echo=TRUE)', or
# Type 'demo(A_PlotVariance)'
# R6 object
A <- Analytical$new()
# plot variance
A$PlotVariance()
# with custom title
A$PlotVariance(title="My Title")
# using set and plot
A$set_oup_params(sigma=25)
A$PlotVariance()
# plot types
A$PlotVariance(title="type=4",type=4)
A$PlotVariance(title="type=5",type=5)
A$PlotVariance(title="type=3 (default)",type=3)
