# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotVarianceToConverge.R",echo=TRUE)', or
# Type 'demo(A_PlotVarianceToConverge)'
# R6 object
A <- Analytical$new()
# plot variance to converge
A$PlotVarianceToConverge()
# with custom title
A$PlotVarianceToConverge(title="My Title")
# using set and plot
A$set_oup_params(sigma=50)
A$set_y_stoch_args(eps=0.15)
A$PlotVarianceToConverge()
