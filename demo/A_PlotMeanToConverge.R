# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotMeanToConverge.R",echo=TRUE)', or
# Type 'demo(A_PlotMeanToConverge)'
# R6 object
A <- Analytical$new()
# plot mean to converge
A$PlotMeanToConverge()
# with custom title
A$PlotMeanToConverge(title="My Title")
# using set and plot
A$set_oup_params(rho=0.3,mu=5)
A$set_y_stoch_args(x=-10,eps=0.15)
A$PlotMeanToConverge()
