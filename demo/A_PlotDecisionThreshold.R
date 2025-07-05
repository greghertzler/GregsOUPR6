# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotDecisionThreshold.R",echo=TRUE)', or
# Type 'demo(A_PlotDecisionThreshold)'
# R6 object
A <- Analytical$new()
# plot decision threshold
A$PlotDecisionThreshold()
# with custom title
A$PlotDecisionThreshold(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=15)
A$set_x_stoch_args(y=-15,phi=1)
A$PlotDecisionThreshold()
