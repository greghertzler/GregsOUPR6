# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotDrift.R",echo=TRUE)', or
# Type 'demo(A_PlotDrift)'
# R6 object
A <- Analytical$new()
# plot drift
A$PlotDrift()
# with custom title
A$PlotDrift(title="My Title")
# using set and plot
A$set_oup_params(rho=0.1,mu=5)
A$PlotDrift()
