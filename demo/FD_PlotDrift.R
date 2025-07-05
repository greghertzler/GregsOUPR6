# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_PlotDrift.R",echo=TRUE)', or
# Type 'demo(FD_PlotDrift)'
# R6 object
FD <- FiniteDifference$new()
# plot
FD$PlotDrift()
# with custom title
FD$PlotDrift(title="My Title")
# using set and plot
FD$set_oup_params(rho=0.1,mu=5)
FD$PlotDrift()
