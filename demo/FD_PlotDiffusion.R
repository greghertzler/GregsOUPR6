# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_PlotDiffusion.R",echo=TRUE)', or
# Type 'demo(FD_PlotDiffusion)'
# R6 object
FD <- FiniteDifference$new()
# plot
FD$PlotDiffusion()
# with custom title
FD$PlotDiffusion(title="My Title")
# using set and plot
FD$set_oup_params(sigma=45)
FD$PlotDiffusion()
