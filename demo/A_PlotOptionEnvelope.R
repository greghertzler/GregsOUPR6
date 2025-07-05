# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotOptionEnvelope.R",echo=TRUE)', or
# Type 'demo(A_PlotOptionEnvelope)'
# R6 object
A <- Analytical$new()
# plot option envelope
A$PlotOptionEnvelope()
# with custom title
A$PlotOptionEnvelope(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=15)
A$set_x_stoch_args(y=-15,phi=1)
A$PlotOptionEnvelope()
# other plot type
A$PlotOptionEnvelope(title="type=4",type=4)
A$PlotOptionEnvelope(title="type=3 (default)",type=3)
