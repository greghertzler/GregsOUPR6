# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_PlotOptionEnvelope.R",echo=TRUE)', or
# Type 'demo(FD_PlotOptionEnvelope)'
# R6 object
FD <- FiniteDifference$new()
# plot
FD$PlotOptionEnvelope()
# with custom title
FD$PlotOptionEnvelope(title="My Title")
# other plot type
FD$PlotOptionEnvelope(title="type=4",type=4)
FD$PlotOptionEnvelope(title="type=3 (default)",type=3)
