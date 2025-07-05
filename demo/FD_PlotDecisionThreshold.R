# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_PlotDecisionThreshold.R",echo=TRUE)', or
# Type 'demo(FD_PlotDecisionThreshold)'
# R6 object
FD <- FiniteDifference$new()
# plot
FD$PlotDecisionThreshold()
# with custom title
FD$PlotDecisionThreshold(title="My Title")
