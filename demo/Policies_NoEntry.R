# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_NoEntry.R",echo=TRUE)', or
# Type 'demo(Policies_NoEntry)'
# R6 object
A <- Analytical$new()
# default
A$DecisionThreshold(phi=1)
A$PlotDecisionThreshold()
