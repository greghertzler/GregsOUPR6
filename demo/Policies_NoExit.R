# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_NoExit.R",echo=TRUE)', or
# Type 'demo(Policies_NoExit)'
# R6 object
A <- Analytical$new()
A$set_plot_info(type=2)
# default
A$DecisionThreshold(phi=-1)
# salvage value
A$DecisionThreshold(c=-10)
