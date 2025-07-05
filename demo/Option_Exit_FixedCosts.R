# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Exit_FixedCosts.R",echo=TRUE)', or
# Type 'demo(Option_Exit_FixedCosts)'
# R6 object
A <- Analytical$new()
# Exit without fixed costs
Exitwithout <- A$OptionEnvelope(phi=-1,plotit=FALSE)$Ohat
A$PlotDecisionThreshold(title="Exit without Fixed Costs")
# Exit with fixed costs
Exitwith <- A$OptionEnvelope(phi=-1,c=10,plotit=FALSE)$Ohat
A$PlotDecisionThreshold(title="Exit with Fixed Costs")
# Option value of bankruptcy
Exitwithout - Exitwith
