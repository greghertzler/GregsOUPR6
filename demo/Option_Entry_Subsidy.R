# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Entry_Subsidy.R",echo=TRUE)', or
# Type 'demo(Option_Entry_Subsidy)'
# R6 object
A <- Analytical$new()
# Entry without subsidy
Entrywithout <- A$OptionEnvelope(phi=1,plotit=FALSE)$Ohat
A$PlotDecisionThreshold(title="Entry without Subsidy")
# Entry with subsidy
Entrywith <- A$OptionEnvelope(phi=1,b=10,plotit=FALSE)$Ohat
A$PlotDecisionThreshold(title="Entry with Subsidy")
# Option value of subsidy
Entrywith - Entrywithout
