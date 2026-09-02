# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Adoption_Obligation.R",echo=TRUE)', or
# Type 'demo(Adoption_Obligation)'
# R6 object
A <- Analytical$new()
# automatic plots with calculations
A$set_flags(plotit=TRUE)
# entry option
A$Option(phi=1)
# exit option
A$Option(phi=-1)
# obligation
A$Obligation()
# prohibition
A$Obligation(phi=1)
