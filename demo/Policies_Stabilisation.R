# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_Stabilisation.R",echo=TRUE)', or
# Type 'demo(Policies_Stabilisation)'
# R6 object
A <- Analytical$new()
# no stabilisation
A$DecisionThreshold(x=seq(from=0,to=100,by=1),mu=60,sigma=50,phi=1)
# costless stabilisation
A$DecisionThreshold(sigma=25)
# costly stabilisation
A$DecisionThreshold(y=19.4)
