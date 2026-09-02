# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_Bankruptcy.R",echo=TRUE)', or
# Type 'demo(Policies_Bankruptcy)'
# R6 object
A <- Analytical$new()
# automatic plots with calculations
A$set_flags(plotit=TRUE)
# entry with fixed costs
A$DecisionThreshold(y=5,sigma=30,phi=1)
A$undo_clear()
A$DecisionThreshold(mu=10,y=0)
A$undo_undo()
# exit with fixed costs and bankruptcy
A$DecisionThreshold(phi=-1)
# no bankruptcy
A$DecisionThreshold(y=0,c=5)
