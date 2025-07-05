# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_RegimeEntryToExit.R",echo=TRUE)', or
# Type 'demo(Option_RegimeEntryToExit)'
# R6 object
A <- Analytical$new()
# Exit decision threshold
k <- A$DecisionThreshold(phi=-1,plotit=FALSE)$k
A$PlotDecisionThreshold(title="Exit Decision Threshold")
# Entry decision threshold
x <- A$DecisionThreshold(phi=1,plotit=FALSE)$k
A$PlotDecisionThreshold(title="Entry Decision Threshold")
# Time from x to k
A$PassageTimePercentiles(k=k,x=x,Ppct=0.25,plotit=FALSE)
A$PlotPassageTimePercentiles(title="Percentile Times from Entry to Exit",zbeg=-20,zend=20)

