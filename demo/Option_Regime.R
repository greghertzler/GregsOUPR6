# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Regime.R",echo=TRUE)', or
# Type 'demo(Option_Regime)'
# R6 object
A <- Analytical$new()
# Exit decision threshold
A$DecisionThreshold(phi=-1,plotit=FALSE)
A$PlotDecisionThreshold(title="Exit Decision Threshold")
k <- A$DecisionThreshold(phi=-1,plotit=FALSE)
# Entry decision threshold
A$DecisionThreshold(phi=1,plotit=FALSE)
A$PlotDecisionThreshold(title="Entry Decision Threshold")
x <- A$DecisionThreshold(phi=1,plotit=FALSE)
# First Passage from Entry to Exit
A$PassageTimePercentiles(k=k,x=x,Ppct=0.25,plotit=FALSE)
A$PlotPassageTimePercentiles(title="Percentile Times from Entry to Exit",zbeg=-40,zend=20)
# First Passage from Exit to Entry
A$PassageTimePercentiles(k=x,x=k,Ppct=0.25,plotit=FALSE)
A$PlotPassageTimePercentiles(title="Percentile Times from Exit to Entry")
