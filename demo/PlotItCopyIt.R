# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/PlotItCopyIt.R",echo=TRUE)', or
# Type 'demo(PlotItCopyIt.R)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
MC <- OUP$get_MonteCarlo()
# just numbers, no plot or copy
A$DecisionThreshold()
# automatic plots to accompany numbers
A$set_flags(plotit=TRUE)
A$DecisionThreshold()
FD$DecisionThreshold()
# automatic copy to clipboard from both calculate and plot functions
A$set_flags(plotit=FALSE,copyit=TRUE)
A$DecisionThreshold()
# overwrite clipboard with tab delimited text
ML$PlotEstimates()
