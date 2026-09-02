# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_VisitingTimeModeMedianMean.R",echo=TRUE)', or
# Type 'demo(MC_VisitingTimeModeMedianMean)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotVisitingTimeModeMedianMean()
# convergence
am <- A$PassageTimeModeMedianMean(k=20,x=-15,omega=0)[[1]][[3]]
mcm <- MC$VisitingTimeModeMedianMean(paths=1)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=10)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=100)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=1000)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=10000)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=100000)[[1]][3,1]
am-mcm
mcm <- MC$VisitingTimeModeMedianMean(paths=1000000)[[1]][3,1]
am-mcm
# plot types
MC$axes_t_stoch()
MC$PlotVisitingTimeModeMedianMean(title="type=-1",type=-1)
MC$PlotVisitingTimeModeMedianMean(title="type=0 (default)",type=0)
