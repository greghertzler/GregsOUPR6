# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_VisitingTimePercentiles.R",echo=TRUE)', or
# Type 'demo(MC_VisitingTimePercentiles)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotVisitingTimePercentiles()
# convergence
am <- A$PassageTimeModeMedianMean(k=20,x=-15,omega=0)[[1]][[2]]
mcm <- MC$VisitingTimePercentiles(paths=1)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=10)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=100)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=1000)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=10000)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=100000)[[1]][2,1]
am-mcm
mcm <- MC$VisitingTimePercentiles(paths=1000000)[[1]][2,1]
am-mcm
# plot types
MC$axes_t_stoch()
MC$PlotVisitingTimePercentiles(title="type=-1",type=-1)
MC$PlotVisitingTimePercentiles(title="type=0 (default)",type=0)
