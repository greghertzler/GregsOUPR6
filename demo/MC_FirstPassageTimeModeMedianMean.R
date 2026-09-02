# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_FirstPassageTimeModeMedianMean.R",echo=TRUE)', or
# Type 'demo(MC_FirstPassageTimeModeMedianMean)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotFirstPassageTimeModeMedianMean()
# convergence
am <- A$PassageTimeModeMedianMean(k=20,x=-15,omega=1)[[1]][[3]]
mcm <- MC$FirstPassageTimeModeMedianMean(paths=1)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=10)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=100)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=1000)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=10000)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=100000)[[1]][3,1]
am-mcm
mcm <- MC$FirstPassageTimeModeMedianMean(paths=1000000)[[1]][3,1]
am-mcm
# plot types
MC$axes_t_stoch()
MC$PlotFirstPassageTimeModeMedianMean(title="type=-1",type=-1)
MC$PlotFirstPassageTimeModeMedianMean(title="type=0 (default)",type=0)
