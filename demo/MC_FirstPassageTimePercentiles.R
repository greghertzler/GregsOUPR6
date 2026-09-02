# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_FirstPassageTimePercentiles.R",echo=TRUE)', or
# Type 'demo(MC_FirstPassageTimePercentiles)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotFirstPassageTimePercentiles()
# convergence
am <- A$PassageTimeModeMedianMean(k=20,x=-15,omega=1)[[1]][[2]]
mcm <- MC$FirstPassageTimePercentiles(paths=1)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=10)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=100)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=1000)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=10000)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=100000)[[1]][2,1]
am-mcm
mcm <- MC$FirstPassageTimePercentiles(paths=1000000)[[1]][2,1]
am-mcm
# plot types
MC$axes_t_stoch()
MC$PlotFirstPassageTimePercentiles(title="type=-1",type=-1)
MC$PlotFirstPassageTimePercentiles(title="type=0 (default)",type=0)
