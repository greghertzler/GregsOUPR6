# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_Mean.R",echo=TRUE)', or
# Type 'demo(MC_Mean)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotMean()
# convergence
am <- A$Mean()[[1]]
mcm <- MC$Mean(paths=1)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=10)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=100)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=1000)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=10000)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=100000)[[1]]
dif <- am-mcm
max(dif)
min(dif)
mcm <- MC$Mean(paths=1000000)[[1]]
dif <- am-mcm
max(dif)
min(dif)
# plot types
MC$PlotMean(title="type=-1",type=-1)
MC$PlotMean(title="type=0 (default)",type=0)
# darker heat map
MC$PlotMean(pmax=0.12,title="Darker Heat Map for Paths",type=-1)
