# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_Variance.R",echo=TRUE)', or
# Type 'demo(MC_Variance)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotVariance()
# convergence
av <- A$Variance()[[1]]
mcv <- MC$Variance(paths=1)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=10)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=100)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=1000)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=10000)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=100000)[[1]]
dif <- av-mcv
max(dif)
min(dif)
mcv <- MC$Variance(paths=1000000)[[1]]
dif <- av-mcv
max(dif)
min(dif)
# plot types
MC$PlotVariance(title="type=-1",type=-1)
MC$PlotVariance(title="type=0 (default)",type=0)
# darker heat map
MC$PlotVariance(pmax=0.12,title="Darker Heat Map for Paths",type=-1)
