# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_BoundedPaths.R",echo=TRUE)', or
# Type 'demo(A_BoundedPaths)'
# R6 object
MC <- MonteCarlo$new()
# automatic plots with calculations
MC$set_flags(plotit=TRUE)
# default
MC$BoundedPaths()
# more paths
MC$BoundedPaths(paths=500)
# different seed for random numbers
MC$BoundedPaths(seed=123)
# no automatic plots
MC$set_flags(plotit=FALSE)
# Integral equation and 4th order Runge-Kutta convergence
# skip is same in both to reuse pseudo-random numbers
options(digits=10)
fie <- MC$BoundedPaths(sigma=50)[[1]]
frk <- MC$BoundedPaths(method=4)[[1]]
difidx <- !is.na(fie) & !is.na(frk)
dif <- fie[difidx]-frk[difidx]
max(dif)
min(dif)
fie <- MC$BoundedPaths(skip=2)[[1]]
frk <- MC$BoundedPaths()[[1]]
difidx <- !is.na(fie) & !is.na(frk)
dif <- fie[difidx]-frk[difidx]
max(dif)
min(dif)
# Wiener Process as plot
MC$set_oup_params(rho=0,sigma=1)
MC$PlotBoundedPaths(last=100,title="Wiener Process")
# Ornstein-Uhlenbeck Process as plot
MC$set_oup_params(rho=0.5,sigma=15)
MC$set_path_args(method=1,skip=1)
MC$PlotBoundedPaths(title="Ornstein-Uhlenbeck Process")
# custom labels
MC$PlotBoundedPaths(title="MyTitle",xaxis="MyxAxis",yaxis="MyyAxis")
# one path
MC$PlotBoundedPaths(title="One Path",first=2,last=2)
# one point
MC$PlotBoundedPaths(title="One Path One Time",tbeg=5,tend=5)
# first 10 out of 500 paths
MC$PlotBoundedPaths(first=1,last=10)
# last 10 out of 500 paths
MC$PlotBoundedPaths(first=491,last=500)
# plot types
MC$PlotBoundedPaths(title="type=-3",type=-3)
MC$PlotBoundedPaths(title="type=-2",type=-2)
MC$PlotBoundedPaths(title="type=-1",type=-1)
MC$PlotBoundedPaths(title="type=0 (default)",type=0)
