# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_ForwardPaths.R",echo=TRUE)', or
# Type 'demo(A_ForwardPaths)'
# R6 object
MC <- MonteCarlo$new()
# automatic plots with calculations
MC$set_flags(plotit=TRUE)
# default
MC$ForwardPaths()
# more paths
MC$ForwardPaths(paths=500)
# different seed for random numbers
MC$ForwardPaths(seed=123)
# not automatic plots
MC$set_flags(plotit=FALSE)
# Integral equation and 4th order Runge-Kutta convergence
# skip is same in both to reuse pseudo-random numbers
options(digits=10)
fie <- MC$ForwardPaths(sigma=50)[[1]]
frk <- MC$ForwardPaths(method=4)[[1]]
dif <- fie-frk
max(dif)
min(dif)
fie <- MC$ForwardPaths(skip=2)[[1]]
frk <- MC$ForwardPaths()[[1]]
dif <- fie-frk
max(dif)
min(dif)
# Wiener Process as plot
MC$set_oup_params(rho=0,sigma=1)
MC$PlotForwardPaths(last=100,title="Wiener Process")
# Ornstein-Uhlenbeck Process as plot
MC$set_oup_params(rho=0.5,sigma=15)
MC$set_path_args(method=1,skip=1)
MC$PlotForwardPaths(title="Ornstein-Uhlenbeck Process")
# custom labels
MC$PlotForwardPaths(title="MyTitle",xaxis="MyxAxis",yaxis="MyyAxis")
# one path
MC$PlotForwardPaths(title="One Path",first=5,last=5)
# one point
MC$PlotForwardPaths(title="One Path One Time",tbeg=5,tend=5)
# first 10 out of 500 paths
MC$PlotForwardPaths(first=1,last=10)
# last 10 out of 500 paths
MC$PlotForwardPaths(first=491,last=500)
# plot types
MC$PlotForwardPaths(title="type=-3",type=-3)
MC$PlotForwardPaths(title="type=-2",type=-2)
MC$PlotForwardPaths(title="type=-1",type=-1)
MC$PlotForwardPaths(title="type=0 (default)",type=0)
