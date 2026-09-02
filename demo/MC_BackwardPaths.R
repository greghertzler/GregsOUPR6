# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_BackwardPaths.R",echo=TRUE)', or
# Type 'demo(A_BackwardPaths)'
# R6 object
MC <- MonteCarlo$new()
# automatic plots with calculations
MC$set_flags(plotit=TRUE)
# default
MC$BackwardPaths()
# more paths
MC$BackwardPaths(paths=500)
# different seed for random numbers
MC$BackwardPaths(seed=123)
# no automatic plots
MC$set_flags(plotit=FALSE)
# Integral equation and 4th order Runge-Kutta convergence
# skip is same in both to reuse pseudo-random numbers
options(digits=10)
bie <- MC$BackwardPaths(sigma=50)[[1]]
brk <- MC$BackwardPaths(method=4)[[1]]
dif <- bie-brk
max(dif)
min(dif)
bie <- MC$BackwardPaths(skip=2)[[1]]
brk <- MC$BackwardPaths()[[1]]
dif <- bie-brk
max(dif)
min(dif)
# Wiener Process as plot
MC$set_oup_params(rho=0,sigma=1)
MC$PlotBackwardPaths(last=100,title="Wiener Process")
# Ornstein-Uhlenbeck Process as plot
MC$set_oup_params(rho=0.5,sigma=15)
MC$set_path_args(method=1,skip=1)
MC$PlotBackwardPaths(title="Ornstein-Uhlenbeck Process")
# custom labels
MC$PlotBackwardPaths(title="MyTitle",xaxis="MyxAxis",yaxis="MyyAxis")
# one path
MC$PlotBackwardPaths(title="One Path",first=5,last=5)
# one point
MC$PlotBackwardPaths(title="One Path One Time",sbeg=5,send=5)
# first 10 out of 500 paths
MC$PlotBackwardPaths(first=1,last=10)
# last 10 out of 500 paths
MC$PlotBackwardPaths(first=491,last=500)
# plot types
MC$PlotBackwardPaths(title="OUP type=-3",type=-3)
MC$PlotBackwardPaths(title="OUP type=-2",type=-2)
MC$PlotBackwardPaths(title="OUP type=-1",type=-1)
MC$PlotBackwardPaths(title="OUP type=0 (default)",type=0)
