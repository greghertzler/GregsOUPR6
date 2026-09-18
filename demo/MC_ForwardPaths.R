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
# Wiener Process as plot
MC$set_oup_params(rho=0,sigma=1)
MC$PlotForwardPaths(last=100,title="Wiener Process")
# Ornstein-Uhlenbeck Process as plot
MC$set_oup_params(rho=0.5,sigma=15)
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
