# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Mean.R",echo=TRUE)', or
# Type 'demo(A_Mean)'
# R6 object
A <- Analytical$new()
# type 'G' to print numbers
G <- A$Mean()
# automatic plot with calculations
A$set_flags(plotit=TRUE)
G <- A$Mean()
# new mean
G <- A$Mean(x=-5,rho=0.1,mu=15)
# t vector manually
G <- A$Mean(t=seq(from=0,to=20,by=0.2))
# horizontal axis automatically
A$axes_y_stoch()
G <- A$Mean()
# using set and plot
A$set_oup_params(rho=0.3,mu=5,sigma=50)
A$set_y_stoch_args(x=-10)
A$PlotMean(title="My Title")
# plot types
A$PlotMean(title="type=-1",type=-1)
A$PlotMean(title="type=1 (click p in legend)",type=1)
A$PlotMean(title="type=2 (click P in legend)",type=2)
A$PlotMean(title="type=0 (default)",type=0)
