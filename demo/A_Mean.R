# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Mean.R",echo=TRUE)', or
# Type 'demo(A_Mean)'
# R6 object
A <- Analytical$new()
# type 'G' to print numbers
G <- A$Mean(plotit=FALSE)
# print and plot
G <- A$Mean()
# new mean
G <- A$Mean(x=-5,rho=0.1,mu=15)
# t vector manually
G <- A$Mean(t=seq(from=0,to=20,by=0.2))
# horizontal axis automatically
A$axes_y_stoch()
G <- A$Mean()
# using set and plot
A$set_oup_params(rho=0.3,mu=5)
A$set_y_stoch_args(x=-10)
A$PlotMean(title="My Title")
# plot types
A$PlotMean(title="type=4 (click p in legend)",type=4)
A$PlotMean(title="type=5 (click P in legend)",type=5)
A$PlotMean(title="type=3 (default)",type=3)
