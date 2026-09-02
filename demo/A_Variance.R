# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Variance.R",echo=TRUE)', or
# Type 'demo(A_Variance)'
# R6 object
A <- Analytical$new()
# type 'H2' to print numbers
H2 <- A$Variance()
# automatic plots with calculations
A$set_flags(plotit=TRUE)
H2 <- A$Variance()
# new Variance
H2 <- A$Variance(rho=1,sigma=10)
# t vector manually
H2 <- A$Variance(t=seq(from=0,to=3,by=0.03))
# horizontal axis automatically
A$axes_y_stoch()
H2 <- A$Variance()
# using set and plot
A$set_oup_params(rho=0.3)
A$set_y_stoch_args(t=seq(from=0,to=10,by=0.1))
A$PlotVariance(title="My Title")
# plot types
A$PlotVariance(title="type=-1",type=-1)
A$PlotVariance(title="type=1 (click p in legend)",type=1)
A$PlotVariance(title="type=2 (click P in legend)",type=2)
A$PlotVariance(title="type=0 (default)",type=0)
