# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Variance.R",echo=TRUE)', or
# Type 'demo(A_Variance)'
# R6 object
A <- Analytical$new()
# type 'H2' to print numbers
H2 <- A$Variance(plotit=FALSE)
# print and plot
H2 <- A$Variance()
# new variance
H2 <- A$Variance(sigma=-30)
# t vector manually
H2 <- A$Variance(t=seq(from=0,to=20,by=0.2))
# horizontal axis automatically
A$axes_y_stoch()
H2 <- A$Variance()
# using set and plot
A$set_oup_params(sigma=25)
A$PlotVariance(title="My Title")
# plot types
A$PlotVariance(title="type=4 (click on p in the legend)",type=4)
A$PlotVariance(title="type=5 (click on P in the legend)",type=5)
A$PlotVariance(title="type=2",type=2)
A$PlotVariance(title="type=3 (default)",type=3)
