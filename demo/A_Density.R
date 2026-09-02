# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Density.R",echo=TRUE)', or
# Type 'demo(A_Density)'
# R6 object
A <- Analytical$new()
# type 'p' to print numbers
p <- A$Density()
# automatic plot with calculations
A$set_flags(plotit=TRUE)
p <- A$Density()
# new density
p <- A$Density(x=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_y_stoch()
p <- A$Density()
# s, t and y manually
p <- A$Density(s=20,t=seq(from=20,to=40,by=0.2),y=seq(from=-100,to=100,by=2))
# vertical axis
A$PlotDensity(pmax=0.03)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5,sigma=50)
A$set_y_stoch_args(s=0,t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=15)
A$PlotDensity(title="My Title")
# plot types
A$PlotDensity(title="type=1",type=1)
A$PlotDensity(title="type=0 (default)",type=0)
