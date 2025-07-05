# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimeVariance.R",echo=TRUE)', or
# Type 'demo(A_PassageTimeVariance)'
# R6 object
A <- Analytical$new()
# type 'variance' to print numbers
variance <- A$PassageTimeVariance(plotit=FALSE)
# print and plot
variance <- A$PassageTimeVariance()
# new variance
variance <- A$PassageTimeVariance(k=5,x=-20,rho=0.1,mu=15,sigma=25)
# z vector manually
variance <- A$PassageTimeVariance(z=seq(from=-30,to=10,by=0.4))
# horizontal axis automatically
A$axes_t_stoch()
variance <- A$PassageTimeVariance()
# using set and plot
A$set_oup_params(mu=5)
A$set_t_stoch_args(x=-10)
A$PlotPassageTimeVariance()
# plot types
A$PlotPassageTimeVariance(title="type=4",type=4)
A$PlotPassageTimeVariance(title="type=3 (default)",type=3)
