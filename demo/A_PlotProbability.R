# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotProbability.R",echo=TRUE)', or
# Type 'demo(A_PlotProbability)'
# R6 object
A <- Analytical$new()
# plot probability
A$PlotProbability()
# with custom title
A$PlotProbability(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=15)
A$set_y_stoch_args(t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=-15,psi=1)
A$PlotProbability()
# plot types
A$PlotProbability(title="type=4",type=4)
A$PlotProbability(title="type=5",type=5)
A$PlotProbability(title="type=2",type=2)
A$PlotProbability(title="type=3 (default)",type=3)
