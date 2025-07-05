# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Probability.R",echo=TRUE)', or
# Type 'demo(A_Probability)'
# R6 object
A <- Analytical$new()
# type 'P' to print numbers
P <- A$Probability(plotit=FALSE)
# print and plot
P <- A$Probability()
# new probability
P <- A$Probability(x=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_y_stoch()
P <- A$Probability()
# t and y vectors manually
P <- A$Probability(t=seq(from=20,to=40,by=0.2),y=seq(from=-60,to=60,by=1.2))
# integral y<=z<infinity
P <- A$Probability(psi=1)
# integral -infinity<z<=y
P <- A$Probability(psi=-1)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_y_stoch_args(t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=15)
A$PlotProbability(title="My Title")
# plot types
A$PlotProbability(title="type=4",type=4)
A$PlotProbability(title="type=5",type=5)
A$PlotProbability(title="type=2",type=2)
A$PlotProbability(title="type=3 (default)",type=3)
