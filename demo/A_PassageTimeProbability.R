# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PassageTimeProbability.R",echo=TRUE)', or
# Type 'demo(A_PassageTimeProbability)'
# R6 object
A <- Analytical$new()
# type 'Pt' to print numbers
Pt <- A$PassageTimeProbability(plotit=FALSE)
# print and plot
Pt <- A$PassageTimeProbability()
# new probability
Pt <- A$PassageTimeProbability(x=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_t_stoch()
Pt <- A$PassageTimeProbability()
# t and z vectors manually
Pt <- A$PassageTimeProbability(t=seq(from=20,to=22,by=0.02),z=seq(from=-20,to=10,by=0.3))
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_t_stoch_args(t=seq(from=0,to=5,by=0.05),z=seq(from=-50,to=50,by=1),omega=0.5)
A$PlotPassageTimeProbability(title="My Title")
# plot types
A$PlotPassageTimeProbability(title="type=4 (click on the legend)",type=4)
A$PlotPassageTimeProbability(title="type=5 (click on the legend)",type=5)
A$PlotPassageTimeProbability(title="type=1 (click on the legend)",type=1)
A$PlotPassageTimeProbability(title="type=2",type=2)
A$PlotPassageTimeProbability(title="type=3 (default)",type=3)
