# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Obligation.R",echo=TRUE)', or
# Type 'demo(A_Obligation)'
# R6 object
A <- Analytical$new()
# type 'B' to print numbers
B <- A$Obligation(plotit=FALSE)
# calculate and plot
B <- A$Obligation()
# new obligation
B <- A$Obligation(y=-15,rho=0.1,mu=15)
# horizontal axes automatically
A$axes_x_stoch()
B <- A$Obligation()
# t scalar, s and x vectors manually
B <- A$Obligation(t=20,s=seq(from=15,to=20,by=0.05),x=seq(from=-50,to=20,by=0.7))
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_x_stoch_args(t=5,s=seq(from=0,to=5,by=0.05),x=seq(from=-50,to=50,by=1),y=5)
A$PlotObligation(title="My Title")
# plot types
A$PlotObligation(title="type=4",type=4)
A$PlotObligation(title="type=5",type=5)
A$PlotObligation(title="type=6",type=6)
A$PlotObligation(title="type=2",type=2)
A$PlotObligation(title="type=3 (default)",type=3)
# rotate plot
A$set_x_stoch_args(phi=1)
A$PlotObligation(title="Rotate to match Option plot")
