# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_DecisionThreshold.R",echo=TRUE)', or
# Type 'demo(A_DecisionThreshold)'
# R6 object
A <- Analytical$new()
# no plot
A$DecisionThreshold(plotit=FALSE)
# print and plot
A$DecisionThreshold()
# new option
A$DecisionThreshold(y=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_x_stoch()
A$DecisionThreshold()
# x vector manually
A$DecisionThreshold(x=seq(from=-50,to=50,by=1))
# entry option
A$DecisionThreshold(phi=1)
# exit option
A$DecisionThreshold(phi=-1)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_x_stoch_args(y=15)
A$PlotDecisionThreshold(title="My Title")
