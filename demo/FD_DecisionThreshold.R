# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_DecisionThreshold.R",echo=TRUE)', or
# Type 'demo(FD_DecisionThreshold)'
# R6 object
FD <- FiniteDifference$new()
# no plot
FD$DecisionThreshold(plotit=FALSE)
# print and plot
FD$DecisionThreshold()
# new option
V <- FD$TerminalValue_Kinked(xo=-15,plotit=FALSE)
FD$DecisionThreshold(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
FD$DecisionThreshold()
# x vector manually
FD$DecisionThreshold(x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1,plotit=FALSE)
FD$DecisionThreshold()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1,plotit=FALSE)
FD$DecisionThreshold()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue(plotit=FALSE)
FD$set_x_stoch_args(x=seq(from=-50,to=50,by=1))
FD$PlotDecisionThreshold(title="My Title")
# exotic options
FD$default_read()
V <- FD$TerminalValue_Linear(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Kinked(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Stepped(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Mitscherlich(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Gompertz(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Logistic(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_Transcendental(plotit=FALSE)
FD$DecisionThreshold()
V <- FD$TerminalValue_YieldIndex(plotit=FALSE)
FD$DecisionThreshold()
