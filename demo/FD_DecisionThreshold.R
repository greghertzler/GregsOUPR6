# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_DecisionThreshold.R",echo=TRUE)', or
# Type 'demo(FD_DecisionThreshold)'
# R6 object
FD <- FiniteDifference$new()
# with automatic plots
FD$set_flags(plotit=TRUE)
FD$DecisionThreshold()
# new option
V <- FD$TerminalValue_Kinked(xo=-15)
FD$DecisionThreshold(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
FD$DecisionThreshold()
# x vector manually
FD$DecisionThreshold(x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1)
FD$DecisionThreshold()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1)
FD$DecisionThreshold()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue()
FD$set_x_stoch_args(x=seq(from=-50,to=50,by=1))
FD$PlotDecisionThreshold(title="My Title")
# exotic options with no automatic plots
FD$set_flags(plotit=FALSE)
FD$undo_undo()
V <- FD$TerminalValue_Linear()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Kinked()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Stepped()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Mitscherlich()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Gompertz()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Logistic()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_Transcendental()
FD$PlotDecisionThreshold()
V <- FD$TerminalValue_YieldIndex()
FD$PlotDecisionThreshold()
