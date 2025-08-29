# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_OptionEnvelope.R",echo=TRUE)', or
# Type 'demo(FD_OptionEnvelope)'
# R6 object
FD <- FiniteDifference$new()
# type 'Ohat' to print numbers
Ohat <- FD$OptionEnvelope(plotit=FALSE)
# print and plot
Ohat <- FD$OptionEnvelope()
# new option
V <- FD$TerminalValue_Kinked(xo=-15,plotit=FALSE)
Ohat <- FD$OptionEnvelope(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
Ohat <- FD$OptionEnvelope()
# x vector manually
Ohat <- FD$OptionEnvelope(x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1,plotit=FALSE)
Ohat <- FD$OptionEnvelope()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1,plotit=FALSE)
Ohat <- FD$OptionEnvelope()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue(plotit=FALSE)
FD$set_x_stoch_args(x=seq(from=-50,to=50,by=1))
FD$PlotOptionEnvelope(title="My Title")
# exotic options
FD$undo_undo()
V <- FD$TerminalValue_Linear(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Kinked(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Stepped(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Mitscherlich(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Gompertz(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Logistic(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_Transcendental(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
V <- FD$TerminalValue_YieldIndex(plotit=FALSE)
Ohat <- FD$OptionEnvelope()
