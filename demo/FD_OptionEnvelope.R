# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_OptionEnvelope.R",echo=TRUE)', or
# Type 'demo(FD_OptionEnvelope)'
# R6 object
FD <- FiniteDifference$new()
# calculate no plot
Ohat <- FD$OptionEnvelope()
# automatic plot with calculation
FD$set_flags(plotit=TRUE)
Ohat <- FD$OptionEnvelope()
# new option
V <- FD$TerminalValue_Kinked(xo=-15)
Ohat <- FD$OptionEnvelope(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
Ohat <- FD$OptionEnvelope()
# x vector manually
Ohat <- FD$OptionEnvelope(x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1)
Ohat <- FD$OptionEnvelope()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1)
Ohat <- FD$OptionEnvelope()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue()
FD$set_x_stoch_args(x=seq(from=-50,to=50,by=1))
FD$PlotOptionEnvelope(title="My Title")
# exotic options
FD$undo_undo()
# no automatic plots
FD$set_flags(plotit=FALSE)
V <- FD$TerminalValue_Linear()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Kinked()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Stepped()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Mitscherlich()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Gompertz()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Logistic()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_Transcendental()
FD$PlotOptionEnvelope()
V <- FD$TerminalValue_YieldIndex()
FD$PlotOptionEnvelope()
# other plot type
FD$PlotOptionEnvelope(title="type=1",type=1)
FD$PlotOptionEnvelope(title="type=0 (default)",type=0)

