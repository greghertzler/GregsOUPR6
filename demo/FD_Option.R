# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Option.R",echo=TRUE)', or
# Type 'demo(FD_Option)'
# R6 object
FD <- FiniteDifference$new()
# calculate no plot
O <- FD$Option()
# automatic plot with calculation
FD$set_flags(plotit=TRUE)
O <- FD$Option()
# new option
V <- FD$TerminalValue_Kinked(xo=-15)
O <- FD$Option(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
O <- FD$Option()
# s and x vectors manually
O <- FD$Option(s=seq(from=15,to=20,by=0.05),x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1)
O <- FD$Option()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1)
O <- FD$Option()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue()
FD$set_x_stoch_args(s=seq(from=0,to=5,by=0.05),x=seq(from=-50,to=50,by=1))
FD$PlotOption(title="My Title")
# plot types
FD$PlotOption(title="type=1",type=1)
FD$PlotOption(title="type=0 (Default)",type=0)
# exotic options
FD$undo_undo()
# no automatic plots
FD$set_flags(plotit=FALSE)
V <- FD$TerminalValue_Linear()
FD$PlotOption()
V <- FD$TerminalValue_Kinked()
FD$PlotOption()
V <- FD$TerminalValue_Stepped()
FD$PlotOption()
V <- FD$TerminalValue_Mitscherlich()
FD$PlotOption()
V <- FD$TerminalValue_Gompertz()
FD$PlotOption()
V <- FD$TerminalValue_Logistic()
FD$PlotOption()
V <- FD$TerminalValue_Transcendental()
FD$PlotOption()
V <- FD$TerminalValue_YieldIndex()
FD$PlotOption()
# custom option
mu <- FD$get_oup_params()$mu
x <- FD$get_x_stoch_args()$x
V <- (x-mu)^2
FD$set_x_stoch_args(V=V)
FD$PlotOption(title="My Custom Option")
