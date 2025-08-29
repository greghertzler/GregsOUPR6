# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Option.R",echo=TRUE)', or
# Type 'demo(FD_Option)'
# R6 object
FD <- FiniteDifference$new()
# type 'O' to print numbers
O <- FD$Option(plotit=FALSE)
# print and plot
O <- FD$Option()
# new option
V <- FD$TerminalValue_Kinked(xo=-15,plotit=FALSE)
O <- FD$Option(rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
FD$axes_x_stoch()
O <- FD$Option()
# s and x vectors manually
O <- FD$Option(s=seq(from=15,to=20,by=0.05),x=seq(from=-100,to=100,by=2))
# entry option
V <- FD$TerminalValue_Kinked(vs=1,plotit=FALSE)
O <- FD$Option()
# exit option
V <- FD$TerminalValue_Kinked(vs=-1,plotit=FALSE)
O <- FD$Option()
# using set and plot
FD$set_oup_params(rho=0.8,mu=-5,sigma=50)
FD$set_V_kinked_args(xo=5)
V <- FD$TerminalValue(plotit=FALSE)
FD$set_x_stoch_args(s=seq(from=0,to=5,by=0.05),x=seq(from=-50,to=50,by=1))
FD$PlotOption(title="My Title")
# plot types
FD$PlotOption(title="type=4",type=4)
FD$PlotOption(title="type=5",type=5)
FD$PlotOption(title="type=2",type=2)
FD$PlotOption(title="type=3 (Default)",type=3)
# exotic options
FD$undo_undo()
V <- FD$TerminalValue_Linear(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Kinked(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Stepped(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Mitscherlich(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Gompertz(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Logistic(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_Transcendental(plotit=FALSE)
O <- FD$Option()
V <- FD$TerminalValue_YieldIndex(plotit=FALSE)
O <- FD$Option()
# custom option
mu <- FD$get_oup_params()$mu
x <- FD$get_x_stoch_args()$x
V <- (x-mu)^2
FD$set_x_stoch_args(V=V)
O <- FD$Option()
