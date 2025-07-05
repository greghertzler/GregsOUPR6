# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_IndexInsurance.R",echo=TRUE)', or
# Type 'demo(Option_IndexInsurance)'
# R6 object
FD <- FiniteDifference$new()
# Terminal value
FD$set_oup_params(rho=0.1,mu=15,sigma=10)
FD$set_x_stoch_args(x=seq(from=-50,to=150,by=1))
V <- FD$TerminalValue_YieldIndex(xo=0,xi=10,xm=20,Vmax=6,Vmin=-4)$V
FD$set_x_stoch_args(V=V)
# Option
FD$PlotOption(title="Index Insurance",xbeg=0,xend=50)
