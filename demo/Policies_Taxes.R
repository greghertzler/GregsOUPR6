# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_Taxes.R",echo=TRUE)', or
# Type 'demo(Policies_Taxes)'
# R6 object
FD <- FiniteDifference$new()
# private installation
FD$TerminalValue_Kinked(x=seq(from=0,to=100,by=1),xo=50,vs=0.05,Vmax=2,Vmin=-0.5,plotit=FALSE)
FD$DecisionThreshold(rho=0.1,mu=60,phi=1)
# government subsides installation
FD$TerminalValue_Kinked(xo=40,Vmax=2.5,Vmin=0,plotit=FALSE)
FD$DecisionThreshold()
# government pays feed-in tariff
FD$TerminalValue_Kinked(xo=50,vs=0.1,Vmax=4,Vmin=-1,plotit=FALSE)
FD$DecisionThreshold()
# no interest on power bills
FD$set_x_stoch_args(r=0)
FD$DecisionThreshold()
# more efficient solar panels
FD$set_x_stoch_args(r=0.05)
FD$TerminalValue_Kinked(xo=40,plotit=FALSE)
FD$DecisionThreshold()
