# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_OtherOptions.R",echo=TRUE)', or
# Type 'demo(Finance_OtherOptions)'
# R6 objects
OUP <- OUProcess$new()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# automatically plot with calculations
FD$set_flags(plotit=TRUE)
# read data and estimate
df<-OUPDataRead("Agric_SA_GMClare")
ML$Estimates(df=df,taucol=1,zcol=3)
# Ornstein-Uhlenbeck Process
FD$set_x_stoch_args(s=seq(from=0,to=1,by=0.01),x=seq(from=0,to=800,by=8))
FD$TerminalValue_Transcendental(xo=100,xi=250,xm=450,Vmax=5,Vmin=0)
FD$TerminalValue_YieldIndex(xo=100,xi=250,xm=450,Vmax=5,Vmin=0)
FD$TerminalValue_YieldIndex(xo=100,xi=250,xm=450,Vmax=4,Vmin=-1)
FD$Option()
