# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_WeatherDerivative.R",echo=TRUE)', or
# Type 'demo(Finance_WeatherDerivative)'
# R6 objects
OUP <- OUProcess$new()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# read data
df<-OUPDataRead("Agric_SA_GMOrroroo")
# set terminal value
FD$set_x_stoch_args(s=seq(from=2025,to=2026,by=0.01),x=seq(from=0,to=500,by=5))
FD$TerminalValue_Kinked(xo=250,vs=-0.02,Vmax=1.5,Vmin=0)
# Ornstein-Uhlenbeck Process
ML$Estimates(df=df,taucol=1,zcol=3)
FD$PlotOption(title="Ornstein-Uhlenbeck Process")
# lower strike
FD$TerminalValue_Kinked(xo=200)
FD$PlotOption(title="Lower Strike")
# Scaled Brownian Motion
ML$Estimates(rhor=0)
FD$PlotOption(title="Scaled Brownian Motion")
