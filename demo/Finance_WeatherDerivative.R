# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_WeatherDerivative.R",echo=TRUE)', or
# Type 'demo(Finance_WeatherDerivative)'
# R6 objects
OUP <- OUProcess$new()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# 2D plots
FD$set_plot_info(type=2)
# Read data and estimate
filePath <- paste0(myDataPath(),"Agric_SA_GMOrroroo.csv")
df<-read.csv(filePath)
ML$Estimates(df=df,tau=1,z=3)
# Ornstein-Uhlenbeck Process
FD$set_x_stoch_args(s=seq(from=2025,to=2026,by=0.01),x=seq(from=0,to=500,by=5))
FD$TerminalValue_Kinked(xo=250,vs=-0.02,Vmax=1.5,Vmin=0)
FD$Option()
FD$TerminalValue_Kinked(xo=200)
FD$Option()
# Scaled Brownian Motion
ML$Estimates(rhor=0)
FD$Option()
