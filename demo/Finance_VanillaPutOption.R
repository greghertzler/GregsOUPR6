# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_VanillaPutOption.R",echo=TRUE)', or
# Type 'demo(Finance_VanillaPutOption)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# read data and estimate
df<-OUPDataRead("Finance_KansasCity_WheatFutures")
A$set_x_stoch_args(s=seq(from=0,to=60,by=0.6),x=seq(from=500,to=600,by=1),t=60,y=540,r=0.0002,phi=-1)
# Ornstein-Uhlenbeck Process
ML$Estimates(df=df,taucol=1,zcol=5)
A$PlotOption(title="Ornstein-Uhlenbeck Process")
# Scaled Brownian Motion
ML$Estimates(rhor=0)
A$PlotOption(title="Scaled Brownian Motion")
# Stationary Process
ML$Estimates(rhor=99)
A$PlotOption(title="Stationary Process")
