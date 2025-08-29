# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_VanillaPutOption.R",echo=TRUE)', or
# Type 'demo(Finance_VanillaPutOption)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# Read data and estimate
df<-read.csv("data/Finance_KansasCity_WheatFutures.csv")
# 2D plots
A$set_plot_info(type=2)
# Ornstein-Uhlenbeck Process
ML$Estimates(df=df,tau=1,z=5,plotit=FALSE)
A$Option(s=seq(from=0,to=60,by=0.6),x=seq(from=500,to=600,by=1),t=60,y=540,r=0.0002,phi=-1,plotit=FALSE)
A$PlotOption(title="Ornstein-Uhlenbeck Process")
# Scaled Brownian Motion
ML$Estimates(rhor=0,plotit=FALSE)
A$PlotOption(title="Scaled Brownian Motion")
# Stationary Process
ML$Estimates(rhor=99,plotit=FALSE)
A$PlotOption(title="Stationary Process")
