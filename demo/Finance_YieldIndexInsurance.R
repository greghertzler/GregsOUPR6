# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_YieldIndexInsurance.R",echo=TRUE)', or
# Type 'demo(Finance_YieldIndexInsurance)'
# R6 objects
OUP <- OUProcess$new()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# 2D plots
FD$set_plot_info(type=2)
# Read data and estimate
filePath <- paste0(myDataPath(),"Agric_SA_GMClare.csv")
df<-read.csv(filePath)
ML$Estimates(df=df,tau=1,z=3)
# Ornstein-Uhlenbeck Process
FD$set_x_stoch_args(s=seq(from=0,to=1,by=0.01),x=seq(from=0,to=800,by=8))
FD$TerminalValue_Transcendental(xo=100,xi=250,xm=450,Vmax=5,Vmin=0)
FD$TerminalValue_YieldIndex(xo=100,xi=250,xm=450,Vmax=5,Vmin=0)
FD$TerminalValue_YieldIndex(xo=100,xi=250,xm=450,Vmax=4,Vmin=-1)
FD$Option()
