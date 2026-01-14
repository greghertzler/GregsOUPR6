# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_StrangleAndStraddle.R",echo=TRUE)', or
# Type 'demo(Finance_StrangleAndStraddle)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# Read data and estimate
df<-OUPReadData("Finance_KansasCity_WheatFutures")
ML$Estimates(df=df,tau=1,z=5,plotit=FALSE)
# Analytical strangle
Aput <- A$Option(s=seq(from=0,to=60,by=0.6),x=seq(from=500,to=600,by=1),t=60,y=540,r=0.0002,phi=-1)[[1]]
Acall <- A$Option(phi=1)[[1]]
Astrangle <- Aput+Acall
# Finite Difference Strangle
FD$TerminalValue_Butterfly(xo=540,xm=540,Vmax=100)
FDstrangle <- FD$Option()[[1]]
err <- FDstrangle-Astrangle
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# Analytical straddle
Acall <- A$Option(y=560)[[1]]
Astrangle <- Aput+Acall
# Finite Difference Strangle
FD$TerminalValue_Butterfly(xm=560)
FDstrangle <- FD$Option()[[1]]
err <- FDstrangle-Astrangle
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
