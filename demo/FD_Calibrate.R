# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Calibrate.R",echo=TRUE)', or
# Type 'demo(FD_Calibrate)'
# R6 objects
FD <- FiniteDifference$new()
A <- Analytical$new()
# calibrate finite difference against analytical
s <- seq(from=10,to=0,by=-0.1)
x <- seq(from=-200,to=200,by=4)
FD$TerminalValue_Kinked(Vmax=10000)
FDopt <- FD$Option(s=s,x=x,mu=-15)[[1]]
Aopt <- A$Option(s=s,x=x,t=10,mu=-15)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# narrow x
x <- seq(from=-100,to=100,by=2)
FDopt <- FD$Option(x=x)[[1]]
Aopt <- A$Option(x=x)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# tweak theta (can be more or less accurate)
FDopt <- FD$Option(theta=0.878)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# change skip (smaller is faster and maybe less accurate)
FDopt <- FD$Option(theta=0.5,skip=5)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 4 times bigger sigma
FDopt <- FD$Option(sigma=60,skip=10)[[1]]
Aopt <- A$Option(sigma=60)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 2 times wider x
x <- seq(from=-200,to=200,by=4)
FDopt <- FD$Option(x=x)[[1]]
Aopt <- A$Option(x=x,t=10)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 2 times as many skips and x nodes
x <- seq(from=-200,to=200,by=2)
FDopt <- FD$Option(x=x,skip=20)[[1]]
Aopt <- A$Option(x=x,t=10)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# check option envelope
FDenv <- FD$OptionEnvelope()[[1]]
Aenv <- A$OptionEnvelope()[[1]]
err <- FDenv-Aenv
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# check decision threshold
FDdec <- FD$DecisionThreshold()
Adec <- A$DecisionThreshold()
errk <- FDdec[[1]]-Adec[[1]]
errOhat <- FDdec[[2]]-Adec[[2]]
message(paste("k difference: ",errk,"  Ohat difference: ",errOhat))
