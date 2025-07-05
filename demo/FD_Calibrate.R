# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Calibrate.R",echo=TRUE)', or
# Type 'demo(FD_Calibrate)'
# R6 objects
FD <- FiniteDifference$new()
A <- Analytical$new()
# calibrate finite difference against analytical
s <- seq(from=10,to=0,by=-0.1)
x <- seq(from=-200,to=200,by=4)
FDopt <- FD$Option(s=s,x=x,plotit=FALSE)[[1]]
Aopt <- A$Option(s=s,x=x,t=10,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# narrow x
x <- seq(from=-100,to=100,by=2)
FDopt <- FD$Option(x=x,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# tweak theta (can be more or less accurate)
FDopt <- FD$Option(theta=0.878,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# change skip (smaller is faster and maybe less accurate)
FDopt <- FD$Option(theta=0.5,skip=5,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 4 times bigger sigma
FDopt <- FD$Option(sigma=60,skip=10,plotit=FALSE)[[1]]
Aopt <- A$Option(sigma=60,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 2 times wider x
x <- seq(from=-200,to=200,by=4)
FDopt <- FD$Option(x=x,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,t=10,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# 2 times as many skips and x nodes
x <- seq(from=-200,to=200,by=2)
n <- length(x)
FDopt <- FD$Option(x=x,skip=20,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,t=10,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# check option envelope
FDenv <- FD$OptionEnvelope(plotit=FALSE)[[1]]
Aenv <- A$OptionEnvelope(plotit=FALSE)[[1]]
err <- FDenv-Aenv
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# check decision threshold
FDdec <- FD$DecisionThreshold(plotit=FALSE)
Adec <- A$DecisionThreshold(plotit=FALSE)
errk <- FDdec[[1]]-Adec[[1]]
errOhat <- FDdec[[2]]-Adec[[2]]
message(paste("k difference: ",errk,"  Ohat difference: ",errOhat))
