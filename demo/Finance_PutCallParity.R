# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_PutCallParity.R",echo=TRUE)', or
# Type 'demo(Finance_PutCallParity)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# Read data and estimate
df<-read.csv("data/Finance_KansasCity_WheatFutures.csv")
ML$Estimates(df=df,tau=1,z=5,plotit=FALSE)
# 2D plots
A$set_plot_info(type=2)
# Call minus put
put <- A$Option(s=seq(from=0,to=60,by=0.6),x=seq(from=500,to=600,by=1),t=60,y=540,r=0.0002,phi=-1)[[1]]
call <- A$Option(phi=1)[[1]]
parity <- call-put
# Obligation
obligation <- A$Obligation(phi=-1)[[1]]
err <- obligation-parity
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
