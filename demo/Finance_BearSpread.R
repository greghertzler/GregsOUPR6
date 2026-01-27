# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_BearSpread.R",echo=TRUE)', or
# Type 'demo(Finance_BearSpread)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
# 2D plots
A$set_plot_info(type=2)
# Read data and estimate
df<-OUPDataRead("Finance_Commodities")
ML$Estimates(df=df,tau=8,z=6)
# Analytical bear spread
puthi <- A$Option(s=seq(from=2025,to=2026,by=0.01),x=seq(from=0,to=150,by=1.5),t=2025,y=70,phi=-1)[[1]]
putlo <- A$Option(y=50)[[1]]
Aspread <- puthi-putlo
# Finite Difference bear spread
FD$TerminalValue_Kinked(xo=70,vs=-1,Vmax=20,Vmin=0)
FDspread <- FD$Option()[[1]]
err <- FDspread-Aspread
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
# Scaled Brownian Motion
ML$Estimates(rhor=0)
FD$Option()
