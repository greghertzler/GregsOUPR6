# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Gompertz.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Gompertz)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Gompertz()
V <- FD$TerminalValue_Gompertz(xi=0,vr=0.1)
V <- FD$TerminalValue_Gompertz(xi=0,vr=0.1,Vmax=2,Vmin=1)
V <- FD$TerminalValue_Gompertz(xi=-5,vr=-0.1,Vmax=2,Vmin=-1)
