# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Mitscherlich.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Mitscherlich)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Mitscherlich()
V <- FD$TerminalValue_Mitscherlich(xo=-10,vr=0.1)
V <- FD$TerminalValue_Mitscherlich(xo=-10,vr=0.1,Vmax=5,Vmin=-5)
V <- FD$TerminalValue_Mitscherlich(xo=10,vr=-0.05,Vmax=5,Vmin=-5)
