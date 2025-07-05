# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Logistic.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Logistic)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Logistic()
V <- FD$TerminalValue_Logistic(xi=-10,vr=0.2)
V <- FD$TerminalValue_Logistic(xi=-10,vr=0.2,Vmax=6,Vmin=-6)
V <- FD$TerminalValue_Logistic(xi=-10,vr=-0.2,Vmax=3,Vmin=-2)
