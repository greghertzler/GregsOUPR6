# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Transcendental.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Transcendental)'
# R6 object
FD <- FiniteDifference$new()
# automatic plot with calculation
FD$set_flags(plotit=TRUE)
V <- FD$TerminalValue_Transcendental()
V <- FD$TerminalValue_Transcendental(xo=-25,xi=-10,xm=5)
V <- FD$TerminalValue_Transcendental(xo=-25,xi=-25,xm=5)
V <- FD$TerminalValue_Transcendental(xo=35,xi=0,xm=-5,Vmax=5,Vmin=-5)
