# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_YieldIndex.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_YieldIndex)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_YieldIndex()
V <- FD$TerminalValue_YieldIndex(xo=-25,xi=-10,xm=5)
V <- FD$TerminalValue_YieldIndex(xo=-25,xi=-25,xm=5,Vmax=9,Vmin=-1)
V <- FD$TerminalValue_YieldIndex(xo=35,xi=0,xm=-5,Vmax=6,Vmin=-1)
