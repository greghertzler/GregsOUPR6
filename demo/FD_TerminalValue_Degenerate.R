# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Degenerate.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Degenerate)'
# R6 object
FD <- FiniteDifference$new()
# automatic plot with calculation
FD$set_flags(plotit=TRUE)
V <- FD$TerminalValue_Degenerate()
V <- FD$TerminalValue_Degenerate(xo=-10)
V <- FD$TerminalValue_Degenerate(xo=-10,Vmax=10,Vmin=5)
V <- FD$TerminalValue_Degenerate(xo=10,Vmax=5,Vmin=-5)
