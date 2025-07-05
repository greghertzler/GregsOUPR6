# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Linear.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Linear)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Linear()
V <- FD$TerminalValue_Linear(xo=-10,vs=2)
