# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Kinked.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Kinked)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Kinked()
V <- FD$TerminalValue_Kinked(xo=10,vs=1)
V <- FD$TerminalValue_Kinked(xo=10,vs=1,Vmax=15,Vmin=-5)
V <- FD$TerminalValue_Kinked(xo=-10,vs=-2,Vmax=15,Vmin=-5)
