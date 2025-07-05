# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Stepped.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Stepped)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Stepped()
V <- FD$TerminalValue_Stepped(xo=-10,vs=1)
V <- FD$TerminalValue_Stepped(xo=-10,vs=1,Vmax=10,Vmin=5)
V <- FD$TerminalValue_Stepped(xo=-10,vs=-1,Vmax=10,Vmin=-5)
