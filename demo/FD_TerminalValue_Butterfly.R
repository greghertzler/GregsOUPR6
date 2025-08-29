# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue_Butterfly.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue_Butterfly)'
# R6 object
FD <- FiniteDifference$new()
V <- FD$TerminalValue_Butterfly()
V <- FD$TerminalValue_Butterfly(xo=-10,xm=10)
V <- FD$TerminalValue_Butterfly(Vmax=15,Vmin=-5)
V <- FD$TerminalValue_Butterfly(vs=-5)
