# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Exit.R",echo=TRUE)', or
# Type 'demo(Option_Exit)'
# R6 object
A <- Analytical$new()
FD <- FiniteDifference$new()
# Exit option eliminates obligation
FD$TerminalValue_Kinked(xo=0,vs=1,plotit=FALSE)
FD$PlotTerminalValue(title="Entry")
FD$TerminalValue_Linear(xo=0,vs=1,plotit=FALSE)
FD$PlotTerminalValue(title="Obligation")
FD$TerminalValue_Kinked(xo=0,vs=-1,plotit=FALSE)
FD$PlotTerminalValue(title="Exit")
# Entry, Obligation, Exit
Entry <- A$Option(phi=1,plotit=FALSE)$O
A$PlotOption(title="Entry")
Obligation <- A$Obligation(phi=-1,plotit=FALSE)$B
A$PlotObligation(title="Obligation")
Exit <- A$Option(phi=-1,plotit=FALSE)$O
A$PlotOption(title="Exit")
# Exit equals Entry minus Obligation
Difference <- Exit-(Entry-Obligation)
message(paste(sep="","Max difference: ",max(Difference),", Min difference: ",min(Difference)))
