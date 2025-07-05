# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Entry.R",echo=TRUE)', or
# Type 'demo(Option_Entry)'
# R6 object
A <- Analytical$new()
FD <- FiniteDifference$new()
# Entry option eliminates prohibition
FD$TerminalValue_Kinked(xo=0,vs=-1,plotit=FALSE)
FD$PlotTerminalValue(title="Exit")
FD$TerminalValue_Linear(xo=0,vs=-1,plotit=FALSE)
FD$PlotTerminalValue(title="Prohibition")
FD$TerminalValue_Kinked(xo=0,vs=1,plotit=FALSE)
FD$PlotTerminalValue(title="Entry")
# Exit, Prohibition, Entry
Exit <- A$Option(phi=-1,plotit=FALSE)$O
A$PlotOption(title="Exit")
Prohibition <- A$Obligation(phi=1,plotit=FALSE)$B
A$PlotObligation(title="Prohibition")
Entry <- A$Option(phi=1,plotit=FALSE)$O
A$PlotOption(title="Entry")
# Entry equals Exit minus Prohibition
Difference <- Entry-(Exit-Prohibition)
message(paste(sep="","Max difference: ",max(Difference),", Min difference: ",min(Difference)))
