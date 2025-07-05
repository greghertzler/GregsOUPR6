# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_Butterfly.R",echo=TRUE)', or
# Type 'demo(Option_Butterfly)'
# R6 object
FD <- FiniteDifference$new()
# construct terminal value
Vleft <- FD$TerminalValue_Kinked(x=seq(from=-30,to=50,by=0.8),xo=0,vs=-1,plotit=FALSE)$V
Vright <- FD$TerminalValue_Kinked(xo=20,vs=1,plotit=FALSE)$V
V <- Vleft+Vright
FD$set_x_stoch_args(V=V)
FD$TerminalValue(name="Butterfly")
# Butterfly option
FD$Option(plotit=FALSE)
FD$PlotOption(title="Butterfly Option")
FD$PlotOptionEnvelope(title="Butterfly Option")
FD$PlotDecisionThreshold(title="Butterfly Exit")
FD$DecisionThreshold(plotit=FALSE,phi=1)
FD$PlotDecisionThreshold(title="Butterfly Enter")
