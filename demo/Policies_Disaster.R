# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_Disaster.R",echo=TRUE)', or
# Type 'demo(Policies_Disaster)'
# R6 object
FD <- FiniteDifference$new()
# 2D plot
FD$set_plot_info(type=2)
# terminal value and option
FD$TerminalValue_Butterfly(xo=-5,xm=5,vs=-5,Vmax=100,plotit=FALSE)
FD$Option(mu=0)
