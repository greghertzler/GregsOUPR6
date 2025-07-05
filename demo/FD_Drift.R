# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Drift.R",echo=TRUE)', or
# Type 'demo(FD_Drift)'
# R6 object
FD <- FiniteDifference$new()
# type 'g' to print numbers
g <- FD$Drift(plotit=FALSE)
# print and plot
g <- FD$Drift()
# new drift
g <- FD$Drift(rho=1.0,mu=15)
# x vector manually
g <- FD$Drift(x=seq(from=-10,to=15,by=0.25))
# horizontal axis automatically
FD$axes_x_stoch()
g <- FD$Drift()
# using set and plot
FD$set_oup_params(rho=0.1,mu=5)
FD$PlotDrift(title="My Title")
