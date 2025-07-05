# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Drift.R",echo=TRUE)', or
# Type 'demo(A_Drift)'
# R6 object
A <- Analytical$new()
# type 'g' to print numbers
g <- A$Drift(plotit=FALSE)
# print and plot
g <- A$Drift()
# new drift
g <- A$Drift(rho=1.0,mu=15)
# z vector manually
g <- A$Drift(z=seq(from=-10,to=15,by=0.25))
# horizontal axis automatically
A$axes_z_stoch()
g <- A$Drift()
# using set and plot
A$set_oup_params(rho=0.1,mu=5)
A$PlotDrift(title="My Title")
