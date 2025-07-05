# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_VarianceToConverge.R",echo=TRUE)', or
# Type 'demo(A_VarianceToConverge)'
# R6 object
A <- Analytical$new()
# type 'H2' to print numbers
H2 <- A$VarianceToConverge(plotit=FALSE)
# print and plot
H2 <- A$VarianceToConverge()
# initial time and slower convergence
H2 <- A$VarianceToConverge(s=0,rho=0.1)
# horizontal axis automatically
A$axes_y_stoch()
H2 <- A$VarianceToConverge()
# change convergence
H2 <- A$VarianceToConverge(eps=0.1)
# using set and plot
A$set_oup_params(sigma=50)
A$set_y_stoch_args(eps=0.15)
A$PlotVarianceToConverge(title="My Title")
