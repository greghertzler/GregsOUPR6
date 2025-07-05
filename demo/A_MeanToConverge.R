# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_MeanToConverge.R",echo=TRUE)', or
# Type 'demo(A_MeanToConverge)'
# R6 object
A <- Analytical$new()
# type 'G' to print numbers
G <- A$MeanToConverge(plotit=FALSE)
# print and plot
G <- A$MeanToConverge()
# initial time and slower convergence
G <- A$MeanToConverge(s=0,rho=0.1)
# horizontal axis automatically
A$axes_y_stoch()
G <- A$MeanToConverge()
# change convergence
G <- A$MeanToConverge(eps=0.1)
# using set and plot
A$set_oup_params(rho=0.3,mu=5)
A$set_y_stoch_args(x=-10,eps=0.15)
A$PlotMeanToConverge(title="My Title")
