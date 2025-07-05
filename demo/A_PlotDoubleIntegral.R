# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotDoubleIntegral.R",echo=TRUE)', or
# Type 'demo(A_PlotDoubleIntegral)'
# R6 object
A <- Analytical$new()
# plot double integral
A$PlotDoubleIntegral()
# with custom title
A$PlotDoubleIntegral(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=15)
A$set_y_stoch_args(t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=-15,psi=1)
A$PlotDoubleIntegral()
# plot types
A$PlotDoubleIntegral(title="type=4",type=4)
A$PlotDoubleIntegral(title="type=5",type=5)
A$PlotDoubleIntegral(title="type=2",type=2)
A$PlotDoubleIntegral(title="type=3 (default)",type=3)
