# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotOption.R",echo=TRUE)', or
# Type 'demo(A_PlotOption)'
# R6 object
A <- Analytical$new()
# plot option
A$PlotOption()
# with custom title
A$PlotOption(title="My Title")
# using set and plot
A$set_oup_params(rho=0.8,mu=15,sigma=50)
A$set_x_stoch_args(t=5,s=seq(from=0,to=5,by=0.05),x=seq(from=-50,to=50,by=1),y=-5,phi=1)
A$PlotOption()
# plot types
A$PlotOption(title="type=4",type=4)
A$PlotOption(title="type=5",type=5)
A$PlotOption(title="type=2",type=2)
A$PlotOption(title="type=3 (default)",type=3)
