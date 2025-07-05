# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Option.R",echo=TRUE)', or
# Type 'demo(A_Option)'
# R6 object
A <- Analytical$new()
# type 'O' to print numbers
O <- A$Option(plotit=FALSE)
# calculate and plot
O <- A$Option()
# new option
O <- A$Option(y=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_x_stoch()
O <- A$Option()
# t scalar, s and x vectors manually
O <- A$Option(t=20,s=seq(from=15,to=20,by=0.05),x=seq(from=-50,to=20,by=0.7))
# integral -infinity<z<=y (entry option)
O <- A$Option(phi=1)
# integral y<=z<infinity (exit option)
O <- A$Option(phi=-1)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5,sigma=50)
A$set_x_stoch_args(t=5,s=seq(from=0,to=5,by=0.05),x=seq(from=-50,to=50,by=1),y=5)
A$PlotOption(title="My Title")
# plot types
A$PlotOption(title="type=4",type=4)
A$PlotOption(title="type=5",type=5)
A$PlotOption(title="type=2",type=2)
A$PlotOption(title="type=3 (default)",type=3)
