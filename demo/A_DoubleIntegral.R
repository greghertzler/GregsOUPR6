# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_DoubleIntegral.R",echo=TRUE)', or
# Type 'demo(A_DoubleIntegral)'
# R6 object
A <- Analytical$new()
# type 'PP' to print numbers
PP <- A$DoubleIntegral(plotit=FALSE)
# print and plot
PP <- A$DoubleIntegral()
# new double integral
PP <- A$DoubleIntegral(x=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_y_stoch()
PP <- A$DoubleIntegral()
# t and y vectors manually
PP <- A$DoubleIntegral(t=seq(from=20,to=40,by=0.2),y=seq(from=-60,to=60,by=1.2))
# integral y<=z<infinity
PP <- A$DoubleIntegral(psi=1)
# integral -infinity<z<=y
PP <- A$DoubleIntegral(psi=-1)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_y_stoch_args(t=seq(from=0,to=5,by=0.05),y=seq(from=-50,to=50,by=1),x=15)
A$PlotDoubleIntegral(title="My Title")
# plot types
A$PlotDoubleIntegral(title="type=4",type=4)
A$PlotDoubleIntegral(title="type=5",type=5)
A$PlotDoubleIntegral(title="type=2",type=2)
A$PlotDoubleIntegral(title="type=3 (default)",type=3)
