# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_OptionEnvelope.R",echo=TRUE)', or
# Type 'demo(A_OptionEnvelope)'
# R6 object
A <- Analytical$new()
# type 'Ohat' to print numbers
Ohat <- A$OptionEnvelope(plotit=FALSE)
# print and plot
Ohat <- A$OptionEnvelope()
# new option
Ohat <- A$OptionEnvelope(y=-15,rho=0.1,mu=15,sigma=25)
# horizontal axes automatically
A$axes_x_stoch()
Ohat <- A$OptionEnvelope()
# x vector manually
Ohat <- A$OptionEnvelope(x=seq(from=-50,to=50,by=1))
# entry option
Ohat <- A$OptionEnvelope(phi=1)
# exit option
Ohat <- A$OptionEnvelope(phi=-1)
# using set and plot
A$set_oup_params(rho=0.8,mu=-5)
A$set_x_stoch_args(s=seq(from=10,to=20,by=0.1),y=15)
A$PlotOptionEnvelope(title="My Title")
# other plot types
A$set_oup_params(rho=0.1)
A$PlotOptionEnvelope(title="type=4",type=4)
A$PlotOptionEnvelope(title="type=5",type=5)
A$PlotOptionEnvelope(title="type=6",type=6)
A$PlotOptionEnvelope(title="type=3 (default)",type=3)
