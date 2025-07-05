# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Diffusion.R",echo=TRUE)', or
# Type 'demo(FD_Diffusion)'
# R6 object
FD <- FiniteDifference$new()
# type 'h2' to print numbers
h2 <- FD$Diffusion(plotit=FALSE)
# print and plot
h2 <- FD$Diffusion()
# new diffusion
h2 <- FD$Diffusion(sigma=-30)
# x vector manually
h2 <- FD$Diffusion(x=seq(from=-10,to=15,by=0.25))
# horizontal axis automatically
FD$axes_x_stoch()
h2 <- FD$Diffusion()
# using set and plot
FD$set_oup_params(sigma=45)
FD$PlotDiffusion(title="My Title")
# other plot type
FD$PlotDiffusion(title="type=2",type=2)
FD$PlotDiffusion(title="type=3 (default)",type=3)
