# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_Diffusion.R",echo=TRUE)', or
# Type 'demo(A_Diffusion)'
# R6 object
A <- Analytical$new()
# type 'h2' to print numbers
h2 <- A$Diffusion(plotit=FALSE)
# print and plot
h2 <- A$Diffusion()
# new diffusion
h2 <- A$Diffusion(sigma=-30)
# z vector manually
h2 <- A$Diffusion(z=seq(from=-10,to=15,by=0.25))
# horizontal axis automatically
A$axes_z_stoch()
h2 <- A$Diffusion()
# using set and plot
A$set_oup_params(sigma=45)
A$PlotDiffusion(title="My Title")
# other plot type
A$PlotDiffusion(title="type=2",type=2)
A$PlotDiffusion(title="type=3 (default)",type=3)
