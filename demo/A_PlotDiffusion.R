# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/A_PlotDiffusion.R",echo=TRUE)', or
# Type 'demo(A_PlotDiffusion)'
# R6 object
A <- Analytical$new()
# plot diffusion
A$PlotDiffusion()
# with custom title
A$PlotDiffusion(title="My Title")
# using set and plot
A$set_oup_params(sigma=45)
A$PlotDiffusion()
# plot types
A$PlotDiffusion(title="type=2",type=2)
A$PlotDiffusion(title="type=3 (default)",type=3)
