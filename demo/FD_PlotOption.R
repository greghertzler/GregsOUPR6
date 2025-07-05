# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_PlotOption.R",echo=TRUE)', or
# Type 'demo(FD_PlotOption)'
# R6 object
FD <- FiniteDifference$new()
# plot
FD$PlotOption()
# with custom title
FD$PlotOption(title="My Title")
# plot types
FD$PlotOption(title="type=4",type=4)
FD$PlotOption(title="type=5",type=5)
FD$PlotOption(title="type=2",type=2)
FD$PlotOption(title="type=3 (default)",type=3)
