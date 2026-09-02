# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_Density.R",echo=TRUE)', or
# Type 'demo(MC_Density)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotDensity()
# convergence
ad <- A$Density()[[1]]
mcd <- MC$Density(paths=1)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=10)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=100)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=1000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=10000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=100000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
mcd <- MC$Density(paths=1000000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
# plot types
MC$PlotDensity(title="type=1",type=1)
MC$PlotDensity(title="type=0 (default)",type=0)
# heat map is the density using colors
MC$PlotDensity(pmax=0.12,title="Heat Map is the Density",type=1)
