# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_Option.R",echo=TRUE)', or
# Type 'demo(MC_Option)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotOption()
# convergence
ao <- A$Option()[[1]]
mco <- MC$Option(paths=1)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=10)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=100)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=1000)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=10000)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=100000)[[1]]
dif <- ao-mco
max(dif)
min(dif)
mco <- MC$Option(paths=1000000)[[1]]
dif <- ao-mco
max(dif)
min(dif)
# double sum backward paths from left to right
MC$Option(phi=1)
# plot types
MC$PlotOption(title="type=1",type=1)
MC$PlotOption(title="type=0 (default)",type=0)
