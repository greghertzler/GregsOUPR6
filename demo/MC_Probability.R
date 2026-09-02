# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_Probability.R",echo=TRUE)', or
# Type 'demo(MC_Probability)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotProbability()
# convergence
ap <- A$Probability()[[1]]
mcp <- MC$Probability(paths=1)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=10)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=100)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=1000)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=10000)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=100000)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
mcp <- MC$Probability(paths=1000000)[[1]]
dif <- ap-mcp
max(dif)
min(dif)
# sum densities from right to left
MC$Probability(psi=1)
# plot types
MC$PlotProbability(title="type=1",type=1)
MC$PlotProbability(title="type=0 (default)",type=0)
