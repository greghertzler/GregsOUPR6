# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_DoubleIntegral.R",echo=TRUE)', or
# Type 'demo(MC_DoubleIntegral)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
#default
MC$PlotDoubleIntegral()
# convergence
adi <- A$DoubleIntegral()[[1]]
mcdi <- MC$DoubleIntegral(paths=1)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=10)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=100)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=1000)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=10000)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=100000)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
mcdi <- MC$DoubleIntegral(paths=1000000)[[1]]
dif <- adi-mcdi
max(dif)
min(dif)
# sum probabilities from right to left
MC$DoubleIntegral(psi=1)
# plot types
MC$PlotDoubleIntegral(title="type=1",type=1)
MC$PlotDoubleIntegral(title="type=0 (default)",type=0)
