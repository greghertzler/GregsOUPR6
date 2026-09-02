# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_FirstPassageTimeProbability.R",echo=TRUE)', or
# Type 'demo(MC_FirstPassageTimeProbability)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotFirstPassageTimeProbability()
# convergence
ad <- A$PassageTimeProbability(x=-30,k=0,sigma=15,omega=1)[[1]]
mcd <- MC$FirstPassageTimeProbability(paths=1)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=10)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=100)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=1000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=10000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=100000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$FirstPassageTimeProbability(paths=1000000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
# plot types
MC$set_plot_args(ptmax=0.6,zbeg=-Inf,zend=Inf)
MC$PlotFirstPassageTimeProbability(title="type=1",type=1)
MC$PlotFirstPassageTimeProbability(title="type=0 (default)",type=0)
