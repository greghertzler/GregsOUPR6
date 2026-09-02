# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_VisitingTimeProbability.R",echo=TRUE)', or
# Type 'demo(MC_VisitingTimeProbability)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotVisitingTimeProbability()
# convergence
ad <- A$PassageTimeProbability(x=-30,k=0,sigma=15,omega=0)[[1]]
mcd <- MC$VisitingTimeProbability(paths=1)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=10)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=100)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=1000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=10000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=100000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeProbability(paths=1000000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
# plot types
MC$set_plot_args(ptmax=0.4,zbeg=-Inf,zend=Inf)
MC$PlotVisitingTimeProbability(title="type=1",type=1)
MC$PlotVisitingTimeProbability(title="type=0 (default)",type=0)
