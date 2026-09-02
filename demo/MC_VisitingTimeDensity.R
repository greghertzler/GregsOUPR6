# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/MC_VisitingTimeDensity.R",echo=TRUE)', or
# Type 'demo(MC_VisitingTimeDensity)'
# R6 objects
OUP <- OUProcess$new()
MC <- OUP$get_MonteCarlo()
A <- OUP$get_Analytical()
# default
MC$PlotVisitingTimeDensity()
# convergence
ad <- A$PassageTimeDensity(x=-30,k=0,sigma=15,omega=0)[[1]]
mcd <- MC$VisitingTimeDensity(paths=1)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=10)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=100)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=1000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=10000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=100000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
mcd <- MC$VisitingTimeDensity(paths=1000000)[[1]]
dif <- ad-mcd
max(dif)
min(dif)
sum(dif)
# plot types
MC$set_plot_args(ptmax=0.4,zbeg=-Inf,zend=Inf)
MC$PlotVisitingTimeDensity(title="type=1",type=1)
MC$PlotVisitingTimeDensity(title="type=0 (default)",type=0)
