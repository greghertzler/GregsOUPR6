# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Resilience_Recovery.R",echo=TRUE)', or
# Type 'demo(Resilience_Recovery)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# read data
df<-myReadData("Ecosys_Kangaroos")
# Euro Kangaroos Visiting Times
oup_params <- ML$Estimates(df,taucol=1,zcol=6)
stddev <- oup_params[[3]]/(2*oup_params[[1]])^0.5
k <- oup_params[[2]]
x <- k-2*stddev
A$set_t_stoch_args(k=k,x=x,z=seq(from=0,to=2000,by=20),omega=0)
A$PlotPassageTimePercentiles(title="Euro Kangaroos")
# Red Kangaroos Visiting Times
oup_params <- ML$Estimates(df,taucol=1,zcol=2)
stddev <- oup_params[[3]]/(2*oup_params[[1]])^0.5
k <- oup_params[[2]]
x <- k-2*stddev
A$set_t_stoch_args(k=k,x=x)
A$PlotPassageTimePercentiles(title="Red Kangaroos")
# Grey Kangaroos Visiting Times
oup_params <- ML$Estimates(df,taucol=1,zcol=4)
stddev <- oup_params[[3]]/(2*oup_params[[1]])^0.5
k <- oup_params[[2]]
x <- k-2*stddev
A$set_t_stoch_args(k=k,x=x)
A$PlotPassageTimePercentiles(title="Grey Kangaroos")
# Euro Kangaroos Visiting Times equal First Passage Time if k=mu
oup_params <- ML$Estimates(df,taucol=1,zcol=6)
stddev <- oup_params[[3]]/(2*oup_params[[1]])^0.5
k <- oup_params[[2]]
x <- k-2*stddev
A$set_t_stoch_args(k=k,x=x,omega=1)
A$PlotPassageTimePercentiles(title="First Passage Times for k=mu")
