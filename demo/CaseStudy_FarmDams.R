# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_FarmDams.R",echo=TRUE)', or
# Type 'demo(CaseStudy_FarmDams)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# set up plots
A$set_plot_info(labels=FALSE)
# read data
df <- OUPDataRead("Agric_NSW_FarmDamsRiverina")
# Baseline water volume
ML$PlotTimeSeries(df=df,taucol=1,zcol=2)
ML$Estimates()
# First Passage Time from mu to 600
x <- ML$get_oup_params()$mu
k <- 600
A$set_t_stoch_args(t=seq(from=0,to=900,by=9),k=k,x=x,z=seq(from=0,to=4000,by=40))
A$PlotPassageTimePercentiles()
# 20% scenario
ML$Estimates(df=df,taucol=1,zcol=3)
x <- ML$get_oup_params()$mu
A$set_t_stoch_args(x=x)
A$PlotPassageTimePercentiles()
# 40% scenario
ML$Estimates(df=df,taucol=1,zcol=4,)
x <- ML$get_oup_params()$mu
A$set_t_stoch_args(x=x)
A$PlotPassageTimePercentiles()
