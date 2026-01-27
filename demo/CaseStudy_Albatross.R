# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_Albatross.R",echo=TRUE)', or
# Type 'demo(CaseStudy_Albatross)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
A$set_plot_info(labels=FALSE)
# read data
df <- OUPDataRead("Ecosys_Albatross")
# Wandering Albatross egg counts
ML$Estimates(df,taucol=1,zcol=2)
# Visiting Time to recover to mu
k <- ML$get_oup_params()$mu
A$set_t_stoch_args(t=seq(from=0,to=7,by=0.07),k=k,x=2,z=seq(from=0,to=15,by=0.15))
A$PlotPassageTimePercentiles()
A$PlotPassageTimePercentiles(type=2)
