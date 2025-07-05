# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Option_OnFutures.R",echo=TRUE)', or
# Type 'demo(Option_OnFutures)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# Wheat futures
df <- read.csv("data/ML_WheatFutures_KansasCity.csv")
# Estimate
ML$Estimates(df=df,taucol=1,zcol=5)
# 2D
A$set_plot_info(type=2)
# European put
A$Option(s=seq(from=0,to=100,by=1),x=seq(from=500,to=600,by=1),t=100,y=520,r=0.000137,phi=-1)
# European call
A$Option(phi=1)
