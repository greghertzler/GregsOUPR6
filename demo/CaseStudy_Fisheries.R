# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_Fisheries.R",echo=TRUE)', or
# Type 'demo(CaseStudy_Fisheries)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
A$set_plot_info(labels=FALSE)
# read data
df<-OUPDataRead("Ecosys_SouthernBluefinTuna")
# Estimate Ecosystem Values and Gross Value Product
ML$Estimates(df,taucol=1,zcol=4)
EV <- A$get_oup_params()$mu
ML$Estimates(df,taucol=1,zcol=3)
GVP <- A$get_oup_params()$mu
EV/GVP
# Gift the ITQs
A$DecisionThreshold(phi=1,plotit=FALSE)
A$axes_x_stoch()
A$DecisionThreshold()
# Rent the ITQs
A$DecisionThreshold(b=-EV)
# Costs of fishing
A$DecisionThreshold(y=15000000)
