# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_ContractsToDeliver.R",echo=TRUE)', or
# Type 'demo(Finance_ContractsToDeliver)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# Read data and estimate
df<-OUPDataRead("Agric_SA_WaiteRotationTrial")
ML$Estimates(df=df,tau=1,z=10)
# expected dry matter
A$Mean(t=seq(from=0,to=1,by=0.01),s=0,x=8717.63)
# design contract with exit
A$DecisionThreshold(y=5522,phi=1,plotit=FALSE)
A$axes_x_stoch()
A$PlotDecisionThreshold()
# exit option
A$DecisionThreshold(phi=-1)
# entry without exit
A$DecisionThreshold(y=5490,phi=1,b=-165.6)
