# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_Agriculture.R",echo=TRUE)', or
# Type 'demo(CaseStudy_Agriculture)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
ML <- OUP$get_MaximumLikelihood()
# matrices for results
coloup <- c("rho","mu","sigma","CV")
coldec <- c("k enter","Ohat enter","k exit","Ohat exit")
rowSA <- c("Clare Wheat","Orroroo Wheat","Hawker Wheat","Clare Sheep","Orroroo Sheep","Hawker Sheep")
rowNSW <- c("Cootamundra Wheat","Temora Wheat","Narrendera Wheat","Cootamundra Sheep","Temora Sheep","Narrendera Sheep")
SAoup <- matrix(0.0,6,4,dimnames = list(rowSA, coloup))
SAdec <- matrix(0.0,6,4,dimnames = list(rowSA, coldec))
NSWoup <- matrix(0.0,6,4,dimnames = list(rowNSW, coloup))
NSWdec <- matrix(0.0,6,4,dimnames = list(rowNSW, coldec))
# South Australia
# Clare
df<-OUPDataRead("Agric_SA_GMClare")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=6)
SAoup[1,1] <- oup_params[[1]]
SAoup[1,2] <- oup_params[[2]]
SAoup[1,3] <- oup_params[[3]]
SAoup[1,4] <- SAoup[1,3]/(2*SAoup[1,1])^0.5/SAoup[1,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
SAdec[1,1] <- decision[[1]]
SAdec[1,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
SAdec[1,3] <- decision[[1]]
SAdec[1,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=7)
SAoup[4,1] <- oup_params[[1]]
SAoup[4,2] <- oup_params[[2]]
SAoup[4,3] <- oup_params[[3]]
SAoup[4,4] <- SAoup[4,3]/(2*SAoup[4,1])^0.5/SAoup[4,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
SAdec[4,1] <- decision[[1]]
SAdec[4,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
SAdec[4,3] <- decision[[1]]
SAdec[4,4] <- decision[[2]]
# Orroroo
df<-OUPDataRead("Agric_SA_GMOrroroo")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=6)
SAoup[2,1] <- oup_params[[1]]
SAoup[2,2] <- oup_params[[2]]
SAoup[2,3] <- oup_params[[3]]
SAoup[2,4] <- SAoup[2,3]/(2*SAoup[2,1])^0.5/SAoup[2,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
SAdec[2,1] <- decision[[1]]
SAdec[2,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
SAdec[2,3] <- decision[[1]]
SAdec[2,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=7)
SAoup[5,1] <- oup_params[[1]]
SAoup[5,2] <- oup_params[[2]]
SAoup[5,3] <- oup_params[[3]]
SAoup[5,4] <- SAoup[5,3]/(2*SAoup[5,1])^0.5/SAoup[5,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
SAdec[5,1] <- decision[[1]]
SAdec[5,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
SAdec[5,3] <- decision[[1]]
SAdec[5,4] <- decision[[2]]
# Hawker
df<-OUPDataRead("Agric_SA_GMHawker")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=6)
SAoup[3,1] <- oup_params[[1]]
SAoup[3,2] <- oup_params[[2]]
SAoup[3,3] <- oup_params[[3]]
SAoup[3,4] <- SAoup[3,3]/(2*SAoup[3,1])^0.5/SAoup[3,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
SAdec[3,1] <- decision[[1]]
SAdec[3,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
SAdec[3,3] <- decision[[1]]
SAdec[3,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=7)
SAoup[6,1] <- oup_params[[1]]
SAoup[6,2] <- oup_params[[2]]
SAoup[6,3] <- oup_params[[3]]
SAoup[6,4] <- SAoup[6,3]/(2*SAoup[6,1])^0.5/SAoup[6,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
SAdec[6,1] <- decision[[1]]
SAdec[6,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
SAdec[6,3] <- decision[[1]]
SAdec[6,4] <- decision[[2]]
# New South Wales
# Cootamundra
df<-OUPDataRead("Agric_NSW_GMCootamundra")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=2)
NSWoup[1,1] <- oup_params[[1]]
NSWoup[1,2] <- oup_params[[2]]
NSWoup[1,3] <- oup_params[[3]]
NSWoup[1,4] <- NSWoup[1,3]/(2*NSWoup[1,1])^0.5/NSWoup[1,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
NSWdec[1,1] <- decision[[1]]
NSWdec[1,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
NSWdec[1,3] <- decision[[1]]
NSWdec[1,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=5)
NSWoup[4,1] <- oup_params[[1]]
NSWoup[4,2] <- oup_params[[2]]
NSWoup[4,3] <- oup_params[[3]]
NSWoup[4,4] <- NSWoup[4,3]/(2*NSWoup[4,1])^0.5/NSWoup[4,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
NSWdec[4,1] <- decision[[1]]
NSWdec[4,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
NSWdec[4,3] <- decision[[1]]
NSWdec[4,4] <- decision[[2]]
# Temora
df<-OUPDataRead("Agric_NSW_GMTemora")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=2)
NSWoup[2,1] <- oup_params[[1]]
NSWoup[2,2] <- oup_params[[2]]
NSWoup[2,3] <- oup_params[[3]]
NSWoup[2,4] <- NSWoup[2,3]/(2*NSWoup[2,1])^0.5/NSWoup[2,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
NSWdec[2,1] <- decision[[1]]
NSWdec[2,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
NSWdec[2,3] <- decision[[1]]
NSWdec[2,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=5)
NSWoup[5,1] <- oup_params[[1]]
NSWoup[5,2] <- oup_params[[2]]
NSWoup[5,3] <- oup_params[[3]]
NSWoup[5,4] <- NSWoup[5,3]/(2*NSWoup[5,1])^0.5/NSWoup[5,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
NSWdec[5,1] <- decision[[1]]
NSWdec[5,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
NSWdec[5,3] <- decision[[1]]
NSWdec[5,4] <- decision[[2]]
# Narrendera
df<-OUPDataRead("Agric_NSW_GMNarrendera")
# Wheat
oup_params <- ML$Estimates(df,taucol=1,zcol=2)
NSWoup[3,1] <- oup_params[[1]]
NSWoup[3,2] <- oup_params[[2]]
NSWoup[3,3] <- oup_params[[3]]
NSWoup[3,4] <- NSWoup[3,3]/(2*NSWoup[3,1])^0.5/NSWoup[3,2]
decision <- A$DecisionThreshold(x=seq(from=-850,to=1150,by=20),phi=1,b=-309)
NSWdec[3,1] <- decision[[1]]
NSWdec[3,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-278)
NSWdec[3,3] <- decision[[1]]
NSWdec[3,4] <- decision[[2]]
# Sheep
oup_params <- ML$Estimates(df,taucol=1,zcol=5)
NSWoup[6,1] <- oup_params[[1]]
NSWoup[6,2] <- oup_params[[2]]
NSWoup[6,3] <- oup_params[[3]]
NSWoup[6,4] <- NSWoup[6,3]/(2*NSWoup[6,1])^0.5/NSWoup[6,2]
decision <- A$DecisionThreshold(phi=1,b=-32)
NSWdec[6,1] <- decision[[1]]
NSWdec[6,2] <- decision[[2]]
decision <- A$DecisionThreshold(phi=-1,c=-29)
NSWdec[6,3] <- decision[[1]]
NSWdec[6,4] <- decision[[2]]
SAoup
SAdec
NSWoup
NSWdec
