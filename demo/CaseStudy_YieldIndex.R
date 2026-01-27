# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_YieldIndex.R",echo=TRUE)', or
# Type 'demo(CaseStudy_YieldIndex)'
# R6 object
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
FD$set_plot_info(labels=FALSE)
# read data
df<-OUPDataRead("Agric_SA_WaiteRotationTrial")
# estimate Apr-Oct rainfall
ML$Estimates(df,taucol=1,zcol=25)
# transcendental yield
FD$TerminalValue_Transcendental(x=seq(from=0,to=900,by=9),xo=0,xi=316,xm=485,Vmax=854)
# yield index
FD$TerminalValue_YieldIndex(x=seq(from=0,to=900,by=9),xo=0,xi=316,xm=485,Vmax=654,Vmin=-200)
#option price
FD$set_plot_info(type=2)
FD$Option(s=seq(from=0,to=1,by=0.01))
# estimate continuous wheat
ML$Estimates(df,taucol=1,zcol=11)
# put option
A$Option(s=seq(from=0,to=1,by=0.01),x=seq(from=230,to=1030,by=8),t=1,654)
