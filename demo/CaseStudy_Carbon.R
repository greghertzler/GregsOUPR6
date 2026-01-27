# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/CaseStudy_Carbon.R",echo=TRUE)', or
# Type 'demo(CaseStudy_Carbon)'
# R6 object
OUP <- OUProcess$new()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
FD$set_plot_info(labels=FALSE)
# read data
df<-OUPDataRead("Climate_CarbonCredits_EECXM")
# Estimate Closing prices of EU carbon credits
ML$Estimates(df,taucol=1,zcol=5)
# keep stochastics, but substitute ACCU spot price for mu
FD$set_oup_params(mu=36.25)
# set annuities
ann100 <- 3.90497
est100 <- 50.3831
ann25 <- 3.73978
est25 <- 70.9525
xo100 <- est100/ann100
xo25 <- est25/ann25
# entry 100 year
FD$TerminalValue_Kinked(x=seq(from=-50,to=100,by=1.5),xo=xo100,vs=ann100,Vmax=10000,Vmin=-est100,plotit=FALSE)
FD$DecisionThreshold()
# exit 100 year
FD$TerminalValue_Kinked(xo=0,vs=-ann100,Vmin=0,plotit=FALSE)
FD$DecisionThreshold()
# entry 25 year
FD$TerminalValue_Kinked(xo=xo25,vs=ann25,Vmin=-est25,plotit=FALSE)
FD$DecisionThreshold()
