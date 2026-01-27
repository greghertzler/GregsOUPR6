# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_HypothesisTests.R",echo=TRUE)', or
# Type 'demo(Estimation_HypothesisTests)'
# R6 object
ML <- MaximumLikelihood$new()
# read simulated data
df<-OUPReadData("OUP_Convergence")
# estimate for rho=0.1
ML$Estimates(df=df,taucol=1,zcol=2,plotit=FALSE)
ML$GoodnessOfFit()
# estimate for rho=0.5
ML$Estimates(df=df,taucol=1,zcol=3,plotit=FALSE)
ML$GoodnessOfFit()
# estimate for rho=2.5
ML$Estimates(df,taucol=1,zcol=4,plotit=FALSE)
ML$GoodnessOfFit()
# read experimental data
df<-OUPReadData("Agric_NSW_SoilHealthHarden")
# estimate for nitrogen burn
ML$Estimates(df=df,taucol=1,zcol=2,plotit=FALSE)
ML$GoodnessOfFit()
# read commodities data
df<-OUPReadData("Finance_Commodities")
# estimate for West Texas intermediate
ML$Estimates(df=df,taucol=1,zcol=7,plotit=FALSE)
ML$GoodnessOfFit()
# scaled brownian motion
ML$Estimates(rhor=0,plotit=FALSE)
ML$LikelihoodRatioTest()
# stationary
ML$Estimates(rhor=99999,plotit=FALSE)
ML$LikelihoodRatioTest()
# 95% upper bound on rho
ML$Estimates(rhor=0.680,plotit=FALSE)
ML$LikelihoodRatioTest()
# mu irrelevant for scaled brownian motion
ML$Estimates(rhor=0,mur=99999,plotit=FALSE)
ML$LikelihoodRatioTest()
ML$Estimates(rhor=0,mur=0,plotit=FALSE)
ML$LikelihoodRatioTest()
# 95% lower bound on sigma
ML$Estimates(sigmar=0.533,plotit=FALSE)
ML$LikelihoodRatioTest()
# infinite upper bound on sigma with stationary process
df<-OUPDataRead("Agric_NSW_SoilHealthHarden")
ML$Estimates(df=df,taucol=1,zcol=2,plotit=FALSE)
ML$Estimates(sigmar=99999,plotit=FALSE)
ML$LikelihoodRatioTest()
ML$Estimates(sigmar=-99999,plotit=FALSE)
ML$LikelihoodRatioTest()
# 95% lower bound of sigma
ML$Estimates(sigmar=79.1,plotit=FALSE)
ML$LikelihoodRatioTest()


