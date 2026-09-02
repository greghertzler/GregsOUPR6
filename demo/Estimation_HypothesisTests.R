# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_HypothesisTests.R",echo=TRUE)', or
# Type 'demo(Estimation_HypothesisTests)'
# R6 object
ML <- MaximumLikelihood$new()
# read simulated data
df<-OUPDataRead("OUP_Convergence")
# estimate for rho=0.1
ML$Estimates(df=df,taucol=1,zcol=2)
ML$GoodnessOfFit()
# estimate for rho=0.5
ML$Estimates(df=df,taucol=1,zcol=3)
ML$GoodnessOfFit()
# estimate for rho=2.5
ML$Estimates(df,taucol=1,zcol=4)
ML$GoodnessOfFit()
# read experimental data
df<-OUPDataRead("Agric_NSW_SoilHealthHarden")
# estimate for nitrogen burn
ML$Estimates(df=df,taucol=1,zcol=2)
ML$GoodnessOfFit()
# read commodities data
df<-OUPDataRead("Finance_Commodities")
# estimate for West Texas intermediate
ML$Estimates(df=df,taucol=1,zcol=7)
ML$GoodnessOfFit()
# scaled brownian motion
ML$Estimates(rhor=0)
ML$LikelihoodRatioTest()
# stationary
ML$Estimates(rhor=99999)
ML$LikelihoodRatioTest()
# 95% upper bound on rho
ML$Estimates(rhor=0.680)
ML$LikelihoodRatioTest()
# mu irrelevant for scaled brownian motion
ML$Estimates(rhor=0,mur=99999)
ML$LikelihoodRatioTest()
ML$Estimates(rhor=0,mur=0)
ML$LikelihoodRatioTest()
# 95% lower bound on sigma
ML$Estimates(sigmar=0.533)
ML$LikelihoodRatioTest()
# infinite upper bound on sigma with stationary process
df<-OUPDataRead("Agric_NSW_SoilHealthHarden")
ML$Estimates(df=df,taucol=1,zcol=2)
ML$Estimates(sigmar=99999)
ML$LikelihoodRatioTest()
ML$Estimates(sigmar=-99999)
ML$LikelihoodRatioTest()
# 95% lower bound of sigma
ML$Estimates(sigmar=79.1)
ML$LikelihoodRatioTest()


