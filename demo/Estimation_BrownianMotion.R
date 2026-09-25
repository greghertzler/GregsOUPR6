# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_BrownianMotion.R",echo=TRUE)', or
# Type 'demo(Estimation_BrownianMotion)'
# R6 object
ML <- MaximumLikelihood$new()
# read simulated data
df<-OUPDataRead("OUP_BrownianMotion")
# test stationary increments e(s)
u<-ML$Estimates(df=df,taucol=1,zcol=2)
r<-ML$Estimates(rhor=800,mur=0,sigmar=40)
ML$LikelihoodRatioTest()
# test z(t)=z(s)+e(s)
u<-ML$Estimates(df=df,taucol=1,zcol=3)
r<-ML$Estimates(rhor=0)
ML$LikelihoodRatioTest()
# test sigma=1
r<-ML$Estimates(rhor=0,sigmar=1)
ML$LikelihoodRatioTest()
