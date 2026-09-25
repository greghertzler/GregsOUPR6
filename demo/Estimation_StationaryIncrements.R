# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_StationaryIncrements.R",echo=TRUE)', or
# Type 'demo(Estimation_StationaryIncrements)'
# R6 object
ML <- MaximumLikelihood$new()
# read simulated data
df<-OUPDataRead("OUP_StationaryIncrements")
# test dqrng::pcg64
u<-ML$Estimates(df=df,taucol=1,zcol=2)
r<-ML$Estimates(rhor=800,mur=0,sigmar=40)
ML$LikelihoodRatioTest()
# test std::mt19937
u<-ML$Estimates(df=df,taucol=1,zcol=3)
r<-ML$Estimates(rhor=800,mur=0,sigmar=40)
ML$LikelihoodRatioTest()
# testr sitmo::prng
u<-ML$Estimates(df,taucol=1,zcol=4)
r<-ML$Estimates(rhor=800,mur=0,sigmar=40)
ML$LikelihoodRatioTest()
# test rnorm
u<-ML$Estimates(df,taucol=1,zcol=5)
r<-ML$Estimates(rhor=800,mur=0,sigmar=40)
ML$LikelihoodRatioTest()
