# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Estimation_MyData.R",echo=TRUE)', or
# Type 'demo(Estimation_MyData)'
# R6 object
ML <- MaximumLikelihood$new()
# read default data
df<-OUPDataRead()
# estimate
ML$Estimates(df)
# plot
ML$PlotEstimates()
