# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_PlotEstimates.R",echo=TRUE)', or
# Type 'demo(ML_PlotEstimates)'
# R6 object
ML <- MaximumLikelihood$new()
# default
ML$PlotEstimates()
# customize plot
ML$PlotEstimates(title="My Estimates",xaxis="xAxis",yaxis="yAxis",tbeg=22,tend=28)
# set time series info (previous title does not persist so set estimation)
ML$set_timeseries_info(dataname="My Data",timename="week",statename="bushels and pecks",estimation="My Estimates")
ML$PlotEstimates()

