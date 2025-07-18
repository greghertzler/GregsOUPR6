# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/ML_PlotTimeSeries.R",echo=TRUE)', or
# Type 'demo(ML_PlotTimeSeries)'
# R6 object
ML <- MaximumLikelihood$new()
# default
ML$PlotTimeSeries()
# customize plot
# set time series info
ML$set_timeseries_info(dataname="My Data",timename="week",statename="bushels and pecks")
ML$PlotTimeSeries(tbeg=22,tend=28)
# override time series info
ML$PlotTimeSeries(title="Title",xaxis="xAxis",yaxis="yAxis")
# revert to time series info
ML$PlotTimeSeries()
