# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/PlotInfo.R",echo=TRUE)', or
# Type 'demo(PlotInfo)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
MC <- OUP$get_MonteCarlo()
# default (type 3)
A$PlotPassageTimeProbability()
# type 4
A$set_plot_info(type=4)
A$PlotPassageTimeProbability()
# type 1
A$set_plot_info(type=1)
A$PlotPassageTimeProbability()
# type 2
A$set_plot_info(type=2)
A$PlotPassageTimeProbability()
# type 3
A$set_plot_info(type=3)
A$PlotPassageTimeProbability()
# vertical axis scales
A$set_plot_info(pmax=0.1,ptmax=0.6)
A$PlotDensity()
A$PlotPassageTimeDensity()
# change font family and size
A$set_plot_info(fontfamily="Arial",fontsize=16)
A$PlotDensity()
FD$PlotOption()
# download using the icon on the plot as 'png' file
ML$PlotTimeSeries()
# as 'svg' file
ML$set_plot_info(fileformat="svg")
A$PlotDensity()
# change file size
FD$set_plot_info(filewidth=140,fileheight=200)
A$PlotDensity()
# change theme, transparent background
A$set_plot_info(theme="light",opaque=0.0)
A$PlotDensity()
# transparent walls and floor
A$set_plot_info(walls=FALSE,floor=FALSE)
A$PlotDensity()
# no title or parameters
A$set_plot_info(labels=FALSE)
A$PlotDensity()
# title anyway does not persist from plot to plot
A$PlotDensity(title="myTitle")
# axis labels do not persist from plot to plot
A$PlotDensity(title="myTitle",xaxis="myX",yaxis="myY",zaxis="myZ")


