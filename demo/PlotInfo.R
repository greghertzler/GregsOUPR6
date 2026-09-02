# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/PlotInfo.R",echo=TRUE)', or
# Type 'demo(PlotInfo)'
# R6 objects
OUP <- OUProcess$new()
A <- OUP$get_Analytical()
FD <- OUP$get_FiniteDifference()
ML <- OUP$get_MaximumLikelihood()
MC <- OUP$get_MonteCarlo()
# default (type 0)
A$PlotPassageTimeProbability()
# learn about plot types
A$get_plot_types()
#set plot type the hard way
A$set_plot_type(type=1,group=6)
A$PlotPassageTimeProbability()
# vertical axis scales
A$set_plot_args(pmax=0.1,ptmax=0.6)
A$PlotDensity()
A$PlotPassageTimeDensity()
# change font family and size
A$set_plot_info(fontfamily="Arial",fontsize=16)
A$PlotDensity()
FD$PlotOption()
# to download as 'png', click the icon on the plot
ML$PlotTimeSeries()
# as 'svg' file (this only works for small 2D plots)
ML$set_plot_info(fileformat="svg")
A$PlotDensity()
# change file size
FD$set_plot_info(filewidth=140,fileheight=200,fileformat="png")
A$PlotDensity()
# change theme, transparent background
dklt <- A$get_plot_info()$plottheme$name
if(dklt == "dark") { A$set_plot_info(theme="light",opaque=0.0) }
if(dklt == "light") { A$set_plot_info(theme="dark",opaque=1.0) }
A$PlotDensity(type=1)
# transparent walls and floor
A$set_plot_info(walls=FALSE,floor=FALSE)
A$PlotDensity()
# no title or parameters
A$set_plot_info(labels=FALSE)
A$PlotDensity()
# title does not persist from plot to plot
A$PlotDensity(title="myTitle")
# axis labels do not persist from plot to plot
A$PlotDensity(title="myTitle",xaxis="myX",yaxis="myY",zaxis="myZ")
# cycle through plot types, 'd' for default, 'n' for next and 'p' for previous
A$PlotPassageTimePercentiles('d')
A$PlotPassageTimePercentiles('n')
A$PlotPassageTimePercentiles('n')
A$PlotPassageTimePercentiles('n')
A$PlotPassageTimePercentiles('n')
A$PlotPassageTimePercentiles('n')
A$PlotPassageTimePercentiles('n')


