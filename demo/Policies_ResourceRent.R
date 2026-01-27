# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Policies_ResourceRent.R",echo=TRUE)', or
# Type 'demo(Policies_ResourceRent)'
# R6 object
A <- Analytical$new()
# Prohibition
A$set_plot_info(type=2)
A$set_x_stoch_args(phi=1)
# zero interest rate
A$Obligation(r=0)
# break even
A$Obligation(y=5)
# resource rent
A$Obligation(b=-10)
# positive interest rate
A$Obligation(y=0,r=0.05,b=0)
# break even
A$Obligation(y=5)
# resource rent
A$Obligation(b=-10)
# entry no discounting
A$DecisionThreshold(y=0,r=0,b=0)
A$DecisionThreshold(y=5)
A$DecisionThreshold(b=-10)
# entry with discounting
A$DecisionThreshold(y=0,r=0.05,b=0)
A$DecisionThreshold(y=5)
A$DecisionThreshold(b=-10)
