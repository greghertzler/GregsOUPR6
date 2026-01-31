# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Adoption_ExitOption.R",echo=TRUE)', or
# Type 'demo(Adoption_ExitOption)'
# R6 object
A <- Analytical$new()
A$set_plot_info(type=2)
# default
A$DecisionThreshold(phi=-1)
# sigma
A$DecisionThreshold(sigma=0)
A$DecisionThreshold(sigma=15)
A$DecisionThreshold(sigma=30)
# mu
A$DecisionThreshold(sigma=15)
A$DecisionThreshold(mu=0)
A$DecisionThreshold(mu=-15)
# rho
A$DecisionThreshold(rho=0.1,mu=15)
A$DecisionThreshold(rho=0.5)
A$DecisionThreshold(rho=1)
A$DecisionThreshold(rho=0.1,mu=-15)
A$DecisionThreshold(rho=0.5)
A$DecisionThreshold(rho=1)
