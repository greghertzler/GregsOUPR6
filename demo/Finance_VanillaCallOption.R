# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_VanillaCallOption.R",echo=TRUE)', or
# Type 'demo(Finance_VanillaCallOption)'
# R6 object
A <- Analytical$new()
# automatically plot with calculations
A$set_flags(plotit=TRUE)
# Brownian Motion
A$Option(rho=0,phi=1)
# Normal Process (sigma^2/2rho=sigma^2)
A$Option(rho=0.5)
# Stationary process
A$Option(rho=99)
