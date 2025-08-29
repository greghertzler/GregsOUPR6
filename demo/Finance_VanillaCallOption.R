# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Finance_VanillaCallOption.R",echo=TRUE)', or
# Type 'demo(Finance_VanillaCallOption)'
# R6 object
A <- Analytical$new()
# 2D plots
A$set_plot_info(type=2)
# Brownian Motion
A$Option(rho=0,phi=1)
A$Option(mu=0)
# Ornstein-Uhlenbeck Process
A$default_read()
A$Option(phi=1)
# Stationary process
A$Option(rho=99)
