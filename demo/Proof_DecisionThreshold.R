# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Proof_DecisionThreshold.R",echo=TRUE)', or
# Type 'demo(Proof_DecisionThreshold)'
# R6 object
A <- Analytical$new()
A$set_oup_params(rho=0,mu=0,sigma=30)
A$set_x_stoch_args(s=seq(from=0,to=20,by=0.1),x=seq(from=-60,to=60,by=0.6))
# transparent walls, floor and no labels
A$set_plot_info(walls=FALSE,floor=FALSE,labels=FALSE)
# obligation and exit
A$set_x_stoch_args(phi=-1)
A$PlotOptionEnvelope(title="Exit Option Envelope",type=1)
A$PlotObligation(title="Obligation (click on Ohat for Exit Option Envelope)")
# prohibition and entry
A$set_x_stoch_args(phi=1)
A$PlotOptionEnvelope(title="Entry Option Envelope")
A$PlotObligation(title="Prohibition (click on Ohat for Entry Option Envelope)")
