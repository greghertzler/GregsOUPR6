# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Proof_DecisionThreshold.R",echo=TRUE)', or
# Type 'demo(Proof_DecisionThreshold)'
# R6 object
A <- Analytical$new()
A$set_oup_params(rho=0,mu=0,sigma=30)
A$set_x_stoch_args(s=seq(from=0,to=20,by=0.1),x=seq(from=-60,to=60,by=0.6))
# dark theme has opaque backgroound to view properly in Viewer
A$set_plot_info(theme="dark",opaque=1,walls=FALSE,floor=FALSE,labels=FALSE)
# obligation and exit
A$set_x_stoch_args(phi=-1)
A$PlotOptionEnvelope(title="Exit Option Envelope",type=6)
A$PlotObligation(title="Obligation with Exit Option Envelope",type=6)
# prohibition and entry
A$set_x_stoch_args(phi=1)
A$PlotOptionEnvelope(title="Entry Option Envelope",type=6)
A$PlotObligation(title="Prohibition with Entry Option Envelope",type=6)
