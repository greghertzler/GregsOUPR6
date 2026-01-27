# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_Calibrate_Plots.R",echo=TRUE)', or
# Type 'demo(FD_Calibrate_Plots)'
library(plotly)
# R6 objects
FD <- FiniteDifference$new()
A <- Analytical$new()

# plotly preliminaries
plot_info <- FD$get_plot_info()
plot_colors <- FD$get_plot_colors()
font <- list(family=plot_info$plotfont$family,size=plot_info$plotfont$size,color=plot_colors$font)
file <- plot_info$plotfile
walls <- FALSE
floor <- FALSE
red <- plot_colors$red
ylw <- plot_colors$ylw
grn <- plot_colors$grn
cyn <- plot_colors$cyn
blu <- plot_colors$blu
mgn <- plot_colors$mgn
gry <- plot_colors$gry
background <- plot_colors$background
reverse <- plot_colors$reverse
bsml <- "<span style='font-size: 10pt;'>"
esml <- "</span>"
bsym <- "<span style='font-family:Symbol;'>"
esym <- "</span>"
lookup <- list(text="Calibrate",yref="container",y=0.95)
spy <- list(x=0,y=-2.5,z=0.5)
view <- list(camera=list(eye=spy),
  xaxis=list(title="<i>x</i>",color=font$color,linewidth=3,gridcolor=gry$d,gridwidth=2,backgroundcolor=gry$b,showbackground=walls,tickmode="auto",nticks=5),
  yaxis=list(title="<i>s</i>",color=font$color,linewidth=3,gridcolor=gry$d,gridwidth=2,backgroundcolor=gry$b,showbackground=walls,tickmode="auto",nticks=5,tickformat=" "),
  zaxis=list(title="err(<i>s,x</i>)",color=font$color,linewidth=3,gridcolor=gry$d,gridwidth=2,backgroundcolor=gry$c,showbackground=floor,range=c(-0.05,0.05),tickmode="auto",nticks=5),
  aspectratio=list(x=1,y=1,z=1))
gradient <- list(c(0,ylw$c),c(0.1,grn$c),c(0.2,cyn$c),c(0.3,blu$c),c(0.4,mgn$c),c(0.5,red$c),c(0.6,mgn$c),c(0.7,blu$c),c(0.8,cyn$c),c(0.9,grn$c),c(1,ylw$c))
rise <- list(x=0,y=-800,z=0)
shine <- list(ambient=0.7,diffuse=0.4,fresnel=0.2,roughness=2.0,specular=0.1)
hover <- list(bgcolor=mgn$e,font=list(color=mgn$a))
bar <- list(len=0.25,thickness=10,y=0.75,yanchor="top")
imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_Calibrate3DSurface")

# calibrate finite difference against analytical
s <- seq(from=10,to=0,by=-0.1)
x <- seq(from=-200,to=200,by=4)
m <- length(s)
n <- length(x)
coordinates <- matrix("",m,n)
FD$TerminalValue_Kinked(Vmax=10000,plotit=FALSE)
FDopt <- FD$Option(s=s,x=x,mu=-15,plotit=FALSE)[[1]]
Aopt <- A$Option(s=s,x=x,t=10,mu=-15,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# narrow x
x <- seq(from=-100,to=100,by=2)
FDopt <- FD$Option(x=x,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# tweak theta (can be more or less accurate)
FDopt <- FD$Option(theta=0.878,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# change skip (smaller is faster and maybe less accurate)
FDopt <- FD$Option(theta=0.5,skip=5,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# 4 times bigger sigma
FDopt <- FD$Option(sigma=60,skip=10,plotit=FALSE)[[1]]
Aopt <- A$Option(sigma=60,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# 2 times wider x
x <- seq(from=-200,to=200,by=4)
FDopt <- FD$Option(x=x,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,t=10,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# 2 times as many skips and x nodes
x <- seq(from=-200,to=200,by=2)
n <- length(x)
coordinates <- matrix("",m,n)
FDopt <- FD$Option(x=x,skip=20,plotit=FALSE)[[1]]
Aopt <- A$Option(x=x,t=10,plotit=FALSE)[[1]]
err <- FDopt-Aopt
message(paste("Max over: ",max(err),"  Max under: ",min(err)))
for(i in 1:m)
{
  for(j in 1:n)
  {
    coordinates[i,j] <- paste(sep="","err(<i>s,x</i>)=",format(err[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
  }
}
oup_params <- FD$get_oup_params()
x_stoch_args <- FD$get_x_stoch_args()
rho <- oup_params[[1]]
mu <- oup_params[[2]]
sigma <- oup_params[[3]]
r <- x_stoch_args[[4]]
theta <- x_stoch_args[[6]]
skip <- x_stoch_args[[7]]
ds <- x_stoch_args[[8]]
dx <- x_stoch_args[[9]]
syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0)
fig <- plot_ly() %>%
  add_trace(.,type="surface",x=x,y=s,z=err,cmid=0,colorbar=bar,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
  config(.,toImageButtonOptions=imageoptions) %>%
  layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
print(fig)

# check option envelope
FDenv <- FD$OptionEnvelope(plotit=FALSE)[[1]]
Aenv <- A$OptionEnvelope(plotit=FALSE)[[1]]
err <- FDenv-Aenv
message(paste("Max over: ",max(err),"  Max under: ",min(err)))

# check decision threshold
FDdec <- FD$DecisionThreshold(plotit=FALSE)
Adec <- A$DecisionThreshold(plotit=FALSE)
errk <- FDdec[[1]]-Adec[[1]]
errOhat <- FDdec[[2]]-Adec[[2]]
message(paste("k difference: ",errk,"  Ohat difference: ",errOhat))
