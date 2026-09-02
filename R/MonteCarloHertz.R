library(R6)
library(Rcpp)
library(plotly)
library(stringr)
library(clipr)

# roxygen ----
#' @title R6 class implementing Monte Carlo simulation of the Ornstein-Uhlenbeck Process
#'
#' @description
#' Monte Carlo simulations of Forward Paths, Bounded Paths, Backward Paths,
#'  Probabilities, Options and Passage Times.  Forward, Backward and Bounded Paths
#'  are simulated by a 4th order Runge-Kutta method and by using the stochastic
#'  integral equation.  Probabilities, Options and Passage Times are various ways
#'  of binning and counting the Paths.
#'
#' @details # Methods:
#'     y stochastic
#'       ForwardPaths
#'       Mean
#'       Variance
#'       Density
#'       Probability
#'       DoubleIntegral
#'     x stochastic
#'       BackwardPaths
#'       Option
#'     t stochastic
#'       ForwardPaths
#'       BoundedPaths
#'       VisitingTimeModeMedianMean
#'       VisitingTimePercentiles
#'       VisitingTimeDensity
#'       VisitingTimeProbability
#'       FirstPassageTimeModeMedianMean
#'       FirstPassageTimePercentiles
#'       FirstPassageTimeDensity
#'       FirstPassageTimeProbability
#'
#' @details # Plots:
#'       PlotForwardPaths
#'       PlotBackwardPaths
#'       PlotBoundedPaths
#'       PlotMean
#'       PlotVariance
#'       PlotDensity
#'       PlotProbability
#'       PlotDoubleIntegral
#'       PlotOption
#'       PlotVisitingTimeModeMedianMean
#'       PlotVisitingTimePercentiles
#'       PlotVisitingTimeDensity
#'       PlotVisitingTimeProbability
#'       PlotFirstPassageTimeModeMedianMean
#'       PlotFirstPassageTimePercentiles
#'       PlotFirstPassageTimeDensity
#'       PlotFirstPassageTimeProbability
#'
#' @details # Arguments of functions:
#'       All arguments are optional in all functions.
#'     OUP parameters
#'       rho:    rate parameter 0<=rho<inf
#'       mu:     location parameter -inf<mu<inf
#'       sigma:  scale parameter -inf<sigma<inf
#'     y stochastic
#'       t:      vector of times s<=t<inf
#'       x:      initial state -inf<x<inf
#'       psi:    <=0 for integral -inf to y,
#'                >0 for integral y to inf
#'     x stochastic
#'       s:      vector of times -inf<s<t
#'       y:      terminal state -inf<y<inf
#'       r:      discount rate -inf<r<inf
#'       phi:    <=0 for exit option,
#'                >0 for entry option
#'     t stochastic
#'       t:      vector of times s<=t<inf
#'       k:      decision threshold -inf<k<int
#'       x:      initial state -inf<x<inf
#'       Ppct:   probability for a percentile 0.01<=Ppct<=0.99
#'     path
#'       paths:  number of paths 1<paths<1,000,000
#'       skip:   subdivide time interval but report at times t 1<=skip<=50
#'       seed:   seed for random number generators -inf<seed<inf
#'       method: 4 for 4th order Runge-Kutta, otherwise integral equation
#'
#' @details # Usage:
#' The MonteCarlo object must first be instantiated before its methods are called.
#'  There are two ways.  The first way instantiates the OUProcess object and
#'  calls a function to get a pointer:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       FD <- OUP$get_FiniteDifference()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'
#' The MonteCarlo object will coordinate arguments to functions with the
#'  Analytial, FiniteDifference and MaximumLikelihood objects.  The second
#'  way instantiates the MonteCarlo object by itself with no coordination:
#'
#'       MC <- MonteCarlol$new()
#'
#' Once the object is instantiated, its methods can be called, to simulate
#'  100,000 Forward Paths, for example:
#'
#'       MC$ForwardPaths(paths=100000)
#'
#' The plot methods can be used to customize the plots, with a title, for example:
#'
#'       A$PlotForwardPaths(title="My Paths Forward")
#'
#' An attempt to plot 100,000 paths would choke the computer, so there are tricks.
#'  One is to select a hundred or so paths for the plot.  Another is to plot heat
#'   maps, just like in a weather report.
#'
#' Other functions and methods are called in the same way.  To see all the
#'  possibilities, check out the demos below.
#'
#' @details # Demos:
#' Demonstration scripts are in files in the 'demo' directory. The most
#'  convenient way to list and run demos are the commands:
#'
#'       OUPDemoList()
#'       OUPDemoRun(<number of demo in list>)
#'
#' Entering the demos by number in the list saves typing.
#'
#' @details # Discussion:
#' Monte Carlo simulation of the Ornstein-Uhlenbeck Process can be done with
#'  either of two methods:  numerically integrating the stochastic differential
#'  equation, or calculating the stochastic integral equation.  The stochastic
#'  differential equation is shocked by a Wiener Process, simulated as
#'  sigma * dt^0.5 * epsilon, where sigma * dt^0.5 is the square root of the
#'  instantaneous variance and epsilon are draws from a standard normal
#'  density. The stochastic integral equation is shocked by the integral of
#'  the Wiener Process, or H * epsilon, where H is the square-root of the
#'  variance over a longer time interval.
#'
#' Numerically integrating the stochastic differential equation uses the Euler,
#'  Marayuma or Runge-Kutta schemes. The Euler and Marayuma schemes are first
#'  and third-order schemes and less accurate than the fourth-order Runge-Kutta
#'  scheme.
#'
#' The fourth-order Runge-Kutta scheme is standard practice.  In tests, shocked
#'  by the same draws from a standard normal density, the paths from the stochastic
#'  integral equation and the fourth-order Runge-Kutta scheme were the same to
#'  within five or six significant digits. To compare the fourth-order Runge-
#'  Kutta scheme with the integral equation:
#'
#'       MC <- MonteCarlo$new()
#'       rk <- MC$ForwardPaths(paths=100000,skip=10,method=4)[[1]]
#'       ie <- MC$ForwardPaths(method=1)[[1]]
#'       dif <- rk-ie
#'       max(dif)
#'       min(dif)
#'       sum(dif)
#'
#' A Forward Path starts from the backward state at the backward time and goes
#'  forward. A single path, sampled from all possible paths, is a Sample Path.
#'  Just like the flea trying to understand the elephant, a sample Path is enough
#'  for Maximum Likelihood Estimation of the Ornstein-Uhlenbeck Process. An
#'  ensemble of paths can be counted to approximate Transition Densities and
#'  Probabilities and Visiting Time Densities and Probabilities. The larger the
#'  ensemble, the better the approximations.
#'
#' Forward Paths are continuous in the state and time which makes them very hard
#'  to count. Instead, Forward Paths are treated as discrete and put into bins.
#'  Counting the number of Forward Paths in each bin and dividing by the total
#'  number of paths approximates Transition Densities. Summing the Transition
#'  Densities approximates Transition Probabilities. Summing again approximates
#'  Double Integrals.  Double Integrals are a curiosity. If time runs backwards,
#'  they become Options.
#'
#' A Forward Path begins from a known state and travels forward into an uncertain
#'  future.  A Backward Path ends with a known state and trudges backward into an
#'  uncertain past. Turning around and travelling back to the future resolves the
#'  uncertainty over time. An example is a Bayesian analysis which begins with a
#'  Diffuse Prior and ends with certainty. Another example is an Option. Simulating
#'  and counting Backward Paths approximates Prior Densities, Prior Probabilities
#'  and Options.  To compare Monte Carlo and Analytical Options:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'       ML$Estimates()
#'       ao <- A$Option()[[1]]
#'       mco <- MC$Option(paths=1000000)[[1]]
#'       dif <- ao-mco
#'       max(dif)
#'       min(dif)
#'       sum(dif)
#'
#' If we count at right angles--in the time direction instead of the state
#'  direction--Forward Paths become Visiting Time Densities and Probabilities.
#'  Bounded Paths become First Passage Time Densities and Probabilities.
#'  Monte Carlo and Analytical Visiting and First Passage Times can also be
#'  compared:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       MC <- OUP$get_MonteCarlo()
#'       ad <- A$PassageTimeDensity(omega=0)[[1]]
#'       mcd <- MC$VisitingTimeDensity(paths=1000000)[[1]]
#'       dif <- ad-mcd
#'       max(dif)
#'       min(dif)
#'       sum(dif)
#'
#' Of course, the question is, 'Why bother?' We have Analytical formulas to do
#'  the counting.  One reason is to explain the formulas.  First Passage Times
#'  make more sense if you plot Bounded Paths and count the number of paths that
#'  have crossed the threshold.  Even in prestigious journal articles, the first
#'  and, possibly, only plot will be a Monte Carlo simulation.
#'
#' Another reason is to validate the formulas.  Although an Analytical formula
#'  may calculate a thousand times faster than a Monte Carlo simulation, arriving
#'  at approximately the same answer both ways is reassuring.

# class ----
MonteCarlo <- R6::R6Class("MonteCarlo",
  portable = FALSE,
  cloneable = FALSE,
  # portable = TRUE,
  # cloneable = TRUE,
#' @import plotly
#' @import stringr
#' @importFrom clipr clipr_available write_clip
#' @export
  #public members ----
  public = list(
    # constructor ----
    #' @description
    #' Create a MonteCarlo object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new MonteCarlo object
    initialize = function(OUP=NULL)
    {
      # pointer to container object ----
      if(!is.null(OUP) && class(OUP)[[1]] == "OUProcess") { private$OUP <- OUP }
      # arguments ----
      private$oup_params <- list(rho=0.5,mu=15,sigma=15)
      xyseq <- seq(from=-30,to=30,by=0.6)
      private$y_stoch_args <- list(t=seq(from=0,to=10,by=0.1),y=xyseq,x=-15,psi=-1)
      private$x_stoch_args <- list(s=seq(from=10,to=0,by=-0.1),x=xyseq,y=0,r=0.05,phi=-1)
      private$t_stoch_args <- list(t=seq(from=0,to=10,by=0.1),k=20,x=-15,omega=1,Ppct=0.75)
      private$path_args <- list(paths=100,skip=1,seed=99999,method=1)
      private$plot_args <- list(pmax=0.06,ptmax=0.6,first=1,last=10,zbeg=-30,zend=30)
      private$syncyxt <- 3
      private$forwardyt <- 3
      # undo ----
      private$undo_args <- list(list(oup_params=private$oup_params,y_stoch_args=private$y_stoch_args,x_stoch_args=private$x_stoch_args,t_stoch_args=private$t_stoch_args,path_args=private$path_args,plot_args=private$plot_args))
      private$undoIx <- 1
      # plot info ----
      plotfont <- list(family="Cambria",size=14)
      plotfile <- list(format="png",width=640,height=480)
      plottheme <- list(name="light",opaque=0.0)
      if(Sys.getenv("RSTUDIO") == "1")
      {
        if(rstudioapi::getThemeInfo()$dark) { plottheme <- list(name="dark",opaque=1.0) }
      }
      plot3D <- list(walls=TRUE,floor=TRUE)
      private$plot_info <- list(plotfont=plotfont,plotfile=plotfile,plottheme=plottheme,plot3D=plot3D,plotlabels=TRUE)
      # flags ----
      private$flags <- list(plotit=FALSE,copyit=FALSE)
      # plot globals ----
      private$plot_colors <- private$rainbow(plottheme$name,plottheme$opaque)
      private$modebar_2D <- list(list("zoom2d", "pan2d", "resetScale2d"),list("toImage"))
      private$modebar_3D <- list(list("zoom3d", "pan3d", "orbitRotation", "tableRotation"),list("resetCameraDefault3d", "resetCameraLastSave3d"),list("hoverClosest3d"),list("toImage"))
      private$plot_types <- c(rep(0,5))
      # ignore warnings from plotly
      options(warn=-1)
    },
    # public set methods ----
    #' @description
    #' Set OUP parameters
    #' @param rho   rate parameter 0<=rho<inf
    #' @param mu    location parameter -inf<mu<inf
    #' @param sigma scale parameter -inf<sigma<inf
    #' @param who   object id of sender
    #' @return list(rho,mu,sigma)
    set_oup_params = function(rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"MC") }
      if(!is.null(rho))
      {
        sca <- private$extract_scalar(rho)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            sca <- 0.0
            message("negative rho set to zero.")
          }
          if(sca != private$oup_params$rho)
          {
            private$oup_params$rho <- sca
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("rho not set.")}
      }
      if(!is.null(mu))
      {
        sca <- private$extract_scalar(mu)
        if(!is.null(sca))
        {
          if(sca != private$oup_params$mu)
          {
            private$oup_params$mu <- sca
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("mu not set.")}
      }
      if(!is.null(sigma))
      {
        sca <- private$extract_scalar(sigma)
        if(!is.null(sca))
        {
          if(sca != private$oup_params$sigma)
          {
            private$oup_params$sigma <- sca
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("sigma not set.")}
      }
      return(private$oup_params)
    },
    #' @description
    #' Set arguments for y as a stochastic state
    #' @param t   vector of m times -inf<t<inf
    #' @param y   vector of n states -inf<y<inf
    #' @param x   initial state -inf<x<inf
    #' @param psi <=0 for integral -inf to y, >0 for integral y to inf
    #' @param who object id of sender
    #' @return list(t,y,x,psi)
    set_y_stoch_args = function(t=NULL,y=NULL,x=NULL,psi=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_y_stoch_args(t,y,x,psi,"MC") }
      if(!is.null(t))
      {
        vec <- private$extract_vector(t,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$y_stoch_args$t))
          {
            if(length(vec) < 102)
            {
              private$y_stoch_args$t <- vec
              private$ystdnorm <- NULL
              private$yforward <- NULL
              private$G <- NULL
              private$H2 <- NULL
              private$p <- NULL
              private$Pneg <- NULL
              private$Ppos <- NULL
              private$PPneg <- NULL
              private$PPpos <- NULL
            }
            else { message("101 times at most.  t not set.") }
          }
        }
        else { message("t not set.") }
      }
      if(!is.null(y))
      {
        vec <- private$extract_vector(y,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$y_stoch_args$y))
          {
            private$y_stoch_args$y <- vec
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
          }
        }
        else { message("y not set.") }
      }
      if(!is.null(x))
      {
        sca <- private$extract_scalar(x)
        if(!is.null(sca))
        {
          if(sca != private$y_stoch_args$x)
          {
            private$y_stoch_args$x <- sca
            private$yforward <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
          }
        }
        else { message("x not set.") }
      }
      if(!is.null(psi))
      {
        sca <- private$extract_scalar(psi)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            if(sca != -1)
            {
              sca <- -1
              message("psi set to -1.")
            }
          }
          else
          {
            if(sca != 1)
            {
              sca <- 1
              message("psi set to 1.")
            }
          }
          if(sca != private$y_stoch_args$psi)
          {
            private$y_stoch_args$psi <- sca
          }
        }
        else { message("psi not set.") }
      }
      private$syncyxt <- 1
      private$forwardyt <- 1
      return(private$y_stoch_args)
    },
    #' @description
    #' Set arguments for x as a stochastic state
    #' @param s   vector of m times -inf<s<t
    #' @param x   vector of n states -inf<x<inf
    #' @param y   terminal state -inf<y<inf
    #' @param r   discount rate -inf<r<inf
    #' @param phi <=0 for exit option, >0 for entry option
    #' @param who object id of sender
    #' @return list(s,x,y,r,phi)
    set_x_stoch_args = function(s=NULL,x=NULL,y=NULL,r=NULL,phi=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,y,r,phi,"MC") }
      if(!is.null(s))
      {
        vec <- private$extract_vector(s,-1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$x_stoch_args$s))
          {
            if(length(vec) < 102)
            {
              private$x_stoch_args$s <- vec
              private$xstdnorm <- NULL
              private$xbackward <- NULL
              private$o <- NULL
              private$Oneg <- NULL
              private$Opos <- NULL
              private$OOneg <- NULL
              private$OOpos <- NULL
            }
            else { message("101 times at most.  s not set.") }
          }
        }
        else { message("s not set.") }
      }
      if(!is.null(x))
      {
        vec <- private$extract_vector(x,1)
        if(!is.null(vec))
        {
          n <- length(vec)
          if(n > 100 || is.null(private$OUP))
          {
            if(!private$vecs_equal(vec,private$x_stoch_args$x))
            {
              private$x_stoch_args$x <- vec
              private$xstdnorm <- NULL
              private$xbackward <- NULL
              private$o <- NULL
              private$Oneg <- NULL
              private$Opos <- NULL
              private$OOneg <- NULL
              private$OOpos <- NULL
            }
          }
          else
          {
            message("x vector must have at least 101 elements")
            message("x not set.")
          }
        }
        else { message("x not set.") }
      }
      if(!is.null(y))
      {
        sca <- private$extract_scalar(y)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$y)
          {
            private$x_stoch_args$y <- sca
            private$xbackward <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
          }
        }
        else { message("y not set.") }
      }
      if(!is.null(r))
      {
        sca <- private$extract_scalar(r)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            message("negative r set to zero.")
            sca <- 0.0
          }
          if(sca != private$x_stoch_args$r)
          {
            private$x_stoch_args$r <- sca
            private$OOneg <- NULL
            private$OOpos <- NULL
          }
        }
        else { message("r not set.") }
      }
      if(!is.null(phi))
      {
        sca <- private$extract_scalar(phi)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            if(sca != -1)
            {
              sca <- -1
              message("phi set to -1.")
            }
          }
          else if(sca > 0)
          {
            if(sca != 1)
            {
              sca <- 1
              message("phi set to 1.")
            }
          }
          if(sca != private$x_stoch_args$phi)
          {
            private$x_stoch_args$phi <- sca
          }
        }
        else { message("phi not set.") }
      }
      private$syncyxt <- 2
      return(private$x_stoch_args)
    },
    #' @description
    #' Set t stochastic arguments
    #' @param t     vector of m times -inf<t<inf
    #' @param k     threshold -inf<k<inf
    #' @param x     initial state -inf<x<inf
    #' @param omega degree of irreversibility 0<=omega<=1
    #' @param Ppct  passage time probability for a percentile  0.01<=Ppct<=0.99
    #' @param who   object id of sender
    #' @return list(t,k,x)
    set_t_stoch_args = function(t=NULL,k=NULL,x=NULL,omega=NULL,Ppct=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_t_stoch_args(t,k,x,omega,Ppct,"MC") }
      if(!is.null(t))
      {
        vec <- private$extract_vector(t,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$t_stoch_args$t))
          {
            if(length(vec) < 102)
            {
              private$t_stoch_args$t <- vec
              private$tstdnorm <- NULL
              private$tforward <- NULL
              private$tbounded <- NULL
              private$tfpt <- NULL
              private$vtmmm <- NULL
              private$vtpct <- NULL
              private$pv <- NULL
              private$Pv <- NULL
              private$fheat <- NULL
              private$fz <- NULL
              private$fptmmm <- NULL
              private$fptpct <- NULL
              private$pf <- NULL
              private$Pf <- NULL
              private$bheat <- NULL
              private$bz <- NULL
            }
            else { message("101 times at most.  t not set.") }
          }
        }
        else { message("t not set.") }
      }
      if(!is.null(k))
      {
        sca <- private$extract_scalar(k)
        if(!is.null(sca))
        {
          if(sca != private$t_stoch_args$k)
          {
            private$t_stoch_args$k <- sca
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("k not set.") }
      }
      if(!is.null(x))
      {
        sca <- private$extract_scalar(x)
        if(!is.null(sca))
        {
          if(sca != private$t_stoch_args$x)
          {
            private$t_stoch_args$x <- sca
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("x not set.") }
      }
      if(!is.null(Ppct))
      {
        sca <- private$extract_scalar(Ppct)
        if(!is.null(sca))
        {
          if(sca < 0.01)
          {
            sca <- 0.01
            message("Ppct has been set to 0.01.")
          }
          else if(sca > 0.99)
          {
            sca <- 0.99
            message("Ppct has been set to 0.99.")
          }
          if(sca != private$t_stoch_args$Ppct)
          {
            private$t_stoch_args$Ppct <- sca
            private$vtpct <- NULL
            private$fptpct <- NULL
          }
        }
        else { message("Ppct not set.") }
      }
      private$syncyxt <- 3
      private$forwardyt <- 3
      return(private$t_stoch_args)
    },
    #' @description
    #' Set path arguments
    #' @param paths  number of paths 1<paths<1,000,000
    #' @param skip   subdivide time intervals but report at times t 1<=skip<=50
    #' @param seed   seed for random number generators -inf<seed<inf
    #' @param method 4 for 4th order Runge-Kutta, otherwise integral equation
    #' @return list(paths,skip,seed,method)
    set_path_args = function(paths=NULL,skip=NULL,seed=NULL,method=NULL)
    {
      if(!is.null(paths))
      {
        sca <- private$extract_integer(paths)
        if(!is.null(sca))
        {
          if(sca != private$path_args$paths)
          {
            if(sca < 1 )
            {
              sca <- 1
              message("paths set to 1.")
            }
            else if (sca > 1000000)
            {
              sca <- 1000000
              message("paths set to 1,000,000.")
            }
            private$path_args$paths <- sca
            private$ystdnorm <- NULL
            private$xstdnorm <- NULL
            private$tstdnorm <- NULL
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("paths not set.") }
      }
      if(!is.null(skip))
      {
        sca <- private$extract_integer(skip)
        if(!is.null(sca))
        {
          if(sca != private$path_args$skip)
          {
            if(sca < 1 )
            {
              sca <- 1
              message("skip set to 1.")
            }
            else if (sca > 50)
            {
              sca <- 50
              message(paste(sep="","skip set to 50."))
            }
            private$path_args$skip <- sca
            private$ystdnorm <- NULL
            private$xstdnorm <- NULL
            private$tstdnorm <- NULL
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("skip not set.") }
      }
      if(!is.null(seed))
      {
        sca <- private$extract_scalar(seed)
        if(!is.null(sca))
        {
          if(sca != private$path_args$seed)
          {
            private$path_args$seed <- sca
            private$ystdnorm <- NULL
            private$xstdnorm <- NULL
            private$tstdnorm <- NULL
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("seed not set.") }
      }
      if(!is.null(method))
      {
        sca <- private$extract_scalar(method)
        if(!is.null(sca))
        {
          if(sca != private$path_args$method)
          {
            private$path_args$method <- sca
            private$yforward <- NULL
            private$xbackward <- NULL
            private$tforward <- NULL
            private$tbounded <- NULL
            private$tfpt <- NULL
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$o <- NULL
            private$Oneg <- NULL
            private$Opos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$vtmmm <- NULL
            private$vtpct <- NULL
            private$pv <- NULL
            private$Pv <- NULL
            private$fheat <- NULL
            private$fz <- NULL
            private$fptmmm <- NULL
            private$fptpct <- NULL
            private$pf <- NULL
            private$Pf <- NULL
            private$bheat <- NULL
            private$bz <- NULL
          }
        }
        else { message("method not set.") }
      }
      return(private$path_args)
    },
    #' @description
    #' Set plot arguments
    #' @param pmax  maximum transition density
    #' @param ptmax maximum visiting time and first passage time densities
    #' @param first first path to plot
    #' @param last  last path to plot
    #' @param zbeg  begin value for state axis
    #' @param zend  end value for state axis
    #' @param who   object id of caller
    #' @return list(pmax,ptmax,first,last,zbeg,zend)
    set_plot_args = function(pmax=NULL,ptmax=NULL,first=NULL,last=NULL,zbeg=NULL,zend=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_args(pmax,ptmax,"MC") }
      if(!is.null(pmax))
      {
        if(is.numeric(pmax))
        {
          if(is.infinite(pmax) || is.nan(pmax)) { private$plot_args$pmax <- NaN }
          else
          {
            sca <- private$extract_scalar(pmax)
            if(sca <= 0) { sca <- NaN }
            private$plot_args$pmax <- sca
          }
        }
        else { message("pmax not set.") }
      }
      if(!is.null(ptmax))
      {
        if(is.numeric(ptmax))
        {
          if(is.infinite(ptmax) || is.nan(ptmax)) { private$plot_args$ptmax <- NaN }
          else
          {
            sca <- private$extract_scalar(ptmax)
            if(sca <= 0) { sca <- NaN }
            private$plot_args$ptmax <- sca
          }
        }
        else { message("ptmax not set.") }
      }
      if(!is.null(first) && !is.null(last))
      {
        p1 <- private$extract_integer(first)
        pn <- private$extract_integer(last)
        if(!is.null(p1) && !is.null(pn))
        {
          hit <- FALSE
          n <- private$path_args$paths
          if(p1 > pn)
          {
            hit <- TRUE
            sca <- p1
            p1 <- pn
            pn <- sca
          }
          if(p1 < 1)
          {
            hit <- TRUE
            p1 <- 1
          }
          else if(p1 > n)
          {
            hit <- TRUE
            p1 <- n
          }
          if(pn < 1)
          {
            hit <- TRUE
            pn <- 1
          }
          else if(pn > n)
          {
            hit <- TRUE
            pn <- n
          }
          if(pn > p1+99)
          {
            hit <- TRUE
            pn <- p1+99
          }
          if(hit) { message(paste(sep="","first:last set to ",p1,":",pn,".")) }
          private$plot_args$first <- p1
          private$plot_args$last <- pn
        }
        else { message("first and last not set.") }
      }
      else if(!is.null(first))
      {
        p1 <- private$extract_integer(first)
        if(!is.null(p1))
        {
          pn <- private$plot_args$last
          if(p1 > pn)
          {
            p1 <- pn
            message(paste(sep="","first set to ",p1,"."))
          }
          else if(p1 < pn-99 && pn-99 > 0)
          {
            p1 <- pn-99
            message(paste(sep="","first set to ",p1,"."))
          }
          if(p1 < 1)
          {
            p1 <- 1
            message("first set to 1.")
          }
          private$plot_args$first <- p1
        }
        else { message("first not set.") }
      }
      else if(!is.null(last))
      {
        pn <- private$extract_integer(last)
        if(!is.null(pn))
        {
          n <- private$path_args$paths
          p1 <- private$plot_args$first
          if(pn < p1)
          {
            pn <- p1
            message(paste(sep="","last set to ",pn,"."))
          }
          else if(pn > p1+99 && p1+99 < n+1)
          {
            pn <- p1+99
            message(paste(sep="","last set to ",pn,"."))
          }
          if(pn > n)
          {
            pn <- n
            message(paste(sep="","last set to ",pn,"."))
          }
          private$plot_args$last <- pn
        }
        else { message("last not set.") }
      }
      if(!is.null(zbeg) && !is.null(zend))
      {
        beg <- private$extract_scalar(zbeg)
        end <- private$extract_scalar(zend)
        if(!is.null(beg) && !is.null(end))
        {
          if(beg > end)
          {
            sca <- beg
            beg <- end
            end <- sca
            message(paste(sep="","zbeg:zend set to ",beg,":",end,"."))
          }
          private$plot_args$zbeg <- beg
          private$plot_args$zend <- end
        }
        if(!is.null(beg)) { private$plot_args$zbeg <- beg }
        else
        {
          private$plot_args$zbeg <- -Inf
          message(paste(sep="","zbeg set to -Inf."))
        }
        if(!is.null(end)) { private$plot_args$zend <- end }
        else
        {
          private$plot_args$zend <- Inf
          message(paste(sep="","zend set to Inf."))
        }
      }
      else if(!is.null(zbeg))
      {
        zbeg <- private$extract_scalar(zbeg)
        if(!is.null(zbeg))
        {
          zend <- private$plot_args$zend
          if(zbeg > zend)
          {
            private$plot_args$zbeg <- zend
            message(paste(sep="","zbeg set to ",zend,"."))
          }
          else{ private$plot_args$zbeg <- zbeg }
        }
        else
        {
          private$plot_args$zbeg <- -Inf
          message(paste(sep="","zbeg set to -Inf."))
        }
      }
      else if(!is.null(zend))
      {
        zend <- private$extract_scalar(zend)
        if(!is.null(zend))
        {
          zbeg <- private$plot_args$zbeg
          if(zend < zbeg)
          {
            private$plot_args$zend <- zbeg
            message(paste(sep="","zend set to ",zbeg,"."))
          }
          else{ private$plot_args$zend <- zend }
        }
        else
        {
          private$plot_args$zend <- Inf
          message(paste(sep="","zend set to Inf."))
        }
      }
      return(private$plot_args)
    },
    #' @description
    #' Set information for plotting
    #' @param fontfamily font family for plot labels
    #' @param fontsize   font size for plot labels
    #' @param fileformat 'png' or 'svg'
    #' @param filewidth  pixel width of 2D plot, pixel width and height of 3D plot
    #' @param fileheight pixel height of 2D plot
    #' @param theme      'light' or 'dark'
    #' @param opaque     transparent to opaque background 0.0<=opaque<=1.0
    #' @param walls      3D walls TRUE or FALSE
    #' @param floor      3D floor TRUE or FALSE
    #' @param labels     title and parameters TRUE or FALSE
    #' @param who        object id of sender
    #' @return list(font,file,theme,3D,labels)
    set_plot_info = function(fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,"MC") }
      if(!is.null(fontfamily))
      {
        chr <- private$extract_character(fontfamily)
        if(!is.null(chr)) { private$plot_info$plotfont$family <- chr }
        else { message("fontfamily not set.") }
      }
      if(!is.null(fontsize))
      {
        sca <- private$extract_scalar(fontsize)
        if(!is.null(sca)) { private$plot_info$plotfont$size <- sca }
        else { message("fontsize not set.") }
      }
      if(!is.null(fileformat))
      {
        chr <- private$extract_character(fileformat)
        if(!is.null(chr))
        {
          if(chr == "png" || chr == "svg") { private$plot_info$plotfile$format <- chr }
          else
          {
            message("fileformat must be 'png' or 'svg'.")
            message("fileformat not set.")
          }
        }
        else { message("fileformat not set.") }
      }
      if(!is.null(filewidth))
      {
        sca <- private$extract_scalar(filewidth)
        if(!is.null(sca)) { private$plot_info$plotfile$width <- sca }
        else { message("filewidth not set.") }
      }
      if(!is.null(fileheight))
      {
        sca <- private$extract_scalar(fileheight)
        if(!is.null(sca)) { private$plot_info$plotfile$height <- sca }
        else { message("fileheight not set.") }
      }
      if(!is.null(theme) || !is.null(opaque))
      {
        if(!is.null(theme))
        {
          chr <- private$extract_character(theme)
          if(!is.null(chr))
          {
            if(chr == "light" || chr == "dark") { private$plot_info$plottheme$name <- chr }
            else
            {
              message("theme not set.")
              message("available themes are: 'light' and 'dark'.")
            }
          }
          else { message("theme not set.") }
        }
        if(!is.null(opaque))
        {
          sca <- private$extract_scalar(opaque)
          if(!is.null(sca))
          {
            if(sca < 0.0)
            {
              sca = 0.0
              message("opaque set to 0.0.")
            }
            else if(sca > 1.0)
            {
              sca = 1.0
              message("opaque set to 1.0.")
            }
            private$plot_info$plottheme$opaque <- sca
          }
          else { message("opaque not set.") }
        }
        private$plot_colors <- private$rainbow(private$plot_info$plottheme$name,private$plot_info$plottheme$opaque)
      }
      if(!is.null(walls))
      {
        bool <- private$extract_boolean(walls)
        if(!is.null(bool)) { private$plot_info$plot3D$walls <- bool  }
        else { message("walls not set.") }
      }
      if(!is.null(floor))
      {
        bool <- private$extract_boolean(floor)
        if(!is.null(bool)) { private$plot_info$plot3D$floor <- bool  }
        else { message("floor not set.") }
      }
      if(!is.null(labels))
      {
        bool <- private$extract_boolean(labels)
        if(!is.null(bool)) { private$plot_info$plotlabels <- bool  }
        else { message("labels not set.") }
      }
      return(private$plot_info)
    },
    #' @description
    #' Set plot types by group
    #' @param type a number identifying the type for a group
    #' @param group identifier for groups of similar plots
    #' @return list(type,group)
    set_plot_type = function(type=NULL,group=NULL)
    {
      group <- private$extract_integer(group)
      if(is.null(group)) { group <- 1 }
      if(!is.null(type))
      {
        if(group == 2) # Mean, Variance
        {
          least <- -1
          most <- 0
        }
        else if(group == 3) # Density, Probability, Double Integral, Option
        {
          least <- 0
          most <- 1
        }
        else if(group == 4) # Visiting Time Mode Median Mean, Visiting Time Percentiles, First Passage Time Mode Median Mean, First Passage Time Percentiles
        {
          least <- -1
          most <- 0
        }
        else if(group ==5) # Visiting Time Density, Visiting Time Probability, First Passage Time Density, First Passage Time Probability
        {
          least <- 0
          most <- 1
        }
        else # Forward Paths, Backward Paths, Bounded Paths
        {
          group <-1
          least <- -3
          most <- 0
        }
        if(is.character(type))
        {
          chr <- str_sub(type,start=1,end=1)
          if(chr == "n" || chr == "N") {  private$plot_types[group] <- private$plot_types[group]+1 }
          else if(chr == "p" || chr == "P") { private$plot_types[group] <- private$plot_types[group]-1 }
          else if(chr == "d" || chr == "D") { private$plot_types[group] <- 0 }
        }
        else if(is.numeric(type)) { private$plot_types[group] <- type }
        if(private$plot_types[group] > most) { private$plot_types[group] <- least }
        else if(private$plot_types[group] < least) { private$plot_types[group] <- most }
      }
      return(list(type=private$plot_types[group],group=group))
    },
    #' @description
    #' Set flags for plotting and copying
    #' @param plotit automatic plot after calculation TRUE or FALSE
    #' @param copyit copy to clipboard TRUE or FALSE
    #' @param who    object id of sender
    #' @return list(plotit,copyit)
    set_flags = function(plotit=NULL,copyit=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_flags(plotit,copyit,"MC") }
      if(!is.null(plotit))
      {
        bool <- private$extract_boolean(plotit)
        if(!is.null(bool)) { private$flags$plotit <- bool  }
        else { message("plotit not set.") }
      }
      if(!is.null(copyit))
      {
        bool <- private$extract_boolean(copyit)
        if(!is.null(bool)) { private$flags$copyit <- bool  }
        else { message("copyit not set.") }
      }
      return(private$flags)
    },
    # public get methods ----
    #' @description
    #' Get all arguments and information
    #' @return list(oup_params,y_stoch_args,x_stoch_args,t_stoch_args,path_args,plot_info)
    get_all = function()
    {
      all <- list(oup_params = private$oup_params,
        y_stoch_args = private$y_stoch_args,
        x_stoch_args = private$x_stoch_args,
        t_stoch_args = private$t_stoch_args,
        path_args = private$path_args,
        plot_args = private$plot_args,
        plot_info = private$plot_info)
      return(all)
    },
    #' @description
    #' Get OUP parameters
    #' @return list(rho,mu,sigma)
    get_oup_params = function() { return(private$oup_params) },
    #' @description
    #' Get arguments for y as a stochastic state
    #' @return list(t,y,x,psi)
    get_y_stoch_args = function() { return(private$y_stoch_args) },
    #' @description
    #' Get arguments for x as a stochastic state
    #' @return list(s,y,r,phi)
    get_x_stoch_args = function() { return(private$x_stoch_args) },
    #' @description
    #' Get arguments for t as a stochastic state
    #' @return list(t,k,x,z,omega,Ppct)
    get_t_stoch_args = function() { return(private$t_stoch_args) },
    #' @description
    #' get path arguments
    #' @return list(paths,skip,seed,method)
    get_path_args = function() { return(private$path_args) },
    #' @description
    #' get plot arguments
    #' @return list(pmax,ptmax,first,last,zbeg,zend)
    get_plot_args = function() { return(private$plot_args) },
    #' @description
    #' Get information for plotting
    #' @return list(font,file,theme,3D,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
    #' @description
    #' Get current types for plot routines
    #' @return (list(types,description))
    get_plot_types = function()
    {
      text <- rbind(c("MonteCarlo groups, types and plots (default type is 0):"),c("  group  types  plots"),c("    1    -3,,0  ForwardPaths BackwardPaths BoundedPaths"),c("    2     -1,0  Mean Variance"),c("    3      0,1  Density Probability DoubleIntegral Option"),c("    4     -1,0  VisitingTimeModeMedianMean VisitingTimePercentiles FirstPassageTimeModeMedianMean FirstPassageTimePercentiles"),c("    5      0,1  VisitingTimeDensity VisitingTimeProbability FirstPassageTimeDensity FirstPassageTimeProbability"))
      return(list(types=private$plot_types,description=text))
    },
    #' @description
    #' Get flags for plotting and copying
    #' @return list(plotit,copyit)
    get_flags = function() { return(private$flags) },
    # public axis and sync methods ----
    #' @description
    #' Scale axes for y stochastic paths
    #' @return NULL
    axes_y_stoch = function()
    {
      # state
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      x <- private$y_stoch_args[[3]]
      m <- length(t)
      if(rho > 0)
      {
        G <- mu+(x-mu)*exp(-rho*(t[m]-t[1]))
        H <- (sigma^2/(2*rho)*(1-exp(-2*rho*(t[m]-t[1]))))^0.5
      }
      else
      {
        G <- x
        H <- sigma*(t[m]-t[1])^0.5
      }
      ydif <- G+H-mu
      if(ydif < abs(x-mu)) { ydif <- abs(x-mu) }
      ydif <- 3*ydif
      if(ydif < 1) { ydif <- 1}
      yscale <- 1
      while(ydif > yscale) { yscale <- 10*yscale }
      ydif <- round(ydif/yscale,1)*yscale
      yby <- ydif/100
      yfrom <- 0.5*(x+mu)-50*yby
      yto <- yfrom+100*yby
      yseq <- seq(from=yfrom,to=yto,by=yby)
      self$set_y_stoch_args(NULL,yseq,NULL,NULL)
      # density
      if(H > 0) { pmax <- 3/(2.506628274631*H) }
      else { pmax <- 1 }
      pscale <- 0.01
      while(pmax > pscale) { pscale <- 10*pscale }
      pmax <- round(pmax/pscale,2)*pscale
      self$set_plot_args(pmax,NULL,NULL,NULL)
      private$syncyxt <- 1
      private$forwardyt <- 1

      return(NULL)
    },
    #' @description
    #' Scale axes for x stochastic paths
    #' @return NULL
    axes_x_stoch = function()
    {
      # state
      backward <- self$BackwardPaths(who="MC")[[1]]
      minmax <- RcppOUPMCMinMax(backward)
      xmin <- minmax[[1]]
      xmax <- minmax[[2]]
      xdif <- xmax-xmin
      if(xdif < 1) { xdif <- 1}
      xscale <- 1
      while(xdif > xscale) { xscale <- 10*xscale }
      xdif <- round(xdif/xscale,1)*xscale
      xby <- xdif/100
      xfrom <- 0.5*(xmin+xmax)-50*xby
      xto <- xfrom+100*xby
      xseq <- seq(from=xfrom,to=xto,by=xby)
      self$set_x_stoch_args(NULL,xseq,NULL,NULL,NULL)

      private$syncyxt <- 2
      return(NULL)
    },
    #' @description
    #' Scale axes for t stochastic paths
    #' @return NULL
    axes_t_stoch = function()
    {
      # time
      pctv <- self$VisitingTimePercentiles(who="MC")[[1]]
      pctf <- self$FirstPassageTimePercentiles(who="MC")[[1]]
      t <- private$t_stoch_args[[1]]
      m <- length(t)
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      s <- t[1]
      if(is.na(pctv[3,1] && is.na(pctf[3,1]))) { tup <- 1.5*(t[m]-s) }
      else if(is.na(pctv[3,1])) { tup <- pctf[3,1]-s }
      else if(is.na(pctf[3,1])) { tup <- pctv[3,1]-s }
      else { tup <- max(pctv[3,1],pctf[3,1])-s }
      if(tup > 0.6*m*dt || tup < 0.4*m*dt)
      {
        tup <- 2*tup
        if(tup < 1) { tup <- 1 }
        tscale <- 1
        while(tup > tscale) { tscale <- 10*tscale }
        tup <- round(tup/tscale,2)*tscale
        tfrom <- s
        tto <- tup+s
        tby <- tup/100
        tseq <- seq(from=tfrom,to=tto,by=tby)
        self$set_t_stoch_args(t=tseq,NULL,NULL,NULL,NULL)
      }
      # density
      mmmv <- self$VisitingTimeModeMedianMean(who="MC")[[1]]
      mmmf <- self$FirstPassageTimeModeMedianMean(who="MC")[[1]]
      ptmax <- 1.2*max(mmmv[1,2],mmmf[1,2])
      pscale <- 0.01
      while(ptmax > pscale) { pscale <- 10*pscale }
      ptmax <- round(ptmax/pscale,2)*pscale
      self$set_plot_args(NULL,ptmax,NULL,NULL)
      private$syncyxt <- 3
      private$forwardyt <- 3

      return(NULL)
    },
    #' @description
    #' Synchronize arguments
    #' @return NULL
    sync_yxt_stoch = function()
    {
      syncyxt <- private$syncyxt
      # y stochastic
      if(syncyxt == 1)
      {
        y_stoch_args <- private$y_stoch_args
        t <- y_stoch_args[[1]]
        y <- y_stoch_args[[2]]
        x <- y_stoch_args[[3]]
        psi <- y_stoch_args[[4]]
        m <- length(t)
        if(m > 1) { ds <- (t[1]-t[m])/(m-1) }
        else { ds <- 0.05 }
        s <- seq(from=t[m],to=t[1],by=ds)
        self$set_x_stoch_args(s,x=y,y=x,NULL,phi=psi)
        self$set_t_stoch_args(t,NULL,x,NULL,NULL)
      }
      # x stochastic
      else if(syncyxt == 2)
      {
        x_stoch_args <- private$x_stoch_args
        s <- x_stoch_args[[1]]
        x <- x_stoch_args[[2]]
        y <- x_stoch_args[[3]]
        phi <- x_stoch_args[[5]]
        m <- length(s)
        if(m > 1) { dt <- (s[1]-s[m])/(m-1) }
        else { dt <- 0.05 }
        t <- seq(from=s[m],to=s[1],by=dt)
        self$set_y_stoch_args(t,y=x,x=y,psi=phi)
        self$set_t_stoch_args(t,k=y,NULL,NULL,NULL)
      }
      # t stochastic
      else
      {
        t_stoch_args <- private$t_stoch_args
        t <- t_stoch_args$t
        k <- t_stoch_args$k
        x <- t_stoch_args$x
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        s <- seq(from=t[m],to=t[1],by=-dt)
        self$set_y_stoch_args(t,NULL,x,NULL)
        self$set_x_stoch_args(s,y=k,NULL,NULL)
      }
      return(NULL)
    },
    # public undo methods ----
    #' @description
    #' Clear undo list and save current arguments to list
    #' @return 1
    undo_clear = function()
    {
      private$undo_args <- NULL
      private$undo_args <- list(list(oup_params=private$oup_params,
        y_stoch_args=private$y_stoch_args,
        x_stoch_args=private$x_stoch_args,
        t_stoch_args=private$t_stoch_args,
        path_args=private$path_args,
        plot_args=private$plot_args))
      private$undoIx <- 1

      return(1)
    },
    #' @description
    #' Save current arguments to undo list
    #' @return number of argument sets
    undo_save = function()
    {
      n <- length(private$undo_args)
      last_undo_args <- private$undo_args[[n]]
      not_equal <- TRUE
      if(private$lists_equal(last_undo_args[[1]],private$oup_params))
      {
        if(private$lists_equal(last_undo_args[[2]],private$y_stoch_args))
        {
          if(private$lists_equal(last_undo_args[[3]],private$x_stoch_args))
          {
            if(private$lists_equal(last_undo_args[[4]],private$t_stoch_args))
            {
              if(private$lists_equal(last_undo_args[[5]],private$path_args))
              {
                if(private$lists_equal(last_undo_args[[6]],private$plot_args))
                {
                  not_equal <- FALSE
                }
              }
            }
          }
        }
      }
      if(not_equal)
      {
        private$undo_args <- append(private$undo_args,
          list(list(oup_params=private$oup_params,
            y_stoch_args=private$y_stoch_args,
            x_stoch_args=private$x_stoch_args,
            t_stoch_args=private$t_stoch_args,
            path_args=private$path_args,
            plot_args=private$plot_args)))
        n <- n+1
        private$undoIx <- n
      }
      return(n)
    },
    #' @description
    #' Replace current arguments from undo list
    #' @param updn    positive to move up, negative to move down
    #' @return c(index of this argument set, number of argument sets)
    undo_undo = function(updn=NULL)
    {
      if(is.null(updn)) { updn <- -1 }
      n <- length(private$undo_args)
      if(updn > 0)
      {
        undoIx <- private$undoIx+1
        if(undoIx > n) { undoIx <- 1 }
      }
      else
      {
        undoIx <- undoIx-1
        if(undoIx < 1) { undoIx <- n }
      }
      these_undo <- private$undo_args[[undoIx]]
      oup <- these_undo[[1]]
      y_stoch <- these_undo[[2]]
      x_stoch <- these_undo[[3]]
      t_stoch <- these_undo[[4]]
      path <- these_undo[[5]]
      plot <- these_undo[[6]]
      self$set_oup_params(oup[[1]],oup[[2]],oup[[3]])
      self$set_y_stoch_args(y_stoch[[1]],y_stoch[[2]],y_stoch[[3]],y_stoch[[4]])
      self$set_x_stoch_args(x_stoch[[1]],x_stoch[[2]],x_stoch[[3]],x_stoch[[4]],x_stoch[[5]])
      self$set_t_stoch_args(t_stoch[[1]],t_stoch[[2]],t_stoch[[3]],t_stoch[[4]],t_stoch[[5]])
      self$set_path_args(path[[1]],path[[2]],path[[3]],path[[4]])
      self$set_plot_args(plot[[1]],plot[[2]],plot[[3]],plot[[4]],plot[[5]],plot[[6]])
      private$undoIx <- undoIx

      return(c(undoIx,n))
    },
    # public calculate methods ----
    #' @description
    #' Calculate and plot forward paths
    #' @param t       vector of m times -inf<t<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(forward(m,paths))
    ForwardPaths = function(t=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      forwardyt <- private$forwardyt
      if(forwardyt == 1) { self$set_y_stoch_args(t,NULL,x,NULL) }
      else { self$set_t_stoch_args(t,NULL,x,NULL,NULL) }
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      if(forwardyt == 1)
      {
        t <- private$y_stoch_args[[1]]
        x <- private$y_stoch_args[[3]]
      }
      else
      {
        t <- private$t_stoch_args[[1]]
        k <- private$t_stoch_args[[2]]
        x <- private$t_stoch_args[[3]]
      }
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      p1 <- private$plot_args[[3]]
      pn <- private$plot_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(forwardyt == 1)
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          m <- length(t)
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }

          private$yforward <- forward
        }
      }
      else
      {
        forward <- private$tforward
        if(is.null(forward))
        {
          m <- length(t)
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }

          private$tforward <- forward
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotForwardPaths()) }
        else if(copyit == TRUE)
        {
          if(forwardyt == 1) { clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Forward Paths",rep("",pn-p1+1)),c("s",t[1],rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,forward[,p1:pn,drop=FALSE])) }
          else { clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Forward Paths",rep("",pn-p1+1)),c("k",k,rep("",pn-p1)),c("s",t[1],rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,forward[,p1:pn,drop=FALSE])) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(forward=forward))
    },
    #' @description
    #' Calculate and plot backward paths
    #' @param s       vector of m times -inf<s<inf
    #' @param y       terminal state -inf<y<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(backward(m,paths))
    BackwardPaths = function(s=NULL,y=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,NULL,y,NULL,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      y <- private$x_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      p1 <- private$plot_args[[3]]
      pn <- private$plot_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      backward <- private$xbackward
      if(is.null(backward))
      {
        m <- length(s)
        if(m > 1) { ds <- (s[1]-s[m])/(m-1) }
        else { ds <- 0.05 }
        stdnorm <- private$xstdnorm
        if(is.null(stdnorm))
        {
          stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
          private$xstdnorm <- stdnorm
        }
        if(method == 4) { backward <- RcppOUPMCBackwardPathRungeKutta(stdnorm,y,m,skip,ds,rho,mu,sigma) }
        else { backward <- RcppOUPMCBackwardPathIntegralEquation(stdnorm,y,m,skip,ds,rho,mu,sigma) }
        private$xbackward <- backward
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotBackwardPaths()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Backward Paths",rep("",pn-p1+1)),c("t",s[1],rep("",pn-p1)),c("y",y,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("s",paste0("path",p1:pn)),cbind(s,backward[,p1:pn,drop=FALSE]))
          private$CopyToClipboard(clip)
        }
      }
      return(list(backward=backward))
    },
    #' @description
    #' Calculate and plot bounded paths
    #' @param t       vector of m times -inf<t<inf
    #' @param x       initial state -inf<x<inf
    #' @param k       threshold or bound -inf<k<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(bounded(m,paths))
    BoundedPaths = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,NULL,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      p1 <- private$plot_args[[3]]
      pn <- private$plot_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      bounded <- private$tbounded
      if(is.null(bounded))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        stdnorm <- private$tstdnorm
        if(is.null(stdnorm))
        {
          stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
          private$tstdnorm <- stdnorm
        }
        if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
        else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
        bounded <- bndfpt[1:m,,drop=FALSE]
        fpt <- bndfpt[m+1,,drop=FALSE]
        private$tbounded <- bounded
        private$tfpt <- fpt
     }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotBoundedPaths()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Bounded Paths",rep("",pn-p1+1)),c("k",k,rep("",pn-p1)),c("s",t[1],rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,bounded[,p1:pn,drop=FALSE]))
          private$CopyToClipboard(clip)
        }
      }
      return(list(bounded=bounded))
    },
    #' @description
    #' Calculate and plot means
    #' @param t       vector of m times -inf<t<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(G(m))
    Mean = function(t=NULL,x=NULL,rho=NULL,mu=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_y_stoch_args(t,NULL,x,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      means <- private$G
      if(is.null(means))
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          m <- length(t)
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$yforward <- forward
        }
        n <- length(y)
        mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
        means <- mvdpd[,1,drop=FALSE]
        variances <- mvdpd[,2,drop=FALSE]
        densities <- mvdpd[,3:(n+2),drop=FALSE]
        probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
        doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
        private$G <- means
        private$H2 <- variances
        private$p <- densities
        if(psi <= 0)
        {
          private$Pneg <- probabilities
          private$PPneg <- doubleintegrals
        }
        else
        {
          private$Ppos <- probabilities
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotMean()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("Mean",""),c("s",t[1]),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","G"),cbind(t,means))
          private$CopyToClipboard(clip)
        }
      }
      return(list(G=means))
    },
    #' @description
    #' Calculate and plot variances
    #' @param t       vector of m times -inf<t<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(H2(m))
    Variance = function(t=NULL,x=NULL,rho=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,NULL,sigma)
      self$set_y_stoch_args(t,NULL,x,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      variances <- private$H2
      if(is.null(variances))
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          m <- length(t)
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$yforward <- forward
        }
        n <- length(y)
        mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
        means <- mvdpd[,1,drop=FALSE]
        variances <- mvdpd[,2,drop=FALSE]
        densities <- mvdpd[,3:(n+2),drop=FALSE]
        probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
        doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
        private$G <- means
        private$H2 <- variances
        private$p <- densities
        if(psi <= 0)
        {
          private$Pneg <- probabilities
          private$PPneg <- doubleintegrals
        }
        else
        {
          private$Ppos <- probabilities
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVariance()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("Variance",""),c("s",t[1]),c("rho",rho),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","H\u00B2"),cbind(t,variances))
          private$CopyToClipboard(clip)
        }
      }
      return(list(H2=variances))
    },
    #' @description
    #' Calculate and plot densities
    #' @param t       vector of m times -inf<t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(p(m,n))
    Density = function(t=NULL,y=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,x,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      m <- length(t)
      n <- length(y)
      densities <- private$p
      if(is.null(densities))
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$yforward <- forward
        }
        mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
        means <- mvdpd[,1,drop=FALSE]
        variances <- mvdpd[,2,drop=FALSE]
        densities <- mvdpd[,3:(n+2),drop=FALSE]
        probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
        doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
        private$G <- means
        private$H2 <- variances
        private$p <- densities
        if(psi <= 0)
        {
          private$Pneg <- probabilities
          private$PPneg <- doubleintegrals
        }
        else
        {
          private$Ppos <- probabilities
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDensity()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",n)),c("Transition Densities",rep("",n)),c("s",t[1],rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("p(t,y)",y),cbind(t,densities))
          private$CopyToClipboard(clip)
        }
      }
      return(list(p=densities))
    },
    #' @description
    #' Calculate and plot probabilities
    #' @param t       vector of m times -inf<t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param x       initial state -inf<x<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(P(m,n))
    Probability = function(t=NULL,y=NULL,x=NULL,psi=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,x,psi)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      m <- length(t)
      n <- length(y)
      # calculate ----
      if(psi <= 0) { probabilities <- private$Pneg }
      else { probabilities <- private$Ppos }
      if(is.null(probabilities))
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$yforward <- forward
        }
        mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
        means <- mvdpd[,1,drop=FALSE]
        variances <- mvdpd[,2,drop=FALSE]
        densities <- mvdpd[,3:(n+2),drop=FALSE]
        probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
        doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
        private$G <- means
        private$H2 <- variances
        private$p <- densities
        if(psi <= 0)
        {
          private$Pneg <- probabilities
          private$PPneg <- doubleintegrals
        }
        else
        {
          private$Ppos <- probabilities
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotProbability()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",n)),c("Transition Probabilities",rep("",n)),c("s",t[1],rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("P(t,y)",y),cbind(t,probabilities))
          private$CopyToClipboard(clip)
        }
      }
      return(list(P=probabilities))
    },
    #' @description
    #' Calculate and plot double integrals
    #' @param t       vector of m times -inf<t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param x       initial state -inf<x<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(PP(m,n))
    DoubleIntegral = function(t=NULL,y=NULL,x=NULL,psi=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,x,psi)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      m <- length(t)
      n <- length(y)
      if(psi <= 0) { doubleintegrals <- private$PPneg  }
      else { doubleintegrals <- private$PPpos }
      if(is.null(doubleintegrals))
      {
        forward <- private$yforward
        if(is.null(forward))
        {
          if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
          else { dt <- 0.05 }
          stdnorm <- private$ystdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$ystdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$yforward <- forward
        }
        mvdpd <- RcppOUPMCForwardCountY(forward,y,psi)
        means <- mvdpd[,1,drop=FALSE]
        variances <- mvdpd[,2,drop=FALSE]
        densities <- mvdpd[,3:(n+2),drop=FALSE]
        probabilities <- mvdpd[,(n+3):(2*n+2),drop=FALSE]
        doubleintegrals <- mvdpd[,(2*n+3):(3*n+2),drop=FALSE]
        private$G <- means
        private$H2 <- variances
        private$p <- densities
        if(psi <= 0)
        {
          private$Pneg <- probabilities
          private$PPneg <- doubleintegrals
        }
        else
        {
          private$Ppos <- probabilities
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDoubleIntegral()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",n)),c("Double Integrals",rep("",n)),c("s",t[1],rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("\u2119(t,y)",y),cbind(t,doubleintegrals))
          private$CopyToClipboard(clip)
        }
      }
      return(list(PP=doubleintegrals))
    },
    #' @description
    #' Calculate and plot option prices
    #' @param s       vector of m times -inf<s<t
    #' @param x       vector of n states -inf<x<inf
    #' @param y       terminal state -inf<y<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(OO(mxn))
    Option = function(s=NULL,x=NULL,y=NULL,r=NULL,phi=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,x,y,r,phi)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      y <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      m <- length(s)
      n <- length(x)
      if(phi <= 0) { options <- private$OOneg }
      else { options <- private$OOpos }
      if(is.null(options))
      {
        if(m > 1) { ds <- (s[1]-s[m])/(m-1) }
        else { ds <- 0.05 }
        backward <- private$xbackward
        if(is.null(backward))
        {
          stdnorm <- private$xstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$xstdnorm <- stdnorm
          }
          if(method == 4) { backward <- RcppOUPMCBackwardPathRungeKutta(stdnorm,y,m,skip,ds,rho,mu,sigma) }
          else { backward <- RcppOUPMCBackwardPathIntegralEquation(stdnorm,y,m,skip,ds,rho,mu,sigma) }
          private$xbackward <- backward
        }
        dpo <- RcppOUPMCBackwardCountX(backward,x,phi,rho,r,ds)
        densities <- dpo[,1:n,drop=FALSE]
        probabilities <- dpo[,(n+1):(2*n),drop=FALSE]
        options <- dpo[,(2*n+1):(3*n),drop=FALSE]
        private$o <- densities
        if(phi <= 0)
        {
          private$Oneg <- probabilities
          private$OOneg <- options
        }
        else
        {
          private$Opos <- probabilities
          private$OOpos <- options
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotOption()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",n)),c("Options",rep("",n)),c("t",s[1],rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("phi",phi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options))
          private$CopyToClipboard(clip)
        }
      }
      return(list(OO=options))
    },
    #' @description
    #' Calculate and plot visiting time mode, median and mean
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(vtmmm(3x3))
    VisitingTimeModeMedianMean = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,0,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      modemedianmean <- private$vtmmm
      if(is.null(modemedianmean))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        forward <- private$tforward
        if(is.null(forward))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$tforward <- forward
        }
        pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$vtmmm <- modemedianmean
        private$vtpct <- percentile
        private$pv <- densities
        private$Pv <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVisitingTimeModeMedianMean()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",2)),c("Visiting Time Mode, Median, Mean",rep("",2)),c("k",k,""),c("s",t[1],""),c("x",x,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("tv","pv","Pv"),modemedianmean)
          private$CopyToClipboard(clip)
        }
      }
      return(list(vtmmm=modemedianmean))
    },
    #' @description
    #' Calculate and plot visiting time percentiles
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param Ppct    probability for a percentile 0.01<=Ppct<=0.99
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(vtpct(3x3))
    VisitingTimePercentiles = function(t=NULL,k=NULL,x=NULL,Ppct=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,0,Ppct)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      percentile <- private$vtpct
      if(is.null(percentile))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        forward <- private$tforward
        if(is.null(forward))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$tforward <- forward
        }
        pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$vtmmm <- modemedianmean
        private$vtpct <- percentile
        private$pv <- densities
        private$Pv <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVisitingTimePercentiles()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",2)),c("Visiting Time Percentiles",rep("",2)),c("k",k,""),c("s",t[1],""),c("x",x,""),c("P%",Ppct,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("t%","pv","Pv"),percentile)
          private$CopyToClipboard(clip)
        }
      }
      return(list(vtpct=percentile))
    },
    #' @description
    #' Calculate and plot visiting time densities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(pv(m))
    VisitingTimeDensity = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,0,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      densities <- private$pv
      if(is.null(densities))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        forward <- private$tforward
        if(is.null(forward))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$tforward <- forward
        }
        pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$vtmmm <- modemedianmean
        private$vtpct <- percentile
        private$pv <- densities
        private$Pv <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVisitingTimeDensity()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("Visiting Time Density",""),c("k",k),c("s",t[1]),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","pv"),cbind(t,pv))
          private$CopyToClipboard(clip)
        }
      }
      return(list(pv=densities))
    },
    #' @description
    #' Calculate and plot visiting time probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(Pv(m))
    VisitingTimeProbability = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,0,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      probabilities <- private$Pv
      if(is.null(probabilities))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        forward <- private$tforward
        if(is.null(forward))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$tforward <- forward
        }
        pctdp <- RcppOUPMCForwardCountT(forward,k,dt,rho,mu,sigma,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$vtmmm <- modemedianmean
        private$vtpct <- percentile
        private$pv <- densities
        private$Pv <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVisitingTimeProbability()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("Visiting Time Probability",""),c("k",k),c("s",t[1]),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","Pv"),cbind(t,Pv))
          private$CopyToClipboard(clip)
        }
      }
      return(list(Pv=probabilities))
    },
    #' @description
    #' Calculate and plot first passage time mode, median and mean
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(fptmmm(3x3))
    FirstPassageTimeModeMedianMean = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,1,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      modemedianmean <- private$fptmmm
      if(is.null(modemedianmean))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        fpt <- private$tfpt
        if(is.null(fpt))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          bounded <- bndfpt[1:m,,drop=FALSE]
          fpt <- bndfpt[m+1,,drop=FALSE]
          private$tbounded <- bounded
          private$tfpt <- fpt
        }
        pctdp <- RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$fptmmm <- modemedianmean
        private$fptpct <- percentile
        private$pf <- densities
        private$Pf <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotFirstPassageTimeModeMedianMean()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",2)),c("First Passage Time Mode, Median and Mean",rep("",2)),c("k",k,""),c("s",t[1],""),c("x",x,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("tf","pf","Pf"),modemedianmean)
          private$CopyToClipboard(clip)
        }
      }
      return(list(fptmmm=modemedianmean))
    },
    #' @description
    #' Calculate and plot first passage time percentiles
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param Ppct    probability for a percentile 0.01<=Ppct<=0.99
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(fptpct(3x3))
    FirstPassageTimePercentiles = function(t=NULL,k=NULL,x=NULL,Ppct=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,1,Ppct)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      percentile <- private$fptpct
      if(is.null(percentile))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        fpt <- private$tfpt
        if(is.null(fpt))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          bounded <- bndfpt[1:m,,drop=FALSE]
          fpt <- bndfpt[m+1,,drop=FALSE]
          private$tbounded <- bounded
          private$tfpt <- fpt
        }
        pctdp <- RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$fptmmm <- modemedianmean
        private$fptpct <- percentile
        private$pf <- densities
        private$Pf <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotFirstPassageTimePercentiles()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",rep("",2)),c("First Passage Time Percentiles",rep("",2)),c("k",k,""),c("s",t[1],""),c("x",x,""),c("P%",Ppct,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("t%","pf","Pf"),percentile)
          private$CopyToClipboard(clip)
        }
      }
      return(list(fptpct=percentile))
    },
    #' @description
    #' Calculate and plot first passage time densities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(pf(m))
    FirstPassageTimeDensity = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,1,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      densities <- private$pf
      if(is.null(densities))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        fpt <- private$tfpt
        if(is.null(fpt))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          bounded <- bndfpt[1:m,,drop=FALSE]
          fpt <- bndfpt[m+1,,drop=FALSE]
          private$tbounded <- bounded
          private$tfpt <- fpt
        }
        pctdp <- RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$fptmmm <- modemedianmean
        private$fptpct <- percentile
        private$pf <- densities
        private$Pf <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotFirstPassageTimeDensity()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("First Passage Time Density",""),c("k",k),c("s",t[1]),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","pf"),cbind(t,pf))
          private$CopyToClipboard(clip)
        }
      }
      return(list(pf=densities))
    },
    #' @description
    #' Calculate and plot first passage time probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param paths   number of paths 1<paths<1,000,000
    #' @param skip    subdivide time interval but report at times t 1<=skip<=50
    #' @param seed    seed for random number generators -inf<seed<inf
    #' @param method  4 for 4th order Runge-Kutta, otherwise integral equation
    #' @param who     object id of caller
    #' @return list(Pf(m))
    FirstPassageTimeProbability = function(t=NULL,k=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,paths=NULL,skip=NULL,seed=NULL,method=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,x,1,NULL)
      self$set_path_args(paths,skip,seed,method)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      probabilities <- private$Pf
      if(is.null(probabilities))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        fpt <- private$tfpt
        if(is.null(fpt))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          bounded <- bndfpt[1:m,,drop=FALSE]
          fpt <- bndfpt[m+1,,drop=FALSE]
          private$tbounded <- bounded
          private$tfpt <- fpt
        }
        pctdp <- RcppOUPMCBoundedCountT(fpt,m,dt,Ppct)
        modemedianmean <- pctdp[1:3,1:3,drop=FALSE]
        percentile <- rbind(pctdp[4,1:3,drop=FALSE],pctdp[2,1:3,drop=FALSE],pctdp[5,1:3,drop=FALSE])
        densities <- pctdp[,4,drop=FALSE]
        probabilities <- pctdp[,5,drop=FALSE]
        private$fptmmm <- modemedianmean
        private$fptpct <- percentile
        private$pf <- densities
        private$Pf <- probabilities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotFirstPassageTimeProbability()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Monte Carlo",""),c("First Passage Time Probability",""),c("k",k),c("s",t[1]),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","Pf"),cbind(t,Pf))
          private$CopyToClipboard(clip)
        }
      }
      return(list(Pf=probabilities))
    },
    # public plot methods ----
    #' @description
    #' Plot forward paths
    #' @param type = -3,...,0, or 'n','p','d' for next, previous, default
    #' @param first  first path to plot
    #' @param last   last path to plot
    #' @param title  text for plot title
    #' @param xaxis  text for x-axis label
    #' @param yaxis  text for y-axis label
    #' @param tbeg   begin value for time axis
    #' @param tend   end value for time axis
    #' @return plot
    PlotForwardPaths = function(type=NULL,first=NULL,last=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,1)[[1]]
      self$set_plot_args(NULL,NULL,first,last,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      if(type < -1.5)
      {
        private$forwardyt <- 1
        t <- private$y_stoch_args[[1]]
        x <- private$y_stoch_args[[3]]
        forward <- private$yforward #protect against recursive call
      }
      else
      {
        private$forwardyt <- 3
        t <- private$t_stoch_args[[1]]
        k <- private$t_stoch_args[[2]]
        x <- private$t_stoch_args[[3]]
        forward <- private$tforward #protect against recursive call
      }
      p1 <- private$plot_args$first
      pn <- private$plot_args$last
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      if(is.null(forward)) { forward <- self$ForwardPaths(who="MC")[[1]] }
      m <- length(t)
      n <- ncol(forward)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        forward <- forward[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      if(p1 > 1 || pn < n)
      {
        forward <- forward[,p1:pn,drop=FALSE]
        n <- ncol(forward)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(type < -1.5) { clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Forward Paths",rep("",pn-p1+1)),c("s",s,rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,forward)) }
        else { clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Forward Paths",rep("",pn-p1+1)),c("k",k,rep("",pn-p1)),c("s",s,rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,forward)) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        if(is.null(title)) { title <- paste(sep="","Forward Paths[",p1,":",pn,"]") }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      if(n < 6) { pathline <- list(color=red$d,width=4) }
      else if(n < 16) { pathline <- list(color=red$d,width=3) }
      else if(n < 31) { pathline <- list(color=red$d,width=2) }
      else { pathline <- list(color=red$d,width=1) }
      threshline <- list(color=gry$c,width=4)
      if(Ixend-Ixbeg < 6) { pathmarker <- list(color=red$b,size=6,symbol="square") }
      else if(Ixend-Ixbeg < 16) { pathmarker <- list(color=red$b,size=5,symbol="square") }
      else if(Ixend-Ixbeg < 31) { pathmarker <- list(color=red$b,size=4,symbol="square") }
      else { pathmarker <- list(color=red$b,size=3,symbol="square") }
      #OUP_MC_ForwardPaths2Dtmark and OUP_MC_ForwardPaths2Dtline
      if(type < -1.5)
      {
        if(labels == TRUE)
        {
          syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>y</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>y</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly()
        #OUP_MC_ForwardPaths2Dtmark
        if(type < -2.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_ForwardPaths2Dtmark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=forward[,j],y=t,mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_ForwardPaths2Dtline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_ForwardPaths2Dtline")
          fig <- plot_ly()
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=forward[,j],y=t,mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      #OUP_MC_ForwardPaths2Dymark and OUP_MC_ForwardPaths2Dyline
      else
      {
        if(labels == TRUE)
        {
          syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>y</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly() %>%
          add_trace(.,type="scattergl",x=c(t[1],t[m]),y=c(k,k),name="<i>k</i>",mode="lines",line=threshline,hoverinfo="x+y")
        #OUP_MC_ForwardPaths2Dymark
        if(type < -0.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_ForwardPaths2Dymark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=t,y=forward[,j],mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_ForwardPaths2Dyline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_ForwardPaths2Dyline")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=t,y=forward[,j],mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      fig <- fig %>% toWebGL()
      return(fig)
    },
     #' @description
    #' Plot backward paths
    #' @param type = -3,...,0, or 'n','p','d' for next, previous, default
    #' @param first  first path to plot
    #' @param last   last path to plot
    #' @param title  text for plot title
    #' @param xaxis  text for x-axis label
    #' @param yaxis  text for y-axis label
    #' @param sbeg   begin value for time axis
    #' @param send   end value for time axis
    #' @param copyit TRUE or FALSE
    #' @return plot
    PlotBackwardPaths = function(type=NULL,first=NULL,last=NULL,title=NULL,xaxis=NULL,yaxis=NULL,sbeg=NULL,send=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,1)[[1]]
      self$set_plot_args(NULL,NULL,first,last,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      y <- private$x_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      p1 <- private$plot_args$first
      pn <- private$plot_args$last
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      backward <- private$backward #protect against recursive call
      if(is.null(backward)) { backward <- self$BackwardPaths(who="MC")[[1]] }
      m <- length(s)
      n <- ncol(backward)
      t <- s[1]
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[2]]
      Ixend <- Inx[[1]]
      if(Ixbeg > 1 || Ixend < m)
      {
        s <- s[Ixbeg:Ixend]
        backward <- backward[Ixbeg:Ixend,,drop=FALSE]
        m <- length(s)
      }
      if(p1 > 1 || pn < n)
      {
        backward <- backward[,p1:pn,drop=FALSE]
        n <- ncol(backward)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Backward Paths",rep("",pn-p1+1)),c("t",t,rep("",pn-p1)),c("y",y,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("s",paste0("path",p1:pn)),cbind(s,backward))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- paste(sep="","Backward Paths[",p1,":",pn,"]") }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      if(n < 6) { pathline <- list(color=red$d,width=4) }
      else if(n < 16) { pathline <- list(color=red$d,width=3) }
      else if(n < 31) { pathline <- list(color=red$d,width=2) }
      else { pathline <- list(color=red$d,width=1) }
      if(Ixend-Ixbeg < 6) { pathmarker <- list(color=red$b,size=4,symbol="square") }
      else if(Ixend-Ixbeg < 16) { pathmarker <- list(color=red$b,size=3,symbol="square") }
      else if(Ixend-Ixbeg < 31) { pathmarker <- list(color=red$b,size=2,symbol="square") }
      else { pathmarker <- list(color=red$b,size=1,symbol="square") }
      #OUP_MC_BackwardPaths2Dsmark and OUP_MC_BackwardPaths2Dsline
      if(type < -1.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly()
        #OUP_MC_BackwardPaths2Dsmark
        if(type < -2.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BackwardPaths2Dsmark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=backward[,j],y=s,mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_BackwardPaths2Dsline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BackwardPaths2Dsline")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=backward[,j],y=s,mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      #OUP_MC_BackwardPaths2Dxmark and OUP_MC_BackwardPaths2Dxline
      else
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>s</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>s</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>x</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly()
        #OUP_MC_BackwardPaths2Dxmark
        if(type < -0.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BackwardPaths2Dxmark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=s,y=backward[,j],mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_BackwardPaths2Dxline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BackwardPaths2Dxline")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=s,y=backward[,j],mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' Plot bounded paths
    #' @param type = -3,...,0, or 'n','p','d' for next, previous, default
    #' @param first  first path to plot
    #' @param last   last path to plot
    #' @param title  text for plot title
    #' @param xaxis  text for x-axis label
    #' @param yaxis  text for y-axis label
    #' @param tbeg   begin value for time axis
    #' @param tend   end value for time axis
    #' @return plot
    PlotBoundedPaths = function(type=NULL,first=NULL,last=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,1)[[1]]
      self$set_plot_args(NULL,NULL,first,last,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      p1 <- private$plot_args$first
      pn <- private$plot_args$last
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      bounded <- private$bounded #protect against recursive call
      if(is.null(bounded)) { bounded <- self$BoundedPaths(who="MC")[[1]] }
      m <- length(t)
      n <- ncol(bounded)
      s <- t[1]
      paths <- n
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        bounded <- bounded[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      if(p1 > 1 || pn < n)
      {
        bounded <- bounded[,p1:pn,drop=FALSE]
        n <- ncol(bounded)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",pn-p1+1)),c("Bounded Paths",rep("",pn-p1+1)),c("k",k,rep("",pn-p1)),c("s",s,rep("",pn-p1)),c("x",x,rep("",pn-p1)),c("rho",rho,rep("",pn-p1)),c("mu",mu,rep("",pn-p1)),c("sigma",sigma,rep("",pn-p1)),c("paths",paths,rep("",pn-p1)),c("skip",skip,rep("",pn-p1)),c("seed",seed,rep("",pn-p1)),c("method",method,rep("",pn-p1)),c("t",paste0("path",p1:pn)),cbind(t,bounded))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- paste(sep="","Bounded Paths[",p1,":",pn,"]") }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      if(n < 6) { pathline <- list(color=red$d,width=4) }
      else if(n < 16) { pathline <- list(color=red$d,width=3) }
      else if(n < 31) { pathline <- list(color=red$d,width=2) }
      else { pathline <- list(color=red$d,width=1) }
      threshline <- list(color=gry$c,width=4)
      if(Ixend-Ixbeg < 6) { pathmarker <- list(color=red$b,size=6,symbol="square") }
      else if(Ixend-Ixbeg < 16) { pathmarker <- list(color=red$b,size=5,symbol="square") }
      else if(Ixend-Ixbeg < 31) { pathmarker <- list(color=red$b,size=4,symbol="square") }
      else { pathmarker <- list(color=red$b,size=3,symbol="square") }
      #OUP_MC_BoundedPaths2Dtmark and OUP_MC_BoundedPaths2Dtline
      if(type < -1.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>y</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>y</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly() %>%
          add_trace(.,type="scattergl",x=c(k,k),y=c(t[1],t[m]),name="<i>k</i>",mode="lines",line=threshline,hoverinfo="x+y")
        #OUP_MC_BoundedPaths2Dtmark
        if(type < -2.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BoundedPaths2Dtmark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=bounded[,j],y=t,mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_BoundedPaths2Dtline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BoundedPaths2Dtline")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=bounded[,j],y=t,mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      #OUP_MC_BoundedPaths2Dymark and OUP_MC_BoundedPaths2Dyline
      else
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>y</i>" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        fig <- plot_ly() %>%
          add_trace(.,type="scattergl",x=c(t[1],t[m]),y=c(k,k),name="<i>k</i>",mode="lines",line=threshline,hoverinfo="x+y")
        #OUP_MC_BoundedPaths2Dymark
        if(type < -0.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BoundedPaths2Dymark")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=t,y=bounded[,j],mode="markers",marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        #OUP_MC_BoundedPaths2Dyline
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_BoundedPaths2Dyline")
          j <- 0
          while(j < n)
          {
            j <- j+1
            fig <- add_trace(fig,type="scattergl",x=t,y=bounded[,j],mode="lines+markers",line=pathline,marker=pathmarker,hoverinfo="text",hovertext=paste0("path[",j+p1-1,"]"))
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot means
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotMean = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,2)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      cyn <- private$plot_colors$cyn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      means <- private$G #protect against recursive call
      if(is.null(means)) { means <- self$Mean(who="MC")[[1]] }
      densities <- private$p #no plot or copy
      if(is.null(densities)) { densities <- self$Density(who="MC")[[1]] }
      m <- length(t)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("Mean",""),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","G"),cbind(t,means))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Mean" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      if(labels == TRUE)
      {
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      lookdown <- list(text=xaxis)
      if(is.null(yaxis)) { yaxis <- "<i>G</i>(<i>t</i>|<i>s,x</i>)" }
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      meanline <- list(color=cyn$d,width=4)
      horzline <- list(color=gry$d,width=1)
      #OUP_MC_Means2Dpaths
      if(type < -0.5)
      {
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Means2Dpaths")
        legendpos <- list(orientation="h",x=0.5,y=1.0,xanchor="center")
        fig <- plot_ly() %>%
          add_trace(.,type="heatmap",x=t,y=y,z=t(densities),name="paths",showscale=FALSE,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="x+y",legendgroup="G") %>%
          add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),mode="lines",line=horzline,hoverinfo="x+y",legendgroup="G",showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      #OUP_MC_Means2D
      else
      {
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Means2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),mode="lines",line=horzline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",hoverinfo="x+y",line=meanline) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot variances
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotVariance = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,2)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      means <- private$G #no plot or copy
      if(is.null(means)) { means <- self$Mean(who="MC")[[1]] }
      variances <- private$H2 #protect against recursive call
      if(is.null(variances)) { variances <- self$Variance(who="MC")[[1]] }
      densities <- private$p #no plot or copy
      if(is.null(densities)) { densities <- self$Density(who="MC")[[1]] }
      sqrts <- variances^0.5
      meansplus <- means+sqrts
      meansminus <- means-sqrts
      varsigma2 <- sigma^2/(2*rho)
      m <- length(t)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        variances <- variances[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        meansplus <- meansplus[Ixbeg:Ixend]
        meansminus <- meansminus[Ixbeg:Ixend]
        m <- length(t)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("Variance",""),c("s",s),c("rho",rho),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","H\u00B2"),cbind(t,variances))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        if(is.null(title)) { title <- "Variance" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      #OUP_MC_Variances2Dpaths
      if(type < -0.5)
      {
        if(labels == TRUE)
        {
          syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>G</i>(<i>t</i>|<i>s,x</i>)&plusmn;<i>H</i>(<i>t</i>|<i>s</i>)" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        meanline <- list(color=cyn$d,width=4)
        varianceline <- list(color=mgn$d,width=4)
        horzline <- list(color=gry$d,width=1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Variances2Dpaths")
        legendpos <- list(orientation="h",x=0.5,y=1.0,xanchor="center")
        fig <- plot_ly() %>%
          add_trace(.,type="heatmap",x=t,y=y,z=t(densities),name="paths",showscale=FALSE,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter",x=t,y=meansplus,name="<i>G</i>(<i>t</i>)&plusmn;<i>H</i>(<i>t</i>)",mode="lines",line=varianceline,hoverinfo="x+y",legendgroup="G+H") %>%
          add_trace(.,type="scatter",x=t,y=meansminus,name="<i>G</i>(<i>t</i>)&plusmn;<i>H</i>(<i>t</i>)",mode="lines",line=varianceline,hoverinfo="x+y",legendgroup="G+H",showlegend=FALSE) %>%
          add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="x+y",legendgroup="G+H",showlegend=FALSE) %>%
          add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),mode="lines",line=horzline,hoverinfo="x+y",legendgroup="G+H",showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      #OUP_MC_Variances2D
      else
      {
        if(labels == TRUE)
        {
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        if(is.null(yaxis)) { yaxis <- "<i>H</i><sup>2</sup>(<i>t</i>|<i>s</i>)" }
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        varianceline <- list(color=mgn$d,width=4)
        horzline <- list(color=gry$d,width=1)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Variance2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=t,y=variances,name="<i>H</i><sup>2</sup>(<i>t</i>)",mode="lines",line=varianceline) %>%
          add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(varsigma2,varsigma2),mode="lines",line=horzline,hoverinfo="x+y",showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot transition densities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotDensity = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      densities <- private$p #protect against recursive call
      if(is.null(densities)) { densities <- self$Density(who="MC")[[1]] }
      m <- length(t)
      n <- length(y)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",n)),c("Transition Densities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("p(t,y)",y),cbind(t,densities))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Transition Densities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_Density2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>y</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>y</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(is.nan(pmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,pmax)) }
        densityline <- list(color=blu$d,width=4)
        densitymarker <- list(color=blu$d,width=2,symbol="square")
        markeropacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Density2D")
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=y,y=densities[i,],mode="lines",line=densityline,opacity=markeropacity,hoverinfo="x+y")
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter",x=y,y=densities[i,],mode="markers",marker=densitymarker,opacity=markeropacity,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_Density3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>p</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        if(x < mu) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(-0.03*pmax,1.03*pmax),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        densityline <- list(color=blu$e,width=8)
        densitymarker <- list(color=blu$e,size=4,symbol="square")
        gradient <- list(c(0,blu$c),c(1,blu$c))
        markeropacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_Density3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=y,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=densities,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE)
        skip<- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="pt")
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=densities[i,],mode="markers",marker=densitymarker,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="pt",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=densities,name="<i>p</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot transition probabilities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotProbability = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(psi <= 0) { probabilities <- private$Pneg } #protect against recursive call
      else { probabilities <- private$Ppos }
      if(is.null(probabilities)) { probabilities <- self$Probability(who="MC")[[1]] }
      densities <- private$p
      m <- length(t)
      n <- length(y)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        probabilities <- probabilities[Ixbeg:Ixend,,drop=FALSE]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        probabilities <- probabilities[,Ixbeg:Ixend,drop=FALSE]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",n)),c("Transition Probabilities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("P(t,y)",y),cbind(t,probabilities))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>y</i>=",esym,format(psi,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Transition Probabilities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_Probability2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>y</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>y</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(psi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,1)) }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,1),side="right") }
        probabilityline <- list(color=grn$d,width=4)
        probabilitymarker <- list(color=grn$d,width=2,symbol="square")
        markeropacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Probability2D")
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=y,y=probabilities[i,],mode="lines",line=probabilityline,opacity=markeropacity,hoverinfo="x+y")
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter",x=y,y=probabilities[i,],mode="markers",marker=probabilitymarker,opacity=markeropacity,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_Probability3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(probabilities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        if(psi > 0) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(psi == 0) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview, zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        probabilityline <- list(color=grn$e,width=8)
        probabilitymarker <- list(color=grn$e,size=4,symbol="square")
        gradient <- list(c(0,grn$c),c(1,grn$c))
        markeropacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_Probability3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=y,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=densities,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE)
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="Pt",showlegend=TRUE)
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=probabilities[i,],mode="markers",marker=probabilitymarker,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="Pt",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=probabilities,name="<i>P</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot double integrals of transition densities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotDoubleIntegral = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      x <- private$y_stoch_args[[3]]
      psi <- private$y_stoch_args[[4]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(psi <= 0) { doubleintegrals <- private$PPneg  } #protect against recursive call
      else { doubleintegrals <- private$PPpos }
      if(is.null(doubleintegrals)) { doubleintegrals <- self$DoubleIntegral(who="MC")[[1]] }
      densities <- private$p
      m <- length(t)
      n <- length(y)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[Ixbeg:Ixend,,drop=FALSE]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[,Ixbeg:Ixend,drop=FALSE]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",n)),c("Double Integrals",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("\u2119(t,y)",y),cbind(t,doubleintegrals))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(type < -1) { message("Types available for Double Integral: -1,0,1,2.  Showing: type=-1") }
      else if(type > 2) { message("Types available for Double Integral: -1,0,1,2.  Showing: type=2")}
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>y</i>=",esym,format(psi,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Double Integrals" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_DoubleIntegral2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>y</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>y</i><br>" }
        if(is.null(yaxis)) { yaxis <- "\u2119(<i>t,y</i>|<i>s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(psi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
        doubleintegralline <- list(color=red$d,width=4)
        doubleintegralmarker <- list(color=red$d,width=2,symbol="square")
        markeropacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_DoubleIntegral2D")
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=y,y=doubleintegrals[i,],mode="lines",line=doubleintegralline,opacity=markeropacity,hoverinfo="x+y")
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter",x=y,y=doubleintegrals[i,],mode="markers",marker=doubleintegralmarker,opacity=markeropacity,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_DoubleIntegral3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "\u2119(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","\u2119(<i>t,y</i>)=",format(doubleintegrals[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        if(psi > 0) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(psi == 0) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        PPmax <- max(doubleintegrals[1,1],doubleintegrals[1,n])
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,range=c(-0.03*PPmax,1.03*PPmax),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        doubleintegralline <- list(color=red$e,width=8)
        doubleintegralmarker <- list(color=red$e,size=4,symbol="square")
        gradient <- list(c(0,red$c),c(1,red$c))
        markeropacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_DoubleIntegral3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=y,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=densities,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE)
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=doubleintegrals[i,],name="\u2119(<i>t,y</i>)",mode="lines",line=doubleintegralline,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="PPt",showlegend=TRUE)
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter3d",x=y,y=rep(t[i],n),z=doubleintegrals[i,],mode="markers",marker=doubleintegralmarker,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="PPt",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=doubleintegrals,name="\u2119(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot options
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param sbeg    begin value for time axis
    #' @param send    end value for time axis
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotOption = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,sbeg=NULL,send=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      y <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(phi <= 0) { options <- private$OOneg } #protect against recursive call
      else { options <- private$OOpos }
      if(is.null(options)) { options <- self$Option(who="MC")[[1]] }
      densities <- private$o
      m <- length(s)
      n <- length(x)
      t <- s[1]
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m || Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        options <- options[Ixend:Ixbeg,,drop=FALSE]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend,drop=FALSE]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",n)),c("Options",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("phi",phi,rep("",n-1)),c("paths",paths,rep("",n-1)),c("skip",skip,rep("",n-1)),c("seed",seed,rep("",n-1)),c("method",method,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Options" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_Option2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        if(is.null(yaxis)) { yaxis <- "\uD835\uDD46(<i>s,x</i>|<i>t,y</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(phi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        optionline <- list(color=red$d,width=4)
        optionmarker <- list(color=red$d,width=2,symbol="square")
        markeropacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_Option2D")
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=options[i,],mode="lines",line=optionline,opacity=markeropacity,hoverinfo="x+y")
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter",x=x,y=options[i,],mode="markers",marker=optionmarker,opacity=markeropacity,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_Option3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis)) { zaxis <- "\uD835\uDD46(<i>s,x</i>|<i>t,y</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
          }
        }
        if(phi > 0) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        OOmax <- max(options[1,1],options[1,n])
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*x[1]-0.03*x[n],1.03*x[n]-0.03*x[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*s[m]-0.03*s[1],1.03*s[1]-0.03*s[m]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,range=c(-0.03*OOmax,1.03*OOmax),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        optionline <- list(color=red$e,width=8)
        optionmarker <- list(color=red$e,size=4,symbol="square")
        gradient <- list(c(0,red$c),c(1,red$c))
        markeropacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_Option3DSurface")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=x,y=s,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=densities,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE)
        skip <- as.integer((m-1)/10)
        if(skip < 1) { skip <- 1 }
        i <- 1
        fig <- add_trace(fig,type="scatter3d",x=x,y=rep(s[i],n),z=options[i,],name="\uD835\uDD46(<i>s,x</i>)",mode="lines",line=optionline,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="OOt",showlegend=TRUE)
        while(i < m)
        {
          i <- i+skip
          markeropacity <- markeropacity-0.07
          fig <- add_trace(fig,type="scatter3d",x=x,y=rep(s[i],n),z=options[i,],mode="markers",marker=optionmarker,opacity=markeropacity,hoverinfo="text",text=coordinates[i,],legendgroup="OOt",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot visiting time mode, median and mean
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotVisitingTimeModeMedianMean = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,0,NULL)
      self$set_plot_args(NULL,ptmax,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      blu <- private$plot_colors$blu
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      vtmmm <- private$vtmmm #protect against recursive call
      if(is.null(vtmmm)) { vtmmm <- self$VisitingTimeModeMedianMean(who="MC")[[1]] }
      pv <- private$pv #no copy or plot
      if(is.null(pv)) { pv <- self$VisitingTimeDensity(who="MC")[[1]] }
      Pv <- private$Pv #no copy or plot
      if(is.null(Pv)) { Pv <- self$VisitingTimeProbability(who="MC")[[1]] }
      m <- length(t)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pv <- pv[Ixbeg:Ixend]
        Pv <- Pv[Ixbeg:Ixend]
        m <- length(t)
      }
      pvmin <- min(pv)
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",2)),c("Visiting Time Mode, Median and Mean",rep("",2)),c("k",k,""),c("s",s,""),c("x",x,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("tv","pv","Pv"),vtmmm)
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Visiting Time Mode, Median and Mean"  }
      }
      else if(is.null(title)) { title <- ""  }
      lookup <- list(text=title,yref="container",y=0.95)
      if(labels == TRUE)
      {
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      densitymarker <- list(color=blu$d,size=5,symbol="square")
      probabilitymarker <- list(color=grn$d,size=5,symbol="square")
      meandashline <- list(color=cyn$d,dash="longdash",width=3)
      meandotline <- list(color=cyn$d,dash="dot",width=1)
      mediandashline <- list(color=grn$d,dash="dash",width=3)
      mediandotline <- list(color=grn$d,dash="dot",width=1)
      modedashline <- list(color=blu$d,dash="dot",width=3)
      modedotline <- list(color=blu$d,dash="dot",width=1)
      fig <- plot_ly()
      # OUP_MC_VisitingTimeModeMedianMean2DDensity
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else if(pvmin < 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(pvmin,ptmax)) }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,ptmax)) }
        legendpos <- list(x=1.0,y=0.9,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimeModeMedianMean2DDensity")
        fig <- add_trace(fig,type="scatter",x=c(vtmmm[3,1],vtmmm[3,1]),y=c(0,vtmmm[3,2]),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[3,1]),y=c(vtmmm[3,2],vtmmm[3,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtmmm[2,1],vtmmm[2,1]),y=c(0,vtmmm[2,2]),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[2,1]),y=c(vtmmm[2,2],vtmmm[2,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtmmm[1,1],vtmmm[1,1]),y=c(0,vtmmm[1,2]),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[1,1]),y=c(vtmmm[1,2],vtmmm[1,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=pv,name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=densitymarker,hoverinfo="x+y")
      }
      # OUP_MC_VisitingTimeModeMedianMean2DProbability
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimeModeMedianMean2DProbability")
        fig <- add_trace(fig,type="scatter",x=c(vtmmm[3,1],vtmmm[3,1]),y=c(0,vtmmm[3,3]),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[3,1]),y=c(vtmmm[3,3],vtmmm[3,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtmmm[2,1],vtmmm[2,1]),y=c(0,vtmmm[2,3]),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[2,1]),y=c(vtmmm[2,3],vtmmm[2,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtmmm[1,1],vtmmm[1,1]),y=c(0,vtmmm[1,3]),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtmmm[1,1]),y=c(vtmmm[1,3],vtmmm[1,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=Pv,name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=probabilitymarker,hoverinfo="x+y")
      }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot visiting time lower, median and upper percentiles
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotVisitingTimePercentiles = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,0,NULL)
      self$set_plot_args(NULL,ptmax,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      blu <- private$plot_colors$blu
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      vtpct <- private$vtpct #protect against recursive call
      if(is.null(vtpct)) { vtpct <- self$VisitingTimePercentiles(who="MC")[[1]] }
      pv <- private$pv #no copy or plot
      if(is.null(pv)) { pv <- self$VisitingTimeDensity(who="MC")[[1]] }
      Pv <- private$Pv #no copy or plot
      if(is.null(Pv)) { Pv <- self$VisitingTimeProbability(who="MC")[[1]] }
      m <- length(t)
      s <- t[1]
      if(Ppct > 0.5)
      {
        Ppctlower <- 1-Ppct
        Ppcthalf <- 0.5
        Ppctupper <- Ppct
      }
      else
      {
        Ppctlower<- Ppct
        Ppcthalf <- 0.5
        Ppctupper <- 1-Ppct
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pv <- pv[Ixbeg:Ixend]
        Pv <- Pv[Ixbeg:Ixend]
        m <- length(t)
      }
      pvmin <- min(pv)
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",2)),c("Visiting Time Percentiles",rep("",2)),c("k",k,""),c("s",s,""),c("x",x,""),c("P%",Ppct,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("t%","pv","Pv"),vtpct)
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",<i>P</i><sub>%</sub>",bsym,"=",esym,format(Ppct,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Visiting Time Percentiles"  }
      }
      else if(is.null(title)) { title <- ""  }
      lookup <- list(text=title,yref="container",y=0.95)
      if(labels == TRUE)
      {
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
      fig <- plot_ly()
      # OUP_MC_VisitingTimePercentiles2DDensity
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        densitymarker <- list(color=blu$d,size=5,symbol="square")
        upperdashline <- list(color=grn$c,dash="longdash",width=3)
        upperdotline <- list(color=grn$c,dash="dot",width=1)
        halfdashline <- list(color=grn$c,dash="dash",width=3)
        halfdotline <- list(color=grn$c,dash="dot",width=1)
        lowerdashline <- list(color=grn$c,dash="dot",width=3)
        lowerdotline <- list(color=grn$c,dash="dot",width=1)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else if(pvmin < 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(pvmin,ptmax)) }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,ptmax)) }
        legendpos <- list(x=1.0,y=0.9,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimePercentiles2DDensity")
        fig <- add_trace(fig,type="scatter",x=c(vtpct[3,1],vtpct[3,1]),y=c(0,vtpct[3,2]),name=paste(sep="","<i>t</i><sub>",format(Ppctupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[3,1]),y=c(vtpct[3,2],vtpct[3,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtpct[2,1],vtpct[2,1]),y=c(0,vtpct[2,2]),name=paste(sep="","<i>t</i><sub>",format(Ppcthalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[2,1]),y=c(vtpct[2,2],vtpct[2,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtpct[1,1],vtpct[1,1]),y=c(0,vtpct[1,2]),name=paste(sep="","<i>t</i><sub>",format(Ppctlower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[1,1]),y=c(vtpct[1,2],vtpct[1,2]),name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=pv,name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=densitymarker,hoverinfo="x+y")
      }
      # OUP_MC_VisitingTimePercentiles2DProbability
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        probabilitymarker <- list(color=grn$d,size=5,symbol="square")
        upperdashline <- list(color=grn$c,dash="longdash",width=3)
        upperdotline <- list(color=grn$c,dash="dot",width=1)
        halfdashline <- list(color=grn$c,dash="dash",width=3)
        halfdotline <- list(color=grn$c,dash="dot",width=1)
        lowerdashline <- list(color=grn$c,dash="dot",width=3)
        lowerdotline <- list(color=grn$c,dash="dot",width=1)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimePercentiles2DProbability")
        fig <- add_trace(fig,type="scatter",x=c(vtpct[3,1],vtpct[3,1]),y=c(0,vtpct[3,3]),name=paste(sep="","<i>t</i><sub>",format(Ppctupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[3,1]),y=c(vtpct[3,3],vtpct[3,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtpct[2,1],vtpct[2,1]),y=c(0,vtpct[2,3]),name=paste(sep="","<i>t</i><sub>",format(Ppcthalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[2,1]),y=c(vtpct[2,3],vtpct[2,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(vtpct[1,1],vtpct[1,1]),y=c(0,vtpct[1,3]),name=paste(sep="","<i>t</i><sub>",format(Ppctlower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,vtpct[1,1]),y=c(vtpct[1,3],vtpct[1,3]),name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=Pv,name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=probabilitymarker,hoverinfo="x+y")
      }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot visiting time densities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotVisitingTimeDensity = function(type=NULL,pmax=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,0,NULL)
      self$set_plot_args(pmax,ptmax,NULL,NULL,zbeg,zend)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      ptmax <- private$plot_args[[2]]
      zbeg <- private$plot_args[[5]]
      zend <- private$plot_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      pv <- private$pv #protect against recursive call
      if(is.null(pv)) { pv <- self$VisitingTimeDensity(who="MC")[[1]] }
      fheat <- private$ForwardHeat()
      heat <- fheat[[1]]
      z <- fheat[[2]]
      m <- length(t)
      n <- length(z)
      s <- t[1]
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pv <- pv[Ixbeg:Ixend]
        heat <- heat[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        heat <- heat[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      pvmin <- min(pv)
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("Visiting Time Density",""),c("k",k),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","pv"),cbind(t,pv))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Visiting Time Density" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_VisitingTimeDensity2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,rangemode="tozero",ticks="outside") }
        else if(pvmin < 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(pvmin,ptmax),ticks="outside") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(0,ptmax),ticks="outside") }
        densitymarker <- list(color=pv,colorscale=list(c(0,blu$b),c(1,blu$d)))
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimeDensity2D")
        fig <- plot_ly() %>%
          add_trace(.,type="bar",x=t,y=pv,name="<i>p<sub>v</sub></i>(<i>t</i>|<i>x</i>)",marker=densitymarker,orientation="v",hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_VisitingTimeDensity3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        coordinatek <- vector("character",m)
        for(i in 1:m) { coordinatek[i] <- paste(sep="","<i>p<sub>v</sub></i>(<i>t</i>)=",format(pv[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>k</i>=",k) }
        if(x < mu) { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(1.03*pvmin-0.03*ptmax,1.03*ptmax-0.03*pvmin),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        gradient <- list(c(0,blu$d),c(1,blu$b))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        pvmesh <- MeshCurtainChunky(rep(k,m),t,pv,rep(0,m))
        pvmarker <- list(size=4,color=blu$d)
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_VisitingTimeDensity3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=z,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=heat,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="mesh3d",x=pvmesh$xvertex,y=pvmesh$yvertex,z=pvmesh$zvertex,i=pvmesh$ivertex,j=pvmesh$jvertex,k=pvmesh$kvertex,intensity=pvmesh$zvertex,name="<i>p<sub>v</sub></i>(<i>t</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.9,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter3d",x=rep(k,m),y=t,z=pv,name="<i>p<sub>v</sub></i>(<i>t</i>)",mode="markers",marker=pvmarker,opacity=0.0,hoverinfo="text",text=coordinatek,showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot visiting time probabilities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotVisitingTimeProbability = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,0,NULL)
      self$set_plot_args(pmax,NULL,NULL,NULL,zbeg,zend)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      zbeg <- private$plot_args[[5]]
      zend <- private$plot_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      Pv <- private$Pv #protect against recursive call
      if(is.null(Pv)) { Pv <- self$VisitingTimeProbability(who="MC")[[1]] }
      fheat <- private$ForwardHeat()
      heat <- fheat[[1]]
      z <- fheat[[2]]
      m <- length(t)
      n <- length(z)
      s <- t[1]
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        Pv <- Pv[Ixbeg:Ixend]
        heat <- heat[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        heat <- heat[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("Visiting Time Probability",""),c("k",k),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","Pv"),cbind(t,Pv))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "Visiting Time Probability" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_VisitingTimeProbability2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(0,1),ticks="outside")
        probabilitymarker <- list(color=Pv,colorscale=list(c(0,grn$b),c(1,grn$d)))
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_VisitingTimeProbability2D")
        fig <- plot_ly() %>%
          add_trace(.,type="bar",x=t,y=Pv,name="<i>P<sub>v</sub></i>(<i>t</i>|<i>x</i>)",marker=probabilitymarker,orientation="v",hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_VisitingTimeProbability3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P<sub>v</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        coordinatek <- vector("character",m)
        for(i in 1:m) { coordinatek[i] <- paste(sep="","<i>P<sub>v</sub></i>(<i>t</i>)=",format(Pv[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>k</i>=",k) }
        if(x < mu) { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        gradient <- list(c(0,grn$d),c(1,grn$b))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        Pvmesh <- MeshCurtainChunky(rep(k,m),t,Pv,rep(0,m))
        Pvmarker <- list(size=4,color=grn$d)
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_VisitingTimeProbability3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=z,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=heat,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="mesh3d",x=Pvmesh$xvertex,y=Pvmesh$yvertex,z=Pvmesh$zvertex,i=Pvmesh$ivertex,j=Pvmesh$jvertex,k=Pvmesh$kvertex,intensity=Pvmesh$zvertex,name="<i>P<sub>v</sub></i>(<i>t</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.9,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter3d",x=rep(k,m),y=t,z=Pv,name="<i>P<sub>v</sub></i>(<i>t</i>)",mode="markers",marker=Pvmarker,opacity=0.0,hoverinfo="text",text=coordinatek,showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot first passage time mode, median and mean
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotFirstPassageTimeModeMedianMean = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,1,NULL)
      self$set_plot_args(NULL,ptmax,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      blu <- private$plot_colors$blu
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      fptmmm <- private$fptmmm #protect against recursive call
      if(is.null(fptmmm)) { fptmmm <- self$FirstPassageTimeModeMedianMean(who="MC")[[1]] }
      pf <- private$pf #no copy or plot
      if(is.null(pf)) { pf <- self$FirstPassageTimeDensity(who="MC")[[1]] }
      Pf <- private$Pf #no copy or plot
      if(is.null(Pf)) { Pf <- self$FirstPassageTimeProbability(who="MC")[[1]] }
      m <- length(t)
      s <- t[1]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pf <- pf[Ixbeg:Ixend]
        Pf <- Pf[Ixbeg:Ixend]
        m <- length(t)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",2)),c("First Passage Time Mode, Median and Mean",rep("",2)),c("k",k,""),c("s",s,""),c("x",x,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("tf","pf","Pf"),fptmmm)
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "First Passage Time Mode, Median and Mean"  }
      }
      else if(is.null(title)) { title <- ""  }
      lookup <- list(text=title,yref="container",y=0.95)
      if(labels == TRUE)
      {
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      densitymarker <- list(color=blu$d,size=5,symbol="square")
      probabilitymarker <- list(color=grn$d,size=5,symbol="square")
      meandashline <- list(color=cyn$d,dash="longdash",width=3)
      meandotline <- list(color=cyn$d,dash="dot",width=1)
      mediandashline <- list(color=grn$d,dash="dash",width=3)
      mediandotline <- list(color=grn$d,dash="dot",width=1)
      modedashline <- list(color=blu$d,dash="dot",width=3)
      modedotline <- list(color=blu$d,dash="dot",width=1)
      fig <- plot_ly()
      # OUP_MC_FirstPassageTimeModeMedianMean2DDensity
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,ptmax)) }
        legendpos <- list(x=1.0,y=0.9,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimeModeMedianMean2DDensity")
        fig <- add_trace(fig,type="scatter",x=c(fptmmm[3,1],fptmmm[3,1]),y=c(0,fptmmm[3,2]),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[3,1]),y=c(fptmmm[3,2],fptmmm[3,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptmmm[2,1],fptmmm[2,1]),y=c(0,fptmmm[2,2]),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[2,1]),y=c(fptmmm[2,2],fptmmm[2,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptmmm[1,1],fptmmm[1,1]),y=c(0,fptmmm[1,2]),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[1,1]),y=c(fptmmm[1,2],fptmmm[1,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=pf,name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=densitymarker,hoverinfo="x+y")
      }
      # OUP_MC_FirstPassageTimeModeMedianMean2DProbability
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimeModeMedianMean2DProbability")
        fig <- add_trace(fig,type="scatter",x=c(fptmmm[3,1],fptmmm[3,1]),y=c(0,fptmmm[3,3]),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[3,1]),y=c(fptmmm[3,3],fptmmm[3,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptmmm[2,1],fptmmm[2,1]),y=c(0,fptmmm[2,3]),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[2,1]),y=c(fptmmm[2,3],fptmmm[2,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptmmm[1,1],fptmmm[1,1]),y=c(0,fptmmm[1,3]),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptmmm[1,1]),y=c(fptmmm[1,3],fptmmm[1,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=Pf,name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=probabilitymarker,hoverinfo="x+y")
      }
      fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot first passage time lower, median and upper percentiles
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotFirstPassageTimePercentiles = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,1,NULL)
      self$set_plot_args(NULL,ptmax,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      Ppct <- private$t_stoch_args[[5]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      blu <- private$plot_colors$blu
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      fptpct <- private$fptpct #protect against recursive call
      if(is.null(fptpct)) { fptpct <- self$FirstPassageTimePercentiles(who="MC")[[1]] }
      pf <- private$pf #no copy or plot
      if(is.null(pf)) { pf <- self$FirstPassageTimeDensity(who="MC")[[1]] }
      Pf <- private$Pf #no copy or plot
      if(is.null(Pf)) { Pf <- self$FirstPassageTimeProbability(who="MC")[[1]] }
      m <- length(t)
      s <- t[1]
      if(Ppct > 0.5)
      {
        Ppctlower <- 1-Ppct
        Ppcthalf <- 0.5
        Ppctupper <- Ppct
      }
      else
      {
        Ppctlower<- Ppct
        Ppcthalf <- 0.5
        Ppctupper <- 1-Ppct
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pf <- pf[Ixbeg:Ixend]
        Pf <- Pf[Ixbeg:Ixend]
        m <- length(t)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",rep("",2)),c("First Passage Time Percentiles",rep("",2)),c("k",k,""),c("s",s,""),c("x",x,""),c("P%",Ppct,""),c("rho",rho,""),c("mu",mu,""),c("sigma",sigma,""),c("paths",paths,""),c("skip",skip,""),c("seed",seed,""),c("method",method,""),c("t%","pf","Pf"),fptpct)
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",<i>P</i><sub>%</sub>",bsym,"=",esym,format(Ppct,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "First Passage Time Percentiles"  }
      }
      else if(is.null(title)) { title <- ""  }
      lookup <- list(text=title,yref="container",y=0.95)
      if(labels == TRUE)
      {
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
      fig <- plot_ly()
      # OUP_MC_FirstPassageTimePercentiles2DDensity
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        densitymarker <- list(color=blu$d,size=5,symbol="square")
        upperdashline <- list(color=grn$c,dash="longdash",width=3)
        upperdotline <- list(color=grn$c,dash="dot",width=1)
        halfdashline <- list(color=grn$c,dash="dash",width=3)
        halfdotline <- list(color=grn$c,dash="dot",width=1)
        lowerdashline <- list(color=grn$c,dash="dot",width=3)
        lowerdotline <- list(color=grn$c,dash="dot",width=1)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,ptmax)) }
        legendpos <- list(x=1.0,y=0.9,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimePercentiles2DDensity")
        fig <- add_trace(fig,type="scatter",x=c(fptpct[3,1],fptpct[3,1]),y=c(0,fptpct[3,2]),name=paste(sep="","<i>t</i><sub>",format(Ppctupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[3,1]),y=c(fptpct[3,2],fptpct[3,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptpct[2,1],fptpct[2,1]),y=c(0,fptpct[2,2]),name=paste(sep="","<i>t</i><sub>",format(Ppcthalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[2,1]),y=c(fptpct[2,2],fptpct[2,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptpct[1,1],fptpct[1,1]),y=c(0,fptpct[1,2]),name=paste(sep="","<i>t</i><sub>",format(Ppctlower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[1,1]),y=c(fptpct[1,2],fptpct[1,2]),name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=pf,name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=densitymarker,hoverinfo="x+y")
      }
      # OUP_MC_FirstPassageTimePercentiles2DProbability
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        probabilitymarker <- list(color=grn$d,size=5,symbol="square")
        upperdashline <- list(color=grn$c,dash="longdash",width=3)
        upperdotline <- list(color=grn$c,dash="dot",width=1)
        halfdashline <- list(color=grn$c,dash="dash",width=3)
        halfdotline <- list(color=grn$c,dash="dot",width=1)
        lowerdashline <- list(color=grn$c,dash="dot",width=3)
        lowerdotline <- list(color=grn$c,dash="dot",width=1)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(t[1],t[m]))
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimePercentiles2DProbability")
        fig <- add_trace(fig,type="scatter",x=c(fptpct[3,1],fptpct[3,1]),y=c(0,fptpct[3,3]),name=paste(sep="","<i>t</i><sub>",format(Ppctupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[3,1]),y=c(fptpct[3,3],fptpct[3,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptpct[2,1],fptpct[2,1]),y=c(0,fptpct[2,3]),name=paste(sep="","<i>t</i><sub>",format(Ppcthalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[2,1]),y=c(fptpct[2,3],fptpct[2,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(fptpct[1,1],fptpct[1,1]),y=c(0,fptpct[1,3]),name=paste(sep="","<i>t</i><sub>",format(Ppctlower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=c(s,fptpct[1,1]),y=c(fptpct[1,3],fptpct[1,3]),name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=t,y=Pf,name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",mode="markers",marker=probabilitymarker,hoverinfo="x+y")
      }
      fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot first passage time densities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param ptmax   maximum visiting time and first passage time densities
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotFirstPassageTimeDensity = function(type=NULL,pmax=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,1,NULL)
      self$set_plot_args(pmax,ptmax,NULL,NULL,zbeg,zend)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      ptmax <- private$plot_args[[2]]
      zbeg <- private$plot_args[[5]]
      zend <- private$plot_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      pf <- private$pf #protect against recursive call
      if(is.null(pf)) { pf <- self$FirstPassageTimeDensity(who="MC")[[1]] }
      bheat <- private$BoundedHeat()
      heat <- bheat[[1]]
      z <- bheat[[2]]
      m <- length(t)
      n <- length(z)
      s <- t[1]
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        pf <- pf[Ixbeg:Ixend]
        heat <- heat[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        heat <- heat[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("First Passage Time Density",""),c("k",k),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","pf"),cbind(t,pf))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "First Passage Time Density" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_FirstPassageTimeDensity2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,rangemode="tozero",ticks="outside") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(0,ptmax),ticks="outside") }
        densitymarker <- list(color=pf,colorscale=list(c(0,blu$b),c(1,blu$d)))
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimeDensity2D")
        fig <- plot_ly() %>%
          add_trace(.,type="bar",x=t,y=pf,name="<i>p<sub>f</sub></i>(<i>t</i>|<i>x</i>)",marker=densitymarker,orientation="v",hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_FirstPassageTimeDensity3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        coordinatek <- vector("character",m)
        for(i in 1:m) { coordinatek[i] <- paste(sep="","<i>p<sub>f</sub></i>(<i>t</i>)=",format(pf[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>k</i>=",k) }
        if(x < mu) { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(-0.03*ptmax,1.03*ptmax),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        gradient <- list(c(0,blu$e),c(1,blu$b))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        pfmesh <- MeshCurtainChunky(rep(k,m),t,pf,rep(0,m))
        pfmarker <- list(size=4,color=blu$d)
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_FirstPassageTimeDensity3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=z,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=heat,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="mesh3d",x=pfmesh$xvertex,y=pfmesh$yvertex,z=pfmesh$zvertex,i=pfmesh$ivertex,j=pfmesh$jvertex,k=pfmesh$kvertex,intensity=pfmesh$zvertex,name="<i>p<sub>f</sub></i>(<i>t</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.9,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter3d",x=rep(k,m),y=t,z=pf,name="<i>p<sub>f</sub></i>(<i>t</i>)",mode="markers",marker=pfmarker,opacity=0.0,hoverinfo="text",text=coordinatek,showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot first passage time probabilities
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
    #' @param pmax    maximum transition density
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotFirstPassageTimeProbability = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_t_stoch_args(NULL,NULL,NULL,1,NULL)
      self$set_plot_args(pmax,NULL,NULL,NULL,zbeg,zend)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      x <- private$t_stoch_args[[3]]
      paths <- private$path_args[[1]]
      skip <- private$path_args[[2]]
      seed <- private$path_args[[3]]
      method <- private$path_args[[4]]
      pmax <- private$plot_args[[1]]
      zbeg <- private$plot_args[[5]]
      zend <- private$plot_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      Pf <- private$Pf #protect against recursive call
      if(is.null(Pf)) { Pf <- self$FirstPassageTimeProbability(who="MC")[[1]] }
      fheat <- private$BoundedHeat()
      heat <- fheat[[1]]
      z <- fheat[[2]]
      m <- length(t)
      n <- length(z)
      s <- t[1]
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        Pf <- Pf[Ixbeg:Ixend]
        heat <- heat[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        heat <- heat[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Monte Carlo",""),c("First Passage Time Probability",""),c("k",k),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("sigma",sigma),c("paths",paths),c("skip",skip),c("seed",seed),c("method",method),c("t","Pf"),cbind(t,Pf))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",paths",bsym,"=",esym,format(paths,big.mark=","),")",esml)
        if(is.null(title)) { title <- "First Passage Time Probability" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_MC_FirstPassageTimeProbability2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(0,1),ticks="outside")
        probabilitymarker <- list(color=Pf,colorscale=list(c(0,grn$b),c(1,grn$d)))
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_MC_FirstPassageTimeProbability2D")
        fig <- plot_ly() %>%
          add_trace(.,type="bar",x=t,y=Pf,name="<i>P<sub>f</sub></i>(<i>t</i>|<i>x</i>)",marker=probabilitymarker,orientation="v",hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_MC_FirstPassageTimeProbability3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P<sub>f</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        coordinatek <- vector("character",m)
        for(i in 1:m) { coordinatek[i] <- paste(sep="","<i>P<sub>f</sub></i>(<i>t</i>)=",format(Pf[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>k</i>=",k) }
        if(x < mu) { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        gradient <- list(c(0,grn$e),c(1,grn$b))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        if(reverse)
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$a),c(0.01,red$a),c(0.02,red$b),c(0.03,red$c),c(0.04,red$d),c(0.05,red$e),c(1,red$e)) }
          else { heatgradient <- list(c(0,gry$a),c(0.1*pmax,red$a),c(0.2*pmax,red$b),c(0.3*pmax,red$c),c(0.4*pmax,red$d),c(0.5*pmax,red$e),c(1,red$e)) }
        }
        else
        {
          if(is.nan(pmax)) { heatgradient <- list(c(0,gry$e),c(0.01,red$e),c(0.02,red$d),c(0.03,red$c),c(0.04,red$b),c(0.05,red$a),c(1,red$a)) }
          else { heatgradient <- list(c(0,gry$e),c(0.1*pmax,red$e),c(0.2*pmax,red$d),c(0.3*pmax,red$c),c(0.4*pmax,red$b),c(0.5*pmax,red$a),c(1,red$a)) }
        }
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        Pfmesh <- MeshCurtainChunky(rep(k,m),t,Pf,rep(0,m))
        Pfmarker <- list(size=4,color=grn$d)
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_MC_FirstPassageTimeProbability3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=z,y=t,z=matrix(0.0,m,n),name="paths",showscale=FALSE,lighting=shine,lightposition=rise,surfacecolor=heat,colorscale=heatgradient,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="mesh3d",x=Pfmesh$xvertex,y=Pfmesh$yvertex,z=Pfmesh$zvertex,i=Pfmesh$ivertex,j=Pfmesh$jvertex,k=Pfmesh$kvertex,intensity=Pfmesh$zvertex,name="<i>P<sub>f</sub></i>(<i>t</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.9,hoverinfo="skip",showlegend=TRUE) %>%
          add_trace(.,type="scatter3d",x=rep(k,m),y=t,z=Pf,name="<i>P<sub>f</sub></i>(<i>t</i>)",mode="markers",marker=Pfmarker,opacity=0.0,hoverinfo="text",text=coordinatek,showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    OUP = NULL,
    # private attributes ----
    oup_params = NULL,
    y_stoch_args = NULL,
    x_stoch_args = NULL,
    t_stoch_args = NULL,
    path_args = NULL,
    plot_args = NULL,
    undo_args = NULL,
    plot_info = NULL,
    flags = NULL,
    # private output fields ----
    yforward = NULL,
    xbackward = NULL,
    tforward = NULL,
    tbounded = NULL,
    G = NULL,
    H2 = NULL,
    p = NULL,
    Pneg = NULL,
    Ppos = NULL,
    PPneg = NULL,
    PPpos = NULL,
    o = NULL,
    Oneg = NULL,
    Opos = NULL,
    OOneg = NULL,
    OOpos = NULL,
    vtmmm = NULL,
    vtpct = NULL,
    pv = NULL,
    Pv = NULL,
    fptmmm = NULL,
    fptpct = NULL,
    pf = NULL,
    Pf = NULL,
    # private globals ----
    undoIx = NULL,
    syncyxt = NULL,
    forwardyt = NULL,
    ystdnorm = NULL,
    xstdnorm = NULL,
    tstdnorm = NULL,
    tfpt = NULL,
    fheat = NULL,
    fz = NULL,
    bheat = NULL,
    bz = NULL,
    plot_colors = NULL,
    modebar_2D = NULL,
    modebar_3D = NULL,
    plot_types = NULL,
    # private colors ----
    rainbow = function(name,opaque)
    {
      if(name == "light")
      {
        red <- list(a="rgb(231,184,184)",b="rgb(211,124,124)",c="rgb(192, 64, 64)",d="rgb(131, 44, 44)",e="rgb( 71, 24, 24)")
        ylw <- list(a="rgb(234,234,161)",b="rgb(219,219, 96)",c="rgb(191,191, 44)",d="rgb(126,126, 29)",e="rgb( 61, 61, 14)")
        grn <- list(a="rgb(194,221,194)",b="rgb(142,193,142)",c="rgb( 90,165, 90)",d="rgb( 62,113, 62)",e="rgb( 34, 61, 34)")
        cyn <- list(a="rgb(178,237,237)",b="rgb(113,222,222)",c="rgb( 47,207,207)",d="rgb( 33,142,142)",e="rgb( 18, 77, 77)")
        blu <- list(a="rgb(227,227,248)",b="rgb(161,161,234)",c="rgb( 96, 96,219)",d="rgb( 44, 44,191)",e="rgb( 29, 29,126)")
        mgn <- list(a="rgb(237,178,237)",b="rgb(222,113,222)",c="rgb(207, 47,207)",d="rgb(142, 33,142)",e="rgb( 77, 18, 77)")
        gry <- list(a="rgb(255,255,255)",b="rgb(192,192,192)",c="rgb(128,128,128)",d="rgb( 64, 64, 64)",e="rgb(  0,  0,  0)")
        background <- paste(sep="","rgba(255,255,255,",opaque,")")
        font <- "rgb(0,0,0)"
        reverse <- FALSE
      }
      else if(name == "dark")
      {
        red <- list(e="rgb(231,184,184)",d="rgb(211,124,124)",c="rgb(192, 64, 64)",b="rgb(131, 44, 44)",a="rgb( 71, 24, 24)")
        ylw <- list(e="rgb(234,234,161)",d="rgb(219,219, 96)",c="rgb(191,191, 44)",b="rgb(126,126, 29)",a="rgb( 61, 61, 14)")
        grn <- list(e="rgb(194,221,194)",d="rgb(142,193,142)",c="rgb( 90,165, 90)",b="rgb( 62,113, 62)",a="rgb( 34, 61, 34)")
        cyn <- list(e="rgb(178,237,237)",d="rgb(113,222,222)",c="rgb( 47,207,207)",b="rgb( 33,142,142)",a="rgb( 18, 77, 77)")
        blu <- list(e="rgb(227,227,248)",d="rgb(161,161,234)",c="rgb( 96, 96,219)",b="rgb( 44, 44,191)",a="rgb( 29, 29,126)")
        mgn <- list(e="rgb(237,178,237)",d="rgb(222,113,222)",c="rgb(207, 47,207)",b="rgb(142, 33,142)",a="rgb( 77, 18, 77)")
        gry <- list(e="rgb(255,255,255)",d="rgb(192,192,192)",c="rgb(128,128,128)",b="rgb( 64, 64, 64)",a="rgb(  0,  0,  0)")
        background <- paste(sep="","rgba(0,0,0,",opaque,")")
        font <- "rgb(255,255,255)"
        reverse <- TRUE
      }
      return(list(red=red,ylw=ylw,grn=grn,cyn=cyn,blu=blu,mgn=mgn,gry=gry,background=background,font=font,reverse=reverse))
    },
    # private input methods ----
    extract_scalar = function(input)
    {
      sca <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        if(is.numeric(input[1]))
        {
          input[1] <- as.numeric(input[1])
          if(!any(is.na(input[1]))) { sca <- input[1] }
          else { message(paste(sep="",input[1]," is not a real number:")) }
        }
        else { message(paste(sep="",input[1]," is not a real number:")) }
      }
      return(sca)
    },
    extract_vector = function(input,updown=0)
    {
      vec <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        OK <- TRUE
        n <- length(input)
        i <- 0
        while(i < n && OK == TRUE)
        {
          i <- i+1
          if(is.numeric(input[i]))
          {
            input[i] <- as.numeric(input[i])
            if(any(is.na(input[i]))) { OK <- FALSE }
          }
          else { OK <- FALSE}
        }
        if(OK == TRUE)
        {
          if(updown > 0) { vec <- sort(input) }
          else if(updown == 0) { vec <- input }
          else { vec <- sort(input,decreasing=TRUE) }
        }
        else { message(paste(sep="","vector[",i,"]=",input[i]," is not a real number.")) }
      }
      return(vec)
    },
    extract_character = function(input)
    {
      chr <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        if(!is.numeric(input[1])) { chr <- as.character(input[1]) }
        else { message(paste(sep="",input[1]," is a number.")) }
      }
      return(chr)
    },
    extract_integer = function(input)
    {
      int <- NULL
      sca <- extract_scalar(input)
      if(!is.null(sca)) { int <- as.integer(sca) }
      return(int)
    },
    extract_boolean = function(input)
    {
      bool <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        if(input[1] == TRUE || input[1] == FALSE) { bool <- input[1] }
        else { message(paste("non-boolean input:",input)) }
      }
      return(bool)
    },
    extract_rgb = function(input)
    {
      clr <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        if(is.character(input))
        {
          rgb <- input[1]
          n <- str_length(rgb)
          if(n > 6)
          {
            if(str_sub(rgb,1,1) == "#")
            {
              red <- octal2decimal(str_sub(rgb,2,3))
              grn <- octal2decimal(str_sub(rgb,4,5))
              blu <- octal2decimal(str_sub(rgb,6,7))
              if(!is.null(red) && !is.null(grn) && !is.null(blu)) { clr <- paste(sep="","rgb(",red,",",grn,",",blu,")") }
              else { message("unrecognized hexadecimal color:",rgb)}
            }
            else if(n > 9)
            {
              if(tolower(str_sub(rgb,1,4)) == "rgb(" && str_sub(rgb,n,n) == ")")
              {
                OK <- TRUE
                hit <- FALSE
                k <- 4
                j <- k
                while(OK == TRUE && hit == FALSE && k < n)
                {
                  k <- k+1
                  token <- str_sub(rgb,k,k)
                  if(token == "," || token == ")")
                  {
                    hit <- TRUE
                    if(k-j == 1) { OK <- FALSE }
                  }
                  else if(token != "0" && token != "1" && token != "2" && token != "3" && token != "4" && token != "5" && token != "6" && token != "7" && token != "8" && token != "9") { OK <- FALSE }
                }
                if(OK == TRUE)
                {
                  red <- as.integer(str_sub(rgb,j+1,k-1))
                  if(red >= 0 && red <= 255)
                  {
                    hit <- FALSE
                    j <- k
                    while(OK == TRUE && hit == FALSE && k < n)
                    {
                      k <- k+1
                      token <- str_sub(rgb,k,k)
                      if(token == "," || token == ")")
                      {
                        hit <- TRUE
                        if(k-j == 1) { OK <- FALSE }
                      }
                      else if(token != "0" && token != "1" && token != "2" && token != "3" && token != "4" && token != "5" && token != "6" && token != "7" && token != "8" && token != "9") { OK <- FALSE }
                    }
                    if(OK == TRUE)
                    {
                      grn <- as.integer(str_sub(rgb,j+1,k-1))
                      if(grn >= 0 && grn <= 255)
                      {
                        hit <- FALSE
                        j <- k
                        while(OK == TRUE && hit == FALSE && k < n)
                        {
                          k <- k+1
                          token <- str_sub(rgb,k,k)
                          if(token == "," || token == ")")
                          {
                            hit <- TRUE
                            if(k-j == 1) { OK <- FALSE }
                          }
                          else if(token != "0" && token != "1" && token != "2" && token != "3" && token != "4" && token != "5" && token != "6" && token != "7" && token != "8" && token != "9") { OK <- FALSE }
                        }
                        if(OK == TRUE)
                        {
                          blu <- as.integer(str_sub(rgb,j+1,k-1))
                          if(blu >= 0 && blu <= 255) { clr <- paste(sep="","rgb(",red,",",grn,",",blu,")") }
                          else { message(paste("blu in rgb color must be integer from 0 to 255:",rgb)) }
                        }
                        else { message(paste("blu in rgb color must be integer from 0 to 255:",rgb)) }
                      }
                      else { message(paste("grn in rgb color must be integer from 0 to 255:",rgb)) }
                    }
                    else { message(paste("grn in rgb color must be integer from 0 to 255:",rgb)) }
                  }
                  else { message(paste("red in rgb color must be integer from 0 to 255:",rgb)) }
                }
                else { message(paste("red in rgb color must be integer from 0 to 255:",rgb)) }
              }
              else { message(paste("unrecognizable rgb color:",rgb)) }
            }
            else { message(paste("not enough characters for rgb color:",rgb)) }
          }
          else { message(paste("not enough characters for hexadecimal or rgb color:",rgb)) }
        }
        else { message(paste("non-character input:",input)) }
      }
      return(clr)
    },
    lists_equal = function(list1,list2)
    {
      allequal <- TRUE
      n1 <- length(list1)
      n2 <- length(list2)
      if(n1 == n2)
      {
        j <- 0
        while(j < n1 && allequal == TRUE )
        {
          j <- j+1
          allequal <- vecs_equal(list1[[j]],list2[[j]])
        }
      }
      else { allequal <- FALSE }
      return(allequal)
    },
    vecs_equal = function(vec1,vec2)
    {
      tolerance <- 1e-12
      allequal <- TRUE
      n1 <- length(vec1)
      n2 <- length(vec2)
      if(n1 == n2)
      {
        j <- 0
        while(j < n1 && allequal == TRUE )
        {
          j <- j+1
          if(abs(vec1[j]-vec2[j]) > tolerance) { allequal <- FALSE }
        }
      }
      else { allequal <- FALSE }
      return(allequal)
    },
    octal2decimal = function(octal)
    {
      dec <- NULL
      n <- str_length(octal)
      if(n > 0)
      {
        byte <- toupper(str_sub(octal,n,n))
        small <- 16
        if(byte == "0") { small <- 0 }
        else if(byte == "1") { small <- 1 }
        else if(byte == "2") { small <- 2 }
        else if(byte == "3") { small <- 3 }
        else if(byte == "4") { small <- 4 }
        else if(byte == "5") { small <- 5 }
        else if(byte == "6") { small <- 6 }
        else if(byte == "7") { small <- 7 }
        else if(byte == "8") { small <- 8 }
        else if(byte == "9") { small <- 9 }
        else if(byte == "A") { small <- 10 }
        else if(byte == "B") { small <- 11 }
        else if(byte == "C") { small <- 12 }
        else if(byte == "D") { small <- 13 }
        else if(byte == "E") { small <- 14 }
        else if(byte == "F") { small <- 15 }
        if(small < 16)
        {
          if(n > 1)
          {
            byte <- toupper(str_sub(octal,n-1,n-1))
            big <- 16
            if(byte == "0") { big <- 0 }
            else if(byte == "1") { big <- 1 }
            else if(byte == "2") { big <- 2 }
            else if(byte == "3") { big <- 3 }
            else if(byte == "4") { big <- 4 }
            else if(byte == "5") { big <- 5 }
            else if(byte == "6") { big <- 6 }
            else if(byte == "7") { big <- 7 }
            else if(byte == "8") { big <- 8 }
            else if(byte == "9") { big <- 9 }
            else if(byte == "A") { big <- 10 }
            else if(byte == "B") { big <- 11 }
            else if(byte == "C") { big <- 12 }
            else if(byte == "D") { big <- 13 }
            else if(byte == "E") { big <- 14 }
            else if(byte == "F") { big <- 15 }
            if(big < 16) { dec <- 16*big+small }
          }
          else { dec <- small }
        }
      }
      return(dec)
    },
    # private plot methods ----
    MeshCurtainChunky = function(x,y,ztop,zbottom)
    {
      m <- min(length(x),length(y),length(ztop),length(zbottom))
      xv <- vector("double",5*m)
      yv <- vector("double",5*m)
      zv <- vector("double",5*m)
      iv <- vector("double",4*m)
      jv <- vector("double",4*m)
      kv <- vector("double",4*m)
      # coordinates
      xv[1] <- x[1]
      xv[2] <- 0.5*(x[1]+x[2])
      xv[3] <- 0.5*(x[1]+x[2])
      xv[4] <- x[1]
      xv[5] <- 0.75*x[1]+0.25*x[2]
      yv[1] <- y[1]
      yv[2] <- 0.5*(y[1]+y[2])
      yv[3] <- 0.5*(y[1]+y[2])
      yv[4] <- y[1]
      yv[5] <- 0.75*y[1]+0.25*y[2]
      zv[1] <- zbottom[1]
      zv[2] <- zbottom[1]
      zv[3] <- ztop[1]
      zv[4] <- ztop[1]
      zv[5] <- 0.5*(ztop[1]+zbottom[1])
      i <- 1
      while(i < m-1)
      {
        i <- i+1
        xv[5*i-4] <- 0.5*(x[i-1]+x[i])
        xv[5*i-3] <- 0.5*(x[i]+x[i+1])
        xv[5*i-2] <- 0.5*(x[i]+x[i+1])
        xv[5*i-1] <- 0.5*(x[i-1]+x[i])
        xv[5*i] <- x[i]
        yv[5*i-4] <- 0.5*(y[i-1]+y[i])
        yv[5*i-3] <- 0.5*(y[i]+y[i+1])
        yv[5*i-2] <- 0.5*(y[i]+y[i+1])
        yv[5*i-1] <- 0.5*(y[i-1]+y[i])
        yv[5*i] <- y[i]
        zv[5*i-4] <- zbottom[i]
        zv[5*i-3] <- zbottom[i]
        zv[5*i-2] <- ztop[i]
        zv[5*i-1] <- ztop[i]
        zv[5*i] <- 0.5*(ztop[i]+zbottom[i])
      }
      if(m > 1)
      {
        xv[5*m-4] <- 0.5*(x[m-1]+x[m])
        xv[5*m-3] <- x[m]
        xv[5*m-2] <- x[m]
        xv[5*m-1] <- 0.5*(x[m-1]+x[m])
        xv[5*m] <- 0.25*x[m-1]+0.75*x[m]
        yv[5*m-4] <- 0.5*(y[m-1]+y[m])
        yv[5*m-3] <- y[m]
        yv[5*m-2] <- y[m]
        yv[5*m-1] <- 0.5*(y[m-1]+y[m])
        yv[5*m] <- 0.25*y[m-1]+0.75*y[m]
        zv[5*m-4] <- zbottom[m]
        zv[5*m-3] <- zbottom[m]
        zv[5*m-2] <- ztop[m]
        zv[5*m-1] <- ztop[m]
        zv[5*m] <- 0.5*(ztop[m]+zbottom[m])
      }
      # triangles, right-hand side indexes are zero-based
      i <- 0
      while(i < m)
      {
        i <- i+1
        iv[4*i-3] <- 5*i-1
        iv[4*i-2] <- 5*i-1
        iv[4*i-1] <- 5*i-1
        iv[4*i] <- 5*i-1
        jv[4*i-3] <- 5*i-5
        jv[4*i-2] <- 5*i-4
        jv[4*i-1] <- 5*i-3
        jv[4*i] <- 5*i-2
        kv[4*i-3] <- 5*i-4
        kv[4*i-2] <- 5*i-3
        kv[4*i-1] <- 5*i-2
        kv[4*i] <- 5*i-5
      }
      return(list(xvertex=xv,yvertex=yv,zvertex=zv,ivertex=iv,jvertex=jv,kvertex=kv))
    },
    MeshWallChunky = function(xleft,xright,yfront,yback,ztop,zbottom)
    {
      m <- min(length(xleft),length(xright),length(yfront),length(yback),length(ztop),length(zbottom))
      xv <- vector("double",14*m)
      yv <- vector("double",14*m)
      zv <- vector("double",14*m)
      iv <- vector("double",24*m)
      jv <- vector("double",24*m)
      kv <- vector("double",24*m)
      i <- 0
      while(i < m)
      {
        i <- i+1
        # coordinates
        xv[14*i-13] <- xleft[i]
        xv[14*i-12] <- xright[i]
        xv[14*i-11] <- xleft[i]
        xv[14*i-10] <- xright[i]
        xv[14*i-9] <- xleft[i]
        xv[14*i-8] <- xright[i]
        xv[14*i-7] <- xleft[i]
        xv[14*i-6] <- xright[i]
        xv[14*i-5] <- 0.5*(xleft[i]+xright[i])
        xv[14*i-4] <- 0.5*(xleft[i]+xright[i])
        xv[14*i-3] <- xleft[i]
        xv[14*i-2] <- xright[i]
        xv[14*i-1] <- 0.5*(xleft[i]+xright[i])
        xv[14*i] <- 0.5*(xleft[i]+xright[i])
        yv[14*i-13] <- yfront[i]
        yv[14*i-12] <- yfront[i]
        yv[14*i-11] <- yback[i]
        yv[14*i-10] <- yback[i]
        yv[14*i-9] <- yfront[i]
        yv[14*i-8] <- yfront[i]
        yv[14*i-7] <- yback[i]
        yv[14*i-6] <- yback[i]
        yv[14*i-5] <- yfront[i]
        yv[14*i-4] <- yback[i]
        yv[14*i-3] <- 0.5*(yfront[i]+yback[i])
        yv[14*i-2] <-  0.5*(yfront[i]+yback[i])
        yv[14*i-1] <-  0.5*(yfront[i]+yback[i])
        yv[14*i] <-  0.5*(yfront[i]+yback[i])
        zv[14*i-13] <- zbottom[i]
        zv[14*i-12] <- zbottom[i]
        zv[14*i-11] <- zbottom[i]
        zv[14*i-10] <- zbottom[i]
        zv[14*i-9] <- ztop[i]
        zv[14*i-8] <- ztop[i]
        zv[14*i-7] <- ztop[i]
        zv[14*i-6] <- ztop[i]
        zv[14*i-5] <- 0.5*(zbottom[i]+ztop[i])
        zv[14*i-4] <- 0.5*(zbottom[i]+ztop[i])
        zv[14*i-3] <- 0.5*(zbottom[i]+ztop[i])
        zv[14*i-2] <- 0.5*(zbottom[i]+ztop[i])
        zv[14*i-1] <- zbottom[i]
        zv[14*i] <- ztop[i]
        # triangles, right-hand side indexes are zero-based
        iv[24*i-23] <- 8+14*(i-1)
        iv[24*i-22] <- 8+14*(i-1)
        iv[24*i-21] <- 8+14*(i-1)
        iv[24*i-20] <- 8+14*(i-1)
        iv[24*i-19] <- 9+14*(i-1)
        iv[24*i-18] <- 9+14*(i-1)
        iv[24*i-17] <- 9+14*(i-1)
        iv[24*i-16] <- 9+14*(i-1)
        iv[24*i-15] <- 10+14*(i-1)
        iv[24*i-14] <- 10+14*(i-1)
        iv[24*i-13] <- 10+14*(i-1)
        iv[24*i-12] <- 10+14*(i-1)
        iv[24*i-11] <- 11+14*(i-1)
        iv[24*i-10] <- 11+14*(i-1)
        iv[24*i-9] <- 11+14*(i-1)
        iv[24*i-8] <- 11+14*(i-1)
        iv[24*i-7] <- 12+14*(i-1)
        iv[24*i-6] <- 12+14*(i-1)
        iv[24*i-5] <- 12+14*(i-1)
        iv[24*i-4] <- 12+14*(i-1)
        iv[24*i-3] <- 13+14*(i-1)
        iv[24*i-2] <- 13+14*(i-1)
        iv[24*i-1] <- 13+14*(i-1)
        iv[24*i] <- 13+14*(i-1)
        jv[24*i-23] <- 0+14*(i-1)
        jv[24*i-22] <- 1+14*(i-1)
        jv[24*i-21] <- 5+14*(i-1)
        jv[24*i-20] <- 4+14*(i-1)
        jv[24*i-19] <- 2+14*(i-1)
        jv[24*i-18] <- 3+14*(i-1)
        jv[24*i-17] <- 7+14*(i-1)
        jv[24*i-16] <- 6+14*(i-1)
        jv[24*i-15] <- 0+14*(i-1)
        jv[24*i-14] <- 2+14*(i-1)
        jv[24*i-13] <- 6+14*(i-1)
        jv[24*i-12] <- 4+14*(i-1)
        jv[24*i-11] <- 1+14*(i-1)
        jv[24*i-10] <- 3+14*(i-1)
        jv[24*i-9] <- 7+14*(i-1)
        jv[24*i-8] <- 5+14*(i-1)
        jv[24*i-7] <- 0+14*(i-1)
        jv[24*i-6] <- 1+14*(i-1)
        jv[24*i-5] <- 3+14*(i-1)
        jv[24*i-4] <- 2+14*(i-1)
        jv[24*i-3] <- 4+14*(i-1)
        jv[24*i-2] <- 5+14*(i-1)
        jv[24*i-1] <- 7+14*(i-1)
        jv[24*i] <- 6+14*(i-1)
        kv[24*i-23] <- 1+14*(i-1)
        kv[24*i-22] <- 5+14*(i-1)
        kv[24*i-21] <- 4+14*(i-1)
        kv[24*i-20] <- 0+14*(i-1)
        kv[24*i-19] <- 3+14*(i-1)
        kv[24*i-18] <- 7+14*(i-1)
        kv[24*i-17] <- 6+14*(i-1)
        kv[24*i-16] <- 2+14*(i-1)
        kv[24*i-15] <- 2+14*(i-1)
        kv[24*i-14] <- 6+14*(i-1)
        kv[24*i-13] <- 4+14*(i-1)
        kv[24*i-12] <- 0+14*(i-1)
        kv[24*i-11] <- 3+14*(i-1)
        kv[24*i-10] <- 7+14*(i-1)
        kv[24*i-9] <- 5+14*(i-1)
        kv[24*i-8] <- 1+14*(i-1)
        kv[24*i-7] <- 1+14*(i-1)
        kv[24*i-6] <- 3+14*(i-1)
        kv[24*i-5] <- 2+14*(i-1)
        kv[24*i-4] <- 0+14*(i-1)
        kv[24*i-3] <- 5+14*(i-1)
        kv[24*i-2] <- 7+14*(i-1)
        kv[24*i-1] <- 6+14*(i-1)
        kv[24*i] <- 4+14*(i-1)
      }
      return(list(xvertex=xv,yvertex=yv,zvertex=zv,ivertex=iv,jvertex=jv,kvertex=kv))
    },
    ForwardHeat = function()
    {
      heat <- private$fheat
      z <- private$fz
      if(is.null(heat) || is.null(z))
      {
        rho <- private$oup_params[[1]]
        mu <- private$oup_params[[2]]
        sigma <- private$oup_params[[3]]
        t <- private$t_stoch_args[[1]]
        k <- private$t_stoch_args[[2]]
        x <- private$t_stoch_args[[3]]
        paths <- private$path_args[[1]]
        skip <- private$path_args[[2]]
        seed <- private$path_args[[3]]
        method <- private$path_args[[4]]
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        forward <- private$tforward
        if(is.null(forward))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { forward <- RcppOUPMCForwardPathRungeKutta(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          else { forward <- RcppOUPMCForwardPathIntegralEquation(stdnorm,x,m,skip,dt,rho,mu,sigma) }
          private$tforward <- forward
        }
        zdif <- 2*max(abs(x-mu),abs(x-k),abs(k-mu))
        if(zdif < 1) { zdif <- 1}
        zscale <- 1
        while(zdif > zscale) { zscale <- 10*zscale }
        zdif <- round(zdif/zscale,1)*zscale
        zby <- zdif/100
        zfrom <- (x+mu+k)/3-50*zby
        zto <- zfrom+100*zby
        z <- seq(from=zfrom,to=zto,by=zby)
        heat <- RcppOUPMCHeatCountZ(forward,z)
        private$fheat <- heat
        private$fz <- z
        private$t_stoch_args$z <- z
      }
      return(list(fheat=heat,fz=z))
    },
    BoundedHeat = function()
    {
      heat <- private$bheat
      z <- private$bz
      if(is.null(heat) || is.null(z))
      {
        rho <- private$oup_params[[1]]
        mu <- private$oup_params[[2]]
        sigma <- private$oup_params[[3]]
        t <- private$t_stoch_args[[1]]
        k <- private$t_stoch_args[[2]]
        x <- private$t_stoch_args[[3]]
        paths <- private$path_args[[1]]
        skip <- private$path_args[[2]]
        seed <- private$path_args[[3]]
        method <- private$path_args[[4]]
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        bounded <- private$tbounded
        if(is.null(bounded))
        {
          stdnorm <- private$tstdnorm
          if(is.null(stdnorm))
          {
            stdnorm <- RcppOUPMCStandardNormal(m,skip,paths,seed)
            private$tstdnorm <- stdnorm
          }
          if(method == 4) { bndfpt <- RcppOUPMCBoundedPathRungeKutta(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          else { bndfpt <- RcppOUPMCBoundedPathIntegralEquation(stdnorm,k,x,m,skip,dt,rho,mu,sigma) }
          bounded <- bndfpt[1:m,,drop=FALSE]
          fpt <- bndfpt[m+1,,drop=FALSE]
          private$tbounded <- bounded
          private$tfpt <- fpt
        }
        zdif <- 2*max(abs(x-mu),abs(x-k),abs(k-mu))
        if(zdif < 1) { zdif <- 1}
        zscale <- 1
        while(zdif > zscale) { zscale <- 10*zscale }
        zdif <- round(zdif/zscale,1)*zscale
        zby <- zdif/100
        zfrom <- (x+mu+k)/3-50*zby
        zto <- zfrom+100*zby
        z <- seq(from=zfrom,to=zto,by=zby)
        heat <- RcppOUPMCHeatCountZ(bounded,z)
        private$bheat <- heat
        private$bz <- z
        private$t_stoch_args$z <- z
      }
      return(list(bheat=heat,bz=z))
    },
    index = function(x,beg,end)
    {
      n <- length(x)
      Ixbeg <- 1
      Ixend <- n
      if(!is.null(beg))
      {
        if(beg == -Inf) { Ixbeg <- 1 }
        else if(beg == Inf) { Ixbeg <- n }
        else
        {
          sca <- private$extract_scalar(beg)
          if(!is.null(sca))
          {
            if(sca < x[1]) { Ixbeg <- 1 }
            else if(sca > x[n]) { Ixbeg <- n }
            else
            {
              hit <- FALSE
              j <- 0
              while(j < n && hit == FALSE)
              {
                j <- j+1
                if(sca <= x[j])
                {
                  hit <- TRUE
                  Ixbeg <- j
                }
              }
            }
          }
          else { message("beg not set.") }
        }
      }
      if(!is.null(end))
      {
        if(end == Inf) { Ixend <- n }
        else if(end == -Inf) { Ixend <- Ixbeg }
        else
        {
          sca <- private$extract_scalar(end)
          if(!is.null(sca))
          {
            if(sca > x[n]) { Ixend <- n }
            else if(sca < x[Ixbeg]) { Ixend <- Ixbeg }
            else
            {
              hit <- FALSE
              j <- n+1
              while(j > 1 && hit == FALSE)
              {
                j <- j-1
                if(sca >= x[j])
                {
                  hit <- TRUE
                  Ixend <- j
                }
              }
            }
          }
          else { message("end not set.") }
        }
      }
      return(list(Ixbeg=Ixbeg,Ixend=Ixend))
    },
    xedni = function(x,beg,end)
    {
      n <- length(x)
      Ixbeg <- n
      Ixend <- 1
      if(!is.null(beg))
      {
        if(beg == -Inf) { Ixbeg <- n }
        else if(beg == Inf) { Ixbeg <- 1 }
        else
        {
          sca <- private$extract_scalar(beg)
          if(!is.null(sca))
          {
            if(sca < x[n]) { Ixbeg <- n }
            else if(sca > x[1]) { Ixbeg <- 1 }
            else
            {
              hit <- FALSE
              j <- n+1
              while(j > 1 && hit == FALSE)
              {
                j <- j-1
                if(sca <= x[j])
                {
                  hit <- TRUE
                  Ixbeg <- j
                }
              }
            }
          }
          else { message("beg not set.") }
        }
      }
      if(!is.null(end))
      {
        if(end == Inf) { Ixend <- 1 }
        else if(end == -Inf) { Ixend <- Ixbeg }
        else
        {
          sca <- private$extract_scalar(end)
          if(!is.null(sca))
          {
            if(sca > x[1]) { Ixend <- 1 }
            else if(sca < x[Ixbeg]) { Ixend <- Ixbeg }
            else
            {
              hit <- FALSE
              j <- 0
              while(j < n && hit == FALSE)
              {
                j <- j+1
                if(sca >= x[j])
                {
                  hit <- TRUE
                  Ixend <- j
                }
              }
            }
          }
          else { message("end not set.") }
        }
      }
      return(list(Ixbeg=Ixbeg,Ixend=Ixend))
    },
    # private clipboard methods ----
    CopyToClipboard = function(clip)
    {
      if(!is.null(private$OUP)) { OUP$CopyToClipboard(clip) }
      else if(interactive() && clipr_available()) { write_clip(clip,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    }
  )
)
