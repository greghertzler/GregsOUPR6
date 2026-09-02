library(R6)
library(plotly)
library(stringr)
library(clipr)

# roxygen ----
#' @title R6 class implementing Analytical formulas for the Ornstein-Uhlenbeck Process
#'
#' @description
#' A set of functions to calculate probabilities, option prices, visiting
#'  times, first passage times and decision thresholds--everything for
#'  simple exit and entry options in a Real Options Analysis.
#'
#' @details # Formulas and Methods:
#'     z stochastic
#'       Drift
#'       Diffusion
#'     y stochastic
#'       Mean
#'       Variance
#'       Density
#'       Probability
#'       DoubleIntegral
#'     x stochastic
#'       Option
#'       OptionEnvelope
#'       DecisionThreshold
#'       Obligation
#'     t stochastic
#'       PassageTimeModeMedianMean
#'       PassageTimePercentiles
#'       PassageTimeDensity
#'       PassageTimeProbability
#'
#' @details # Plots
#'       PlotDrift
#'       PlotDiffusion
#'       PlotMean
#'       PlotVariance
#'       PlotDensity
#'       PlotProbability
#'       PlotDoubleIntegral
#'       PlotOption
#'       PlotOptionEnvelope
#'       PlotDecisionThreshold
#'       PlotObligation
#'       PlotPassageTimeModeMedianMean
#'       PlotPassageTimePercentiles
#'       PlotPassageTimeDensity
#'       PlotPassageTimeProbability
#'
#' @details # Arguments of functions:
#'       All arguments are optional in all functions.
#'     OUP parameters
#'       rho:    rate parameter 0<=rho<inf
#'       mu:     location parameter -inf<mu<inf
#'       sigma:  scale parameter -inf<sigma<inf
#'     z stochastic
#'       z:      vector of states -inf<z<inf
#'     y stochastic
#'       t:      vector of times s<=t<inf
#'       y:      vector of states -inf<y<inf
#'       s:      initial time -inf<s<inf
#'       x:      initial state -inf<x<inf
#'       psi:    <=0 for integral -inf to y,
#'                >0 for integral y to inf
#'       eps:    proportion remaining after convergence 0<=eps<=1
#'     x stochastic
#'       s:      vector of times -inf<s<t
#'       x:      vector of states -inf<x<inf
#'       t:      terminal time -inf<t<inf
#'       y:      terminal state -inf<y<inf
#'       r:      discount rate -inf<r<inf
#'       phi:    <=0 for exit option,
#'                >0 for entry option
#'     t stochastic
#'       t:      vector of times s<=t<inf
#'       k:      decision threshold -inf<k<int
#'       s:      initial time -inf<s<inf
#'       x:      initial state -inf<x<inf
#'       omega:  degree of irreversibility 0<=omega<=1
#'       Ppct:   passage time probability for a percentile 0.01<=Ppct<=0.99
#'
#' @details # Usage:
#' The Analytical object must first be instantiated before its methods are called.
#'  There are two ways.  The first way instantiates the OUProcess object and
#'  calls a function to get a pointer:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       FD <- OUP$get_FiniteDifference()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'
#' The Analytical object will coordinate arguments to functions with the
#'  FiniteDifference, MaximumLikelihood and MonteCarlo objects.  The second
#'  way instantiates the Analytical object by itself with no coordination:
#'
#'       A <- Analytical$new()
#'
#' Once the object is instantiated, its methods can be called, to calculate and
#'  plot a Decision Threshold, for example:
#'
#'       A$DecisionThreshold()
#'
#' The plot methods can be used to customize the plots, with a title, for example:
#'
#'       A$PlotDecisionThreshold(title="My Decision")
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
#' Mathematical formulas and convergence criteria assume exact arithmetic.
#'  Floating-point arithmetic is another matter. It may overflow, underflow
#'  or have extreme cancellation. In principle, everything can be calculated
#'  by Monte Carlo simulation or Finite Difference methods. These solve
#'  continuous equations at discrete nodes, creating another source of error.
#'  Monte Carlo simulation of passage times can be accurate or biased,
#'  sometimes wildly so. The Finite Difference method can be accurate or give
#'  reasonable looking solutions that diverge from the true answers.
#'
#' Even if great care is taken, it is a good idea to calculate a solution in
#'  more than one way. This is one reason for analytical formulas.
#'
#' Another reason is to speed the calculations. Compared with the Finite
#'  Difference method, an analytical formula can calculate option prices
#'  2 to 3 times quicker. Compared with Monte Carlo simulation, an
#'  analytical formula can calculate 800 to 1,000 times quicker.
#'
#' The bottom line of a real options analysis is the decision threshold
#'  and the time to crossing the decision threshold. The many functions available
#'  are used in the calculations, but they all accumulate to the bottom line:
#'
#'       A <- Analytical$new()
#'       A$DecisionThreshold()
#'       A$PassageTimeModeMedianMean()
#'       A$PassageTimePercentiles()
#'
#' Perhaps the better measure of passage times are the percentiles. The mean
#'  does not exist if rho is zero and the variance does not exist if rho is
#'  small. The mode and median always exist.  Percentiles, which include the
#'  median, always exist.  The formulas don't calculate the variance.
#'
#' The passage time calculations are challenging. For example, the passage
#'  time density is too complicated to explain in this short discussion.
#'  But don't be surprised if it is bi-modal or even negative. Care must be
#'  taken in interpreting passage times.

# class ----
#' @import plotly
#' @import stringr
#' @importFrom clipr clipr_available write_clip
#' @export
Analytical <- R6::R6Class("Analytical",
  portable = FALSE,
  cloneable = FALSE,
  # portable = TRUE,
  # cloneable = TRUE,
  #public members ----
  public = list(
    # constructor ----
    #' @description
    #' Create an Analytical object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new Analytical object
    initialize = function(OUP=NULL)
    {
      # pointer to container object ----
      if(!is.null(OUP) && class(OUP)[[1]] == "OUProcess") { private$OUP <- OUP }
      # arguments ----
      private$oup_params <- list(rho=0.5,mu=15,sigma=15)
      xyzseq <- seq(from=-30,to=30,by=0.6)
      private$z_stoch_args <- list(z=xyzseq)
      private$y_stoch_args <- list(t=seq(from=0,to=10,by=0.1),y=xyzseq,s=0,x=-15,psi=-1,eps=0.05)
      private$x_stoch_args <- list(s=seq(from=10,to=0,by=-0.1),x=xyzseq,t=10,y=0,r=0.05,phi=0,b=0,c=0)
      private$t_stoch_args <- list(t=seq(from=0,to=10,by=0.1),k=0,s=0,x=15,z=seq(from=-30,to=30,by=0.6),omega=1,Ppct=0.75)
      private$plot_args <- list(pmax=0.06,ptmax=0.6)
      private$psiphi <- -1
      # undo ----
      private$undo_args <- list(list(oup_params=private$oup_params,z_stoch_args=private$z_stoch_args,y_stoch_args=private$y_stoch_args,x_stoch_args=private$x_stoch_args,t_stoch_args=private$t_stoch_args,plot_args=private$plot_args))
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
      private$plot_types <- c(rep(0,6))
    },
    # public set methods ----
    #' @description
    #' Set OUP parameters
    #' @param rho   rate parameter 0<=rho<inf
    #' @param mu    location parameter -inf<mu<inf
    #' @param sigma scale parameter -inf<sigma<inf
    #' @param who   object id of caller
    #' @return list(rho,mu,sigma)
    set_oup_params = function(rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"A") }
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
            private$g <- NULL
            private$G <- NULL
            private$Gteps <- NULL
            private$H2 <- NULL
            private$H2teps <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
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
            private$g <- NULL
            private$G <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
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
            private$h2 <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
          }
        }
        else { message("sigma not set.")}
      }
      return(private$oup_params)
    },
    #' @description
    #' Set z as a stochastic state
    #' @param z vector of n states -inf<z<inf
    #' @return list(z)
    set_z_stoch_args = function(z=NULL)
    {
      if(!is.null(z))
      {
        vec <- private$extract_vector(z,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(z,private$z_stoch_args$z))
          {
            private$z_stoch_args$z <- vec
            private$g <- NULL
            private$h2 <- NULL
          }
        }
        else { message("z not set.") }
      }
      return(private$z_stoch_args)
    },
    #' @description
    #' Set y as a stochastic state and its arguments
    #' @param t   vector of m times s<=t<inf
    #' @param y   vector of n states -inf<y<inf
    #' @param s   initial time -inf<s<inf
    #' @param x   initial state -inf<x<inf
    #' @param psi <=0 for integral -inf to y, >0 for integral y to inf
    #' @param eps proportion remaining after convergence 0<=eps<=1
    #' @param who object id of caller
    #' @return list(t,y,s,x,psi,eps)
    set_y_stoch_args = function(t=NULL,y=NULL,s=NULL,x=NULL,psi=NULL,eps=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_y_stoch_args(t,y,x,psi,"A") }
      hit <- FALSE
      if(!is.null(t))
      {
        vec <- private$extract_vector(t,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$y_stoch_args$t))
          {
            private$y_stoch_args$t <- vec
            private$G <- NULL
            private$H2 <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            hit <- TRUE
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
            hit <- TRUE
          }
        }
        else { message("y not set.") }
      }
      if(!is.null(s))
      {
        sca <- private$extract_scalar(s)
        if(!is.null(sca))
        {
          if(sca != private$y_stoch_args$s)
          {
            private$y_stoch_args$s <- sca
            private$G <- NULL
            private$Gteps <- NULL
            private$H2 <- NULL
            private$H2teps <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            hit <- TRUE
          }
        }
        else { message("s not set.") }
      }
      if(!is.null(x))
      {
        sca <- private$extract_scalar(x)
        if(!is.null(sca))
        {
          if(sca != private$y_stoch_args$x)
          {
            private$y_stoch_args$x <- sca
            private$G <- NULL
            private$p <- NULL
            private$Pneg <- NULL
            private$Ppos <- NULL
            private$PPneg <- NULL
            private$PPpos <- NULL
            hit <- TRUE
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
            hit <- TRUE
          }
        }
        else { message("psi not set.") }
      }
      if(!is.null(eps))
      {
        sca <- private$extract_scalar(eps)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            sca <- 0
            message("eps has been set to 0.")
          }
          else if(sca > 1)
          {
            sca <- 1
            message("eps has been set to 1.")
          }
          if(sca != private$y_stoch_args$eps)
          {
            private$y_stoch_args$eps <- sca
            private$Gteps <- NULL
            private$H2teps <- NULL
            hit <- TRUE
          }
        }
        else { message("eps not set.") }
      }
      t1 <- private$y_stoch_args$t[1]
      if(private$y_stoch_args$s > t1)
      {
        private$y_stoch_args$s <- t1
        message(paste(sep="","s has been set to ",t1,"."))
        private$G <- NULL
        private$H2 <- NULL
        private$p <- NULL
        private$Pneg <- NULL
        private$Ppos <- NULL
        private$PPneg <- NULL
        private$PPpos <- NULL
      }
      if(hit == TRUE) { private$psiphi <- private$y_stoch_args$psi}
      return(private$y_stoch_args)
    },
    #' @description
    #' Set x as a stochastic state and its arguments
    #' @param s   vector of m times -inf<s<t
    #' @param x   vector of n states -inf<x<inf
    #' @param t   terminal time -inf<t<inf
    #' @param y   terminal state -inf<y<inf
    #' @param r   discount rate -inf<r<inf
    #' @param phi <=0 for exit option, >0 for entry option
    #' @param b   lump-sum benefit for entry option
    #' @param c   lump-sum cost for exit option
    #' @param who object id of caller
    #' @return list(s,x,t,y,r,phi)
    set_x_stoch_args = function(s=NULL,x=NULL,t=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,y,r,phi,"A") }
      hit <- FALSE
      if(!is.null(s))
      {
        vec <- private$extract_vector(s,-1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$x_stoch_args$s))
          {
            private$x_stoch_args$s <- vec
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
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
              private$OOneg <- NULL
              private$OOpos <- NULL
              private$OOhatneg <- NULL
              private$OOhatpos <- NULL
              private$shatneg <- NULL
              private$shatpos <- NULL
              private$dOOdsconvexneg <- NULL
              private$dOOdsconvexpos <- NULL
              private$dOOdsconcaveneg <- NULL
              private$dOOdsconcavepos <- NULL
              private$dOOdspatchneg <- NULL
              private$dOOdspatchpos <- NULL
              private$BCneg <- NULL
              private$BCpos <- NULL
              hit <- TRUE
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
      if(!is.null(t))
      {
        sca <- private$extract_scalar(t)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$t)
          {
            private$x_stoch_args$t <- sca
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
          }
        }
        else { message("t not set.") }
      }
      if(!is.null(y))
      {
        sca <- private$extract_scalar(y)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$y)
          {
            private$x_stoch_args$y <- sca
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
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
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
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
            hit <- TRUE
          }
        }
        else { message("phi not set.") }
      }
      if(!is.null(b))
      {
        sca <- private$extract_scalar(b)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$b)
          {
            private$x_stoch_args$b <- sca
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
          }
        }
        else { message("b not set.") }
      }
      if(!is.null(c))
      {
        sca <- private$extract_scalar(c)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$c)
          {
            private$x_stoch_args$c <- sca
            private$OOneg <- NULL
            private$OOpos <- NULL
            private$OOhatneg <- NULL
            private$OOhatpos <- NULL
            private$shatneg <- NULL
            private$shatpos <- NULL
            private$dOOdsconvexneg <- NULL
            private$dOOdsconvexpos <- NULL
            private$dOOdsconcaveneg <- NULL
            private$dOOdsconcavepos <- NULL
            private$dOOdspatchneg <- NULL
            private$dOOdspatchpos <- NULL
            private$KOOneg <- NULL
            private$KOOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
          }
        }
        else { message("c not set.") }
      }
      sm <- private$x_stoch_args$s[1]
      if(private$x_stoch_args$t < sm)
      {
        private$x_stoch_args$t <- sm
        message(paste(sep="","t has been set to ",sm,"."))
        private$OOneg <- NULL
        private$OOpos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
      }
      if(hit == TRUE) { private$psiphi <- private$x_stoch_args$phi}
      return(private$x_stoch_args)
    },
    #' @description
    #' Set t stochastic arguments
    #' @param t     vector of m times s<=t<inf
    #' @param k     threshold -inf<k<inf
    #' @param s     initial time -inf<s<inf
    #' @param x     initial state -inf<x<inf
    #' @param z     vector of n alternate initial states -inf<z<inf
    #' @param omega degree of irreversibility 0<=omega<=1
    #' @param Ppct  passage time probability for a percentile  0.01<=Ppct<=0.99
    #' @param who   object id of caller
    #' @return list(t,k,s,x,z,omega)
    set_t_stoch_args = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,Ppct=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_t_stoch_args(t,k,x,omega,Ppct,"A") }
      if(!is.null(t))
      {
        vec <- private$extract_vector(t,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$t_stoch_args$t))
          {
            private$t_stoch_args$t <- vec
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
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
            zz <- private$t_stoch_args$z
            n <- length(zz)
            if(n > 1)
            {
              zzBy <- (zz[n]-zz[1])/(n-1)
              if(sca > 0)
              {
                zzFrom <- sca-as.integer((sca-zz[1])/zzBy)*zzBy
                zzTo <- zzFrom+(n-1)*zzBy
              }
              else
              {
                zzTo <- sca+as.integer((zz[n]-sca)/zzBy)*zzBy
                zzFrom <- zzTo-(n-1)*zzBy
              }
              private$t_stoch_args$z <- seq(from=zzFrom,to=zzTo,by=zzBy)
            }
            else { private$t_stoch_args$z <- sca }
            private$t_stoch_args$k <- sca
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
          }
        }
        else { message("k not set.") }
      }
      if(!is.null(s))
      {
        sca <- private$extract_scalar(s)
        if(!is.null(sca))
        {
          if(sca != private$t_stoch_args$s)
          {
            private$t_stoch_args$s <- sca
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
          }
        }
        else { message("s not set.") }
      }
      if(!is.null(x))
      {
        sca <- private$extract_scalar(x)
        if(!is.null(sca))
        {
          if(sca != private$t_stoch_args$x)
          {
            private$t_stoch_args$x <- sca
            private$tmodemedianmean <- NULL
            private$tpercentile <- NULL
            private$ptx <- NULL
            private$Ptx <- NULL
          }
        }
        else { message("x not set.") }
      }
      if(!is.null(z))
      {
        vec <- private$extract_vector(z,1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$t_stoch_args$z))
          {
            kk <- private$t_stoch_args$k
            n <- length(vec)
            if(n > 1)
            {
              zzBy <- (vec[n]-vec[1])/(n-1)
              if(kk > 0)
              {
                zzFrom <- kk-as.integer((kk-vec[1])/zzBy)*zzBy
                zzTo <- zzFrom+(n-1)*zzBy
              }
              else
              {
                zzTo <- kk+as.integer((vec[n]-kk)/zzBy)*zzBy
                zzFrom <- zzTo-(n-1)*zzBy
              }
              private$t_stoch_args$z <- seq(from=zzFrom,to=zzTo,by=zzBy)
            }
            else { private$t_stoch_args$z <- kk }
            private$tmodesmediansmeans <- NULL
            private$tpercentiles <- NULL
            private$pt <- NULL
            private$Pt <- NULL
          }
        }
        else { message("z not set.") }
      }
      if(!is.null(omega))
      {
        sca <- private$extract_scalar(omega)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            sca <- 0
            message("omega has been set to 0.")
          }
          else if(sca > 1)
          {
            sca <- 1
            message("omega has been set to 1.")
          }
          if(sca != private$t_stoch_args$omega)
          {
            private$t_stoch_args$omega <- sca
            private$tmodemedianmean <- NULL
            private$tmodesmediansmeans <- NULL
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
            private$ptx <- NULL
            private$pt <- NULL
            private$Ptx <- NULL
            private$Pt <- NULL
          }
        }
        else { message("omega not set.") }
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
            private$tpercentile <- NULL
            private$tpercentiles <- NULL
          }
        }
        else { message("Ppct not set.") }
      }
      t1 <- private$t_stoch_args$t[1]
      if(private$t_stoch_args$s > t1)
      {
        private$t_stoch_args$s <- t1
        message(paste(sep="","s has been set to ",t1,"."))
        private$tmodemedianmean <- NULL
        private$tmodesmediansmeans <- NULL
        private$tpercentile <- NULL
        private$tpercentiles <- NULL
        private$ptx <- NULL
        private$pt <- NULL
        private$Ptx <- NULL
        private$Pt <- NULL
      }
      return(private$t_stoch_args)
    },
    #' @description
    #' Set plot arguments
    #' @param pmax  maximum transition density
    #' @param ptmax maximum visiting time and first passage time densities
    #' @param who   object id of caller
    #' @return list(pmax,ptmax)
    set_plot_args = function(pmax=NULL,ptmax=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_args(pmax,ptmax,"A") }
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
    #' @param who        object id of caller
    #' @return list(type,font,file,theme,3D)
    set_plot_info = function(fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,"A") }
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
        if(group == 2) # Diffusion
        {
          least <- -1
          most <- 0
        }
        else if(group == 3) # Mean, Variance
        {
          least <- -1
          most <- 2
        }
        else if(group == 4) # Density, Probability, Double Integral, Option, Option Envelope, Obligation
        {
          least <- 0
          most <- 1
        }
        else if(group == 5) # Passage Time Mode Median Mean, Passage Time Percentiles
        {
          least <- -3
          most <- 2
        }
        else if(group == 6) # Passage Time Density, Passage Time Probability
        {
          least <- 0
          most <- 1
        }
        else # Drift, Decision Threshold
        {
          group = 1
          least <- 0
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
    #' @param who    object id of caller
    #' @return list(plotit,copyit)
    set_flags = function(plotit=NULL,copyit=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_flags(plotit,copyit,"A") }
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
    #' @return list(oup_params,z_stoch_args,y_stoch_args,x_stoch_args,t_stoch_args,plot_info)
    get_all = function()
    {
      all <- list(oup_params = private$oup_params,
        z_stoch_args = private$z_stoch_args,
        y_stoch_args = private$y_stoch_args,
        x_stoch_args = private$x_stoch_args,
        t_stoch_args = private$t_stoch_args,
        plot_args = private$plot_args,
        plot_info = private$plot_info)
      return(all)
    },
    #' @description
    #' Get OUP parameters
    #' @return list(rho,mu,sigma)
    get_oup_params = function() { return(private$oup_params) },
    #' @description
    #' Get z as a stochastic state
    #' @return list(z)
    get_z_stoch_args = function() { return(private$z_stoch_args) },
    #' @description
    #' Get y as a stochastic state and its arguments
    #' @return list(t,y,s,x,psi,eps)
    get_y_stoch_args = function() { return(private$y_stoch_args) },
    #' @description
    #' Get x as a stochastic state and its arguments
    #' @return list(s,x,t,y,r,phi)
    get_x_stoch_args = function() { return(private$x_stoch_args) },
    #' @description
    #' Get t stochastic arguments
    #' @return list(t,k,s,x,z,omega,Ppct)
    get_t_stoch_args = function() { return(private$t_stoch_args) },
    #' @description
    #' Get plot arguments
    #' @return list(pmax,ptmax)
    get_plot_args = function() { return(private$plot_args) },
    #' @description
    #' Get information for plotting
    #' @return list(type,font,file,theme,3D,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
    #' @description
    #' Get current types for plot routines
    #' @return (list(types,descriptioon))
    get_plot_types = function()
    {
      text <- rbind(c("Analytical groups, types and plots (default type is 0):"),c("  group  types  plots"),c("    1        0  Drift DecisionThreshold"),c("    2     -1,0  Diffusion"),c("    3    -1,,2  Mean Variance"),c("    4      0,1  Density Probability DoubleIntegral Option OptionEnvelope Obligation"),c("    5    -3,,2  PassageTimeModeMedianMean PassageTimePercentiles"),c("    6      0,1  PassageTimeDensity PassageTimeProbability"))
      return(list(types=private$plot_types,description=text))
    },
    #' @description
    #' Get flags for plotting and copying
    #' @return list(plotit,copyit)
    get_flags = function() { return(private$flags) },
    # public axis and sync methods ----
    #' @description
    #' Scale axes for z stochastic arguments
    #' @return NULL
    axes_z_stoch = function()
    {
      # get
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      # state
      if(rho > 0) { z <- 2*abs(sigma)/(2*rho)^0.5 }
      else { z <- 2*abs(sigma) }
      if(z < 1) { z <- 1}
      zscale <- 1
      while(z > zscale) { zscale <- 10*zscale }
      z <- round(z/zscale,1)*zscale
      zby <- z/50
      zscale <- 1
      while(abs(mu) > zscale) { zscale <- 10*zscale }
      z <- round(mu/zscale,2)*zscale
      zfrom <- z-50*zby
      zto <- zfrom+100*zby
      zseq <- seq(from=zfrom,to=zto,by=zby)

      self$set_z_stoch_args(zseq)

      return(NULL)
    },
    #' @description
    #' Scale axes for y stochastic arguments
    #' @return NULL
    axes_y_stoch = function()
    {
      # get
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      eps <- private$y_stoch_args[[6]]
      # state
      if(rho > 0) { y <- abs(x-mu)+4*sigma/(2*rho)^0.5}
      else { y <- abs(x-mu)+4*sigma }
      if(y < 1) { y <- 1}
      yscale <- 1
      while(y > yscale) { yscale <- 10*yscale }
      y <- round(y/yscale,1)*yscale
      yby <- y/100
      yscale <- 1
      while(y > yscale) { yscale <- 10*yscale }
      y <- round(y/yscale,2)*yscale
      yfrom <- 0.5*(x+mu)-50*yby
      yto <- yfrom+100*yby
      yseq <- seq(from=yfrom,to=yto,by=yby)
      # time
      t <- 100
      if(eps > 0 && rho > 0) { t <- -1.6*log(eps)/rho }
      if(t > 100) { t <- 100 }
      else if(t < 1) { t <- 1}
      tscale <- 1
      while(t > tscale) { tscale <- 10*tscale }
      t <- round(t/tscale,1)*tscale
      tfrom <- s
      tto <- s+t
      tby <- t/100
      tseq <- seq(from=tfrom,to=tto,by=tby)
      self$set_y_stoch_args(tseq,yseq,NULL,NULL,NULL,NULL)

      # density
      if(rho < 0.0000000001) { variance <- sigma^2*t }
      else { variance <- sigma^2/(2*rho)*(1-exp(-2*rho*t)) }
      if(variance > 0.00994718394324347) { pmax <- 2.5/(2*3.14159265358979*variance)^0.5 }
      else { pmax = 10 }
      pscale <- 0.01
      while(pmax > pscale) { pscale <- 10*pscale }
      pmax <- round(pmax/pscale,2)*pscale
      self$set_plot_args(pmax,NULL)

      return(NULL)
    },
    #' @description
    #' Scale axes for x stochastic arguments
    #' @return NULL
    axes_x_stoch = function()
    {
      # get
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      y <- private$x_stoch_args[[4]]
      # state
      k <- self$DecisionThreshold(who="A")[[1]]
      if(is.finite(k))
      {
        x <- 2*abs(k-y)
        if(x < 1) { x <- 1}
        xscale <- 1
        while(x > xscale) { xscale <- 10*xscale }
        x <- round(x/xscale,1)*xscale
        xby <- x/50
        xscale <- 1
        while(abs(k) > xscale) { xscale <- 10*xscale }
        x <- round(k/xscale,2)*xscale
        if(k < y) { xfrom <- x-30*xby }
        else { xfrom <- x-70*xby }
        xto <- xfrom+100*xby
      }
      else
      {
        x <- 2*abs(sigma)
        if(x < 1) { x <- 1}
        xscale <- 1
        while(x > xscale) { xscale <- 10*xscale }
        x <- round(x/xscale,1)*xscale
        xby <- x/50
        xscale <- 1
        while(abs(y) > xscale) { xscale <- 10*xscale }
        x <- round(y/xscale,2)*xscale
        xfrom <- x-50*xby
        xto <- xfrom+100*xby
      }
      xseq <- seq(from=xfrom,to=xto,by=xby)
      # time
      t <- 100
      if(rho > 0) { t <- -1.6*log(0.05)/rho }
      if(t > 100) { t <- 100 }
      else if(t < 1) { t <- 1}
      tscale <- 1
      while(t > tscale) { tscale <- 10*tscale }
      t <- round(t/tscale,1)*tscale
      sfrom <- s[1]
      sto <- sfrom-t
      sby <- -t/100
      sseq <-seq(from=sfrom,to=sto,by=sby)

      self$set_x_stoch_args(sseq,xseq,NULL,NULL,NULL,NULL,NULL,NULL)

      return(NULL)
    },
    #' @description
    #' Scale axes for t stochastic arguments
    #' @return NULL
    axes_t_stoch = function()
    {
      # get
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      omega <- private$t_stoch_args[[6]]
      # state
      if(is.finite(k))
      {
        z <- 2*abs(k-x)
        if(z < 1) { z <- 1}
        zscale <- 1
        while(z > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,1)*zscale
        zby <- z/50
        if(k < x) { zfrom <- k-30*zby }
        else { zfrom <- k-70*zby }
        zto <- zfrom+100*zby
      }
      else
      {
        z <- 2*abs(sigma)
        if(z < 1) { z <- 1}
        zscale <- 1
        while(z > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,1)*zscale
        zby <- z/50
        zscale <- 1
        while(abs(z) > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,2)*zscale
        zfrom <- z-50*zby
        zto <- zfrom+100*zby
      }
      zseq <- seq(from=zfrom,to=zto,by=zby)
      # time
      tpercentiles <- self$PassageTimePercentiles(who="A")[[2]]
      tuppers <- tpercentiles[[3]]
      n <- length(tuppers)
      if(is.infinite(tuppers[1])) { tuppers[1] <- s+1 }
      if(is.infinite(tuppers[n])) { tuppers[n] <- s+1 }
      if(tuppers[n] > tuppers[1]) { t <- (tuppers[n]-s) }
      else { t <- (tuppers[1]-s) }
      if(t < 1) { t <- 1}
      tscale <- 1
      while(t > tscale) { tscale <- 10*tscale }
      t <- round(t/tscale,2)*tscale
      tfrom <- s
      tto <- t+s
      tby <- t/100
      tseq <- seq(from=tfrom,to=tto,by=tby)

      self$set_t_stoch_args(tseq,NULL,NULL,NULL,zseq,NULL,NULL)

      # density
      tmodemedianmean <- self$PassageTimeModeMedianMean(who="A")[[1]]
      tmode <- tmodemedianmean[[1]]
      if(is.finite(tmode))
      {
        ptmax <- 1.2*private$OUPPassageTimeDensity(s,x,tmode,k,omega,rho,mu,sigma,tby)
        pscale <- 0.01
        while(ptmax > pscale) { pscale <- 10*pscale }
        ptmax <- round(ptmax/pscale,2)*pscale
        private$plot_info$plottype$ptmax <- ptmax
        self$set_plot_args(NULL,ptmax)
      }
      return(NULL)
    },
    #' @description
    #' Synchronize states
    #' @return NULL
    sync_zyxt_stoch = function()
    {
      # get
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      if(phi <= 0)
      {
        decision <- private$KOOneg
        if(is.null(decision))
        {
          k <- mu
          if(rho > 0) { x <- k+sigma*(2/rho)^0.5 }
          else { x <= k }
        }
        else
        {
          x <- self$DecisionThreshold(phi=1,who="A")[[1]]
          k <- self$DecisionThreshold(phi=-1,who="A")[[1]]
        }
      }
      else
      {
        decision <- private$KOOpos
        if(is.null(decision))
        {
          k <- mu
          if(rho > 0) { x <- k-sigma*(2/rho)^0.5 }
          else { x <- k }
        }
        else
        {
          x <- self$DecisionThreshold(phi=-1,who="A")[[1]]
          k <- self$DecisionThreshold(phi=1,who="A")[[1]]
        }
      }
      # t state
      if(is.finite(k))
      {
        z <- 2*abs(k-x)
        if(z < 1) { z <- 1}
        zscale <- 1
        while(z > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,1)*zscale
        zby <- z/50
        if(k < x) { zfrom <- k-30*zby }
        else { zfrom <- k-70*zby }
        zto <- zfrom+100*zby
      }
      else
      {
        z <- 2*abs(sigma)
        if(z < 1) { z <- 1}
        zscale <- 1
        while(z > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,1)*zscale
        zby <- z/50
        zscale <- 1
        while(abs(z) > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,2)*zscale
        zfrom <- z-50*zby
        zto <- zfrom+100*zby
      }
      zseq <- seq(from=zfrom,to=zto,by=10*zby)
      xyzseq <- seq(from=zfrom,to=zto,by=zby)

      self$set_t_stoch_args(NULL,k,NULL,x,zseq,NULL,NULL)

      # z state
      self$set_z_stoch_args(xyzseq)

      # y state
      self$set_y_stoch_args(NULL,xyzseq,NULL,NULL,private$psiphi,NULL)

      # x state
      self$set_x_stoch_args(NULL,xyzseq,NULL,NULL,NULL,private$psiphi,NULL,NULL)

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
        z_stoch_args=private$z_stoch_args,
        y_stoch_args=private$y_stoch_args,
        x_stoch_args=private$x_stoch_args,
        t_stoch_args=private$t_stoch_args,
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
        if(private$lists_equal(last_undo_args[[2]],private$z_stoch_args))
        {
          if(private$lists_equal(last_undo_args[[3]],private$y_stoch_args))
          {
            if(private$lists_equal(last_undo_args[[4]],private$x_stoch_args))
            {
              if(private$lists_equal(last_undo_args[[5]],private$t_stoch_args))
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
            z_stoch_args=private$z_stoch_args,
            y_stoch_args=private$y_stoch_args,
            x_stoch_args=private$x_stoch_args,
            t_stoch_args=private$t_stoch_args,
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
      z_stoch <- these_undo[[2]]
      y_stoch <- these_undo[[3]]
      x_stoch <- these_undo[[4]]
      t_stoch <- these_undo[[5]]
      plot <- these_undo[[6]]
      self$set_oup_params(oup[[1]],oup[[2]],oup[[3]])
      self$set_z_stoch_args(z_stoch[[1]])
      self$set_y_stoch_args(y_stoch[[1]],y_stoch[[2]],y_stoch[[3]],y_stoch[[4]],y_stoch[[5]],y_stoch[[6]])
      self$set_x_stoch_args(x_stoch[[1]],x_stoch[[2]],x_stoch[[3]],x_stoch[[4]],x_stoch[[5]],x_stoch[[6]],x_stoch[[7]],x_stoch[[8]])
      self$set_t_stoch_args(t_stoch[[1]],t_stoch[[2]],t_stoch[[3]],t_stoch[[4]],t_stoch[[5]],t_stoch[[6]],t_stoch[[7]])
      self$set_plot_args(plot[[1]],plot[[2]])
      private$undoIx <- undoIx

      return(c(undoIx,n))
    },
    # public calculate methods ----
    #' @description
    #' Calculate, plot and return drifts
    #' @param z       vector of n states
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param who     object id of caller
    #' @return list(g(1xn))
    Drift = function(z=NULL,rho=NULL,mu=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_z_stoch_args(z)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      z <- private$z_stoch_args[[1]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      drifts <- private$g
      if(is.null(drifts))
      {
        drifts <- RcppOUPADrift(z,rho,mu)
        private$g <- drifts
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDrift()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          clip <- rbind(c("Analytical",rep("",n)),c("Drift",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("z",z),c("g",drifts))
          private$CopyToClipboard(clip)
        }
      }
      return(list(g=drifts))
    },
    #' @description
    #' Calculate, plot and return diffusions
    #' @param z       vector of n states
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(h2(1xn))
    Diffusion = function(z=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(NULL,NULL,sigma)
      self$set_z_stoch_args(z)
      sigma <- private$oup_params[[3]]
      z <- private$z_stoch_args[[1]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      diffusions <- private$h2
      if(is.null(diffusions))
      {
        diffusions <- RcppOUPADiffusion(z,sigma)
        private$h2 <- diffusions
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDiffusion()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          clip <- rbind(c("Analytical",rep("",n)),c("Diffusion",rep("",n)),c("sigma",sigma,rep("",n-1)),c("z",z),c("h\u00B2",diffusions))
          private$CopyToClipboard(clip)
        }
      }
      return(list(h2=diffusions))
    },
    #' @description
    #' Calculate and plot means and time for means to converge
    #' @param t       vector of m times s<=t<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param eps     proportion remaining after convergence 0<=eps<=1
    #' @param who     object id of caller
    #' @return list(G(mx1),Gteps)
    Mean = function(t=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,eps=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_y_stoch_args(t,NULL,s,x,NULL,eps)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      eps <- private$y_stoch_args[[6]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      means <- private$G
      timeeps <- private$Gteps
      if(is.null(means) || is.null(timeeps))
      {
        Gt <- RcppOUPAMean(t,s,x,rho,mu,eps)
        m <- length(t)
        means <- Gt[1:m]
        timeeps <- Gt[m+1]
        private$G <- means
        private$Gteps <- timeeps
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotMean()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Analytical",""),c("Mean",""),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("eps",eps),c("Gteps",timeeps),c("t","G"),cbind(t,means))
          private$CopyToClipboard(clip)
        }
      }
      return(list(G=means,Gteps=timeeps))
    },
    #' @description
    #' Calculate and plot variances and time for variances to converge
    #' @param t       vector of m times s<=t<inf
    #' @param s       initial time -inf<s<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param eps     proportion remaining after convergence 0<=eps<=1
    #' @param who     object id of caller
    #' @return list(H2(mx1),H2teps)
    Variance = function(t=NULL,s=NULL,rho=NULL,sigma=NULL,eps=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,NULL,sigma)
      self$set_y_stoch_args(t,NULL,s,NULL,NULL,eps)
      rho <- private$oup_params[[1]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      eps <- private$y_stoch_args[[6]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      variances <- private$H2
      timeeps <- private$H2teps
      if(is.null(variances) || is.null(timeeps))
      {
        H2t <- RcppOUPAVariance(t,s,rho,sigma,eps)
        m <- length(t)
        variances <- H2t[1:m]
        timeeps <- H2t[m+1]
        private$H2 <- variances
        private$H2teps <- timeeps
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotVariance()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Analytical",""),c("Variance",""),c("s",s),c("rho",rho),c("sigma",sigma),c("eps",eps),c("H2teps",timeeps),c("t","H\u00B2"),cbind(t,variances))
          private$CopyToClipboard(clip)
        }
      }
      return(list(H2=variances,H2teps=timeeps))
    },
    #' @description
    #' Calculate and plot transition densities
    #' @param t       vector of m times s<=t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(p(mxn))
    Density = function(t=NULL,y=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,s,x,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      densities <- private$p
      if(is.null(densities))
      {
        densities <- RcppOUPADensity(t,y,s,x,rho,mu,sigma)
        private$p <- densities
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDensity()) }
        else if(copyit == TRUE)
        {
          n <- length(y)
          clip <- rbind(c("Analytical",rep("",n)),c("Transition Densities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("p(t,y)",y),cbind(t,densities))
          private$CopyToClipboard(clip)
        }
      }
      return(list(p=densities))
    },
    #' @description
    #' Calculate and plot transition probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param y       vector of nstates -inf<y<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(P(mxn))
    Probability = function(t=NULL,y=NULL,s=NULL,x=NULL,psi=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,s,x,psi,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(psi <= 0)
      {
        probabilities <- private$Pneg
        if(is.null(probabilities))
        {
          probabilities <- RcppOUPAProbability(t,y,s,x,rho,mu,sigma,psi)
          private$Pneg <- probabilities
        }
      }
      else
      {
        probabilities <- private$Ppos
        if(is.null(probabilities))
        {
          probabilities <- RcppOUPAProbability(t,y,s,x,rho,mu,sigma,psi)
          private$Ppos <- probabilities
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotProbability()) }
        else if(copyit == TRUE)
        {
          n <- length(y)
          clip <- rbind(c("Analytical",rep("",n)),c("Transition Probabilities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("P(t,y)",y),cbind(t,probabilities))
          private$CopyToClipboard(clip)
        }
      }
      return(list(P=probabilities))
    },
    #' @description
    #' Calculate and plot double integrals of transition densities
    #' @param t       vector of m times s<=t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(PP(mxn))
    DoubleIntegral = function(t=NULL,y=NULL,s=NULL,x=NULL,psi=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_y_stoch_args(t,y,s,x,psi,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(psi <= 0 )
      {
        doubleintegrals <- private$PPneg
        if(is.null(doubleintegrals))
        {
          doubleintegrals <- RcppOUPADoubleIntegral(t,y,s,x,rho,mu,sigma,psi)
          private$PPneg <- doubleintegrals
        }
      }
      else
      {
        doubleintegrals <- private$PPpos
        if(is.null(doubleintegrals))
        {
          doubleintegrals <- RcppOUPADoubleIntegral(t,y,s,x,rho,mu,sigma,psi)
          private$PPpos <- doubleintegrals
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDoubleIntegral()) }
        else if(copyit == TRUE)
        {
          n <- length(y)
          clip <- rbind(c("Analytical",rep("",n)),c("Double Integrals",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("\u2119(t,y)",y),cbind(t,doubleintegrals))
          private$CopyToClipboard(clip)
        }
      }
      return(list(PP=doubleintegrals))
    },
    #' @description
    #' Calculate and plot option prices
    #' @param s       vector of m times -inf<s<t
    #' @param x       vector of n states -inf<x<inf
    #' @param t       terminal time -inf<t<inf
    #' @param y       terminal state -inf<y<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(OO(mxn))
    Option = function(s=NULL,x=NULL,t=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,x,t,y,r,phi,b,c)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(phi <= 0)
      {
        options <- private$OOneg
        if(is.null(options))
        {
          options <- RcppOUPAOption(s,x,t,y,rho,mu,sigma,r,phi,b,c)
          private$OOneg <- options
        }
      }
      else
      {
        options <- private$OOpos
        if(is.null(options))
        {
          options <- RcppOUPAOption(s,x,t,y,rho,mu,sigma,r,phi,b,c)
          private$OOpos <- options
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotOption()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Options",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("c",c,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options)) }
          else { clip <- rbind(c("Analytical",rep("",n)),c("Options",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options)) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(OO=options))
    },
    #' @description
    #' Calculate and plot the envelope of option prices
    #' @param x       vector of n states -inf<x<inf
    #' @param y       terminal state -inf<y<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(OOhat(1xn),shat(1xn))
    OptionEnvelope = function(x=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / gets ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(NULL,x,NULL,y,r,phi,b,c)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(phi <= 0)
      {
        OOhat <- private$OOhatneg
        shat <- private$shatneg
        if(is.null(OOhat) || is.null(shat))
        {
          OOs <- RcppOUPAOptionEnvelope(s,x,t,y,rho,mu,sigma,r,phi,b,c)
          OOhat <- OOs[1,]
          shat <- OOs[2,]
          private$OOhatneg <- OOhat
          private$shatneg <- shat
        }
      }
      else
      {
        OOhat <- private$OOhatpos
        shat <- private$shatpos
        if(is.null(OOhat)|| is.null(shat))
        {
          OOs <- RcppOUPAOptionEnvelope(s,x,t,y,rho,mu,sigma,r,phi,b,c)
          OOhat <- OOs[1,]
          shat <- OOs[2,]
          private$OOhatpos <- OOhat
          private$shatpos <- shat
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotOptionEnvelope()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Option Envelope",rep("",n)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("c",c,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat)) }
          else { clip <- rbind(c("Analytical",rep("",n)),c("Option Envelope",rep("",n)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat)) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(OOhat=OOhat,shat=shat))
    },
    #' @description
    #' Calculate and plot the decision threshold
    #' @param x       vector of n states -inf<x<inf
    #' @param y       terminal state -inf<y<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(k,OO)
    DecisionThreshold = function(x=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(NULL,x,NULL,y,r,phi,b,c)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(phi <= 0)
      {
        decision <- private$KOOneg
        if(is.null(decision))
        {
          decision <- RcppOUPADecisionThreshold(y,rho,mu,sigma,r,phi,b,c)
          private$KOOneg <- decision
        }
      }
      else
      {
        decision <- private$KOOpos
        if(is.null(decision))
        {
          decision <- RcppOUPADecisionThreshold(y,rho,mu,sigma,r,phi,b,c)
          private$KOOpos <- decision
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDecisionThreshold()) }
        else if(copyit == TRUE)
        {
          if(phi <= 0) { clip <- rbind(c("Analytical",""),c("Decision Threshold",""),c("y",y),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("c",c),c("k",decision[1]),c("\u00D4",decision[2])) }
          else { clip <- rbind(c("Analytical",""),c("Decision Threshold",""),c("y",y),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("b",b),c("k",decision[1]),c("\u00D4",decision[2])) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(k=decision[1],OO=decision[2]))
    },
    #' @description
    #' Calculate and plot obligations and prohibitions, ie. benefit/cost analysis
    #' @param s       vector of m times -inf<s<t
    #' @param x       vector of n states -inf<x<inf
    #' @param t       terminal time -inf<t<inf
    #' @param y       terminal state -inf<y<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for obligation, >0 for prohibition
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param who     object id of caller
    #' @return list(BC(mxn))
    Obligation = function(s=NULL,x=NULL,t=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,rho=NULL,mu=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_x_stoch_args(s,x,t,y,r,phi,b,c)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(phi <= 0)
      {
        obligations <- private$BCneg
        if(is.null(obligations))
        {
          obligations <- RcppOUPAObligation(s,x,t,y,rho,mu,r,phi,b,c)
          private$BCneg <- obligations
        }
      }
      else
      {
        obligations <- private$BCpos
        if(is.null(obligations))
        {
          obligations <- RcppOUPAObligation(s,x,t,y,rho,mu,r,phi,b,c)
          private$BCpos <- obligations
        }
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotObligation()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Obligation",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("c",c,rep("",n-1)),c("\uD835\uDD39(s,x)",x),cbind(s,obligations)) }
          else { clip <- rbind(c("Analytical",rep("",n)),c("Obligation",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("c",c,rep("",n-1)),c("\u2102(s,x)",x),cbind(s,obligations)) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(BC=obligations))
    },
    #' @description
    #' Calculate and plot passage time modes, medians and means
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<z<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(tmodemedianmean,tmodesmediansmeans(n),)
    PassageTimeModeMedianMean = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(NULL,k,s,x,z,omega,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      tmodemedianmean <- private$tmodemedianmean
      tmodesmediansmeans <- private$tmodesmediansmeans
      if(is.null(tmodemedianmean) || is.null(tmodesmediansmeans))
      {
        n <- length(z)
        tmmmtmmmx <- RcppOUPAPassageTimeModeMedianMean(k,s,x,omega,rho,mu,sigma,z)
        tmode <- tmmmtmmmx[1,n+1,drop=FALSE]
        tmedian <- tmmmtmmmx[2,n+1,drop=FALSE]
        tmean <- tmmmtmmmx[3,n+1,drop=FALSE]
        tmodes <- tmmmtmmmx[1,1:n,drop=FALSE]
        tmedians <- tmmmtmmmx[2,1:n,drop=FALSE]
        tmeans <- tmmmtmmmx[3,1:n,drop=FALSE]
        tmodemedianmean <- list(tmode=tmode,tmedian=tmedian,tmean=tmean)
        private$tmodemedianmean <- tmodemedianmean
        tmodesmediansmeans <- list(tmodes=tmodes,tmedians=tmedians,tmeans=tmeans)
        private$tmodesmediansmeans <- tmodesmediansmeans
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotPassageTimeModeMedianMean()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Mode Median Mean",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("tmode(x)"),tmodemedianmean[[1]],paste0("tmode(z)"),tmodesmediansmeans[[1]]),c(paste0("tmedian(x)"),tmodemedianmean[[2]],paste0("tmedian(z)"),tmodesmediansmeans[[2]]),c(paste0("tmean(x)"),tmodemedianmean[[3]],paste0("tmean(z)"),tmodesmediansmeans[[3]]))
          private$CopyToClipboard(clip)
        }
      }
      return(list(tmodemedianmean=tmodemedianmean,tmodesmediansmeans=tmodesmediansmeans))
    },
    #' @description
    #' Calculate and plot passage time percentiles
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<z<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param Ppct    passage time probability for a percentile 0.01<=Ppct<=0.99
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(tpercentile(3),tpercentiles(3xn))
    PassageTimePercentiles = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,Ppct=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(NULL,k,s,x,z,omega,Ppct)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      Ppct <- private$t_stoch_args[[7]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      tpercentile <- private$tpercentile
      tpercentiles <- private$tpercentiles
      if(is.null(tpercentile) || is.null(tpercentiles))
      {
        n <- length(z)
        tpcttpctx <- RcppOUPAPassageTimePercentiles(k,s,x,omega,Ppct,rho,mu,sigma,z)
        tlower <- tpcttpctx[1,n+1,drop=FALSE]
        tmedian <- tpcttpctx[2,n+1,drop=FALSE]
        tupper <- tpcttpctx[3,n+1,drop=FALSE]
        tlowers <- tpcttpctx[1,1:n,drop=FALSE]
        tmedians <- tpcttpctx[2,1:n,drop=FALSE]
        tuppers <- tpcttpctx[3,1:n,drop=FALSE]
        tpercentile <- list(tlower=tlower,tmedian=tmedian,tupper=tupper)
        private$tpercentile <- tpercentile
        tpercentiles <- list(tlowers=tlowers,tmedians=tmedians,tuppers=tuppers)
        private$tpercentiles <- tpercentiles
      }
      # # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotPassageTimePercentiles()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          if(Ppct > 0.5) { clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Percentiles",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("Ppct",Ppct,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("t",1-Ppct,"(x)"),tpercentile[[1]],paste0("t",1-Ppct,"(z)"),tpercentiles[[1]]),c(paste0("t0.5(x)"),tpercentile[[2]],paste0("t0.5(z)"),tpercentiles[[2]]),c(paste0("t",Ppct,"(x)"),tpercentile[[3]],paste0("t",Ppct,"(z)"),tpercentiles[[3]])) }
          else { clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Percentiles",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("Ppct",Ppct,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("t",Ppct,"(x)"),tpercentile[[1]],paste0("t",Ppct,"(z)"),tpercentiles[[1]]),c(paste0("t0.5(x)"),tpercentile[[2]],paste0("t0.5(z)"),tpercentiles[[2]]),c(paste0("t",1-Ppct,"(x)"),tpercentile[[3]],paste0("t",1-Ppct,"(z)"),tpercentiles[[3]])) }
          private$CopyToClipboard(clip)
        }
      }
      return(list(tpercentile=tpercentile,tpercentiles=tpercentiles))
    },
    #' @description
    #' Calculate and plot passage time densities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<z<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(ptx(mx1),pt(mxn))
    PassageTimeDensity = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,s,x,z,omega,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) || is.null(pt))
      {
        n <- length(z)
        ptptx <- RcppOUPAPassageTimeDensity(t,k,s,x,omega,rho,mu,sigma,z)
        pt <- ptptx[,1:n,drop=FALSE]
        ptx <- ptptx[,n+1,drop=FALSE]
        private$ptx <- ptx
        private$pt <- pt
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotPassageTimeDensity()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Densities",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("t","ptx","pt(t,z)",z),cbind(t,ptx,t,pt))
          private$CopyToClipboard(clip)
        }
      }
      return(list(ptx=ptx,pt=pt))
    },
    #' @description
    #' Calculate and plot passage time probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<z<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(Ptx(1xn),Pt(mxn))
    PassageTimeProbability = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_t_stoch_args(t,k,s,x,z,omega,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) || is.null(Pt))
      {
        n <- length(z)
        PtPtx <- RcppOUPAPassageTimeProbability(t,k,s,x,omega,rho,mu,sigma,z)
        Pt <- PtPtx[,1:n,drop=FALSE]
        Ptx <- PtPtx[,n+1,drop=FALSE]
        private$Ptx <- Ptx
        private$Pt <- Pt
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotPassageTimeProbability()) }
        else if(copyit == TRUE)
        {
          n <- length(z)
          clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Probabilities",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("k",k,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("t","Ptx","Pt(t,z)",z),cbind(t,Ptx,t,Pt))
          private$CopyToClipboard(clip)
        }
      }
      return(list(Ptx=Ptx,Pt=Pt))
    },
    # public plot methods ----
    #' @description
    #' Plot drifts
    #' @param type  = 0
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotDrift = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      z <- private$z_stoch_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      drifts <- private$g #protect against recursive call
      if(is.null(drifts)) { drifts <- self$Drift(who="A")[[1]] }
      n <- length(z)
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        drifts <- drifts[Ixbeg:Ixend]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n)),c("Drift",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("z",z),c("g",drifts))
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_A_Drift2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),")",esml)
        if(is.null(title)) { title <- "Drift" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>z</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>z</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "<i>g</i>(<i>z</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside")
      driftline <- list(color=cyn$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Drift2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=z,y=drifts,mode="lines",line=driftline) %>%
        config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot diffusions
    #' @param type  = -1, 0, or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotDiffusion = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,2)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      z <- private$z_stoch_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      drifts <- private$g #no plot or copy
      if(is.null(drifts)) { drifts <- self$Drift(who="A")[[1]] }
      diffusions <- private$h2 #protect against recursive call
      if(is.null(diffusions)) { diffusions <- self$Diffusion(who="A")[[1]] }
      sqrt <- diffusions^0.5
      driftsplus <- drifts+sqrt
      driftsminus <- drifts-sqrt
      n <- length(z)
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        drifts <- drifts[Ixbeg:Ixend]
        diffusions <- diffusions[Ixbeg:Ixend]
        driftsplus <- driftsplus[Ixbeg:Ixend]
        driftsminus <- driftsminus[Ixbeg:Ixend]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {

      }
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n)),c("Diffusion",rep("",n)),c("sigma",sigma,rep("",n-1)),c("z",z),c("h\u00B2",diffusions))
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_A_Diffusion2Dg and OUP_A_Diffusion2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Diffusion" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>z</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>z</i><br>" }
      }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      #OUP_A_Diffusion2Dg
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>g</i>(<i>z</i>)&plusmn;<i>h</i>" }
        lookleft <- list(text=yaxis)
        horz=list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert=list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        diffusionline <- list(color=mgn$d,width=4)
        driftline <- list(color=cyn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Diffusion2Dg")
        legendpos <- list(orientation="h",x=1.05,y=1.0,xanchor="right")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=z,y=drifts,name="<i>g</i>(<i>z</i>)",mode="lines",line=driftline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=z,y=driftsplus,name="<i>g</i>(<i>z</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,hoverinfo="x+y",legendgroup="g+h") %>%
          add_trace(.,type="scatter",x=z,y=driftsminus,name="<i>g</i>(<i>z</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,hoverinfo="x+y",legendgroup="g+h",showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_Diffusion2D
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>h</i><sup>2</sup>" }
        lookleft <- list(text=yaxis)
        horz=list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert=list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        diffusionline <- list(color=mgn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Diffusion2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=z,y=diffusions,mode="lines",line=diffusionline) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot means
    #' @param type  = -1,...,2 or 'n','p','d' for next, previous, default
    #' @param pmax    maximum scale for vertical axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotMean = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      eps <- private$y_stoch_args[[6]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      means <- private$G #protect against recursive call
      if(is.null(means)) { means <- self$Mean(who="A")[[1]] }
      timeeps <- private$Gteps #protect against recursive call
      if(is.null(timeeps)) { timeeps <- self$Mean(who="A")[[2]] }
      densities <- private$p #no plot or copy
      if(is.null(densities)) { densities <- self$Density(who="A")[[1]] }
      probabilities <- private$P #no plot or copy
      if(is.null(probabilities)) { probabilities <- self$Probability(who="A")[[1]] }
      meanteps <- x
      if(rho > 0) { meanteps <- mu+(x-mu)*exp(-rho*(timeeps-s)) }
      m <- length(t)
      n <- length(y)
      if(n > 1) { dy <- (y[n]-y[1])/(n-1) }
      else { dy <- 1 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        probabilities <- probabilities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        probabilities <- probabilities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",""),c("Mean",""),c("s",s),c("x",x),c("rho",rho),c("mu",mu),c("eps",eps),c("Gteps",timeeps),c("t","G"),cbind(t,means))
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>e</i>=",esym,format(eps,digits=4),")",esml)
      if(labels == TRUE)
      {
        if(is.null(title)) { title <- "Mean" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      #2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>G</i>(<i>t</i>|<i>s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        # OUP_A_MeanToConverge2D
        if(type < -0.5)
        {
          meanline <- list(color=cyn$d,width=4)
          meantepsline <- list(dash="dot",color=cyn$d,width=4)
          tepsline <- list(dash="dot",color=cyn$d,width=4)
          horzline <- list(color=gry$d,width=1)
          if(x < mu) { tepsG <- list(x=timeeps,y=meanteps,text=paste(sep="","&nbsp;<i>t</i><sub>",eps,"</sub>",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>G</i>(<i>t</i><sub>",eps,"</sub>)",bsym,"=",esym,format(meanteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="top",align="left",showarrow=FALSE) }
          else { tepsG <- list(x=timeeps,y=meanteps,text=paste(sep="","&nbsp;<i>t</i><sub>",eps,"</sub",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>G</i>(<i>t</i><sub>",eps,"</sub>)",bsym,"=",esym,format(meanteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE) }
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_MeanToConverge2D")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter",x=t,y=means,mode="lines",line=meanline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),mode="lines",line=horzline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(x,x),mode="lines",line=horzline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],timeeps),y=c(meanteps,meanteps),mode="lines",line=meantepsline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(timeeps,timeeps),y=c(mu,x),mode="lines",line=tepsline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,annotations=tepsG,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_Mean2D
        else
        {
          meanline <- list(color=cyn$d,width=4)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Mean2D")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter",x=t,y=means,mode="lines",line=meanline) %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        meanline <- list(color=cyn$e,width=8)
        if(x < mu) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
        # OUP_A_Mean3DDensity
        if(type < 1.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
          pmeans <- vector("double",m)
          coordinatemeans <- vector("character",m)
          coordinatepmeans <- vector("character",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            pmeans[i] <- private$OUPDensity(s,x,t[i],means[i],rho,mu,sigma,dy)
            coordinatepmeans[i] <- paste(sep="","<i>p</i>(<i>t,G</i>)=",format(pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>p</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>y</i>=",format(y[j],digits=4)) }
          }
          pmesh <- MeshCurtainSmooth(means,t,pmeans,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,range=c(-0.03*pmax,1.03*pmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          pmeanline <- list(color=cyn$d,width=8)
          densityline <- list(color=blu$d,width=8)
          gradientpmeans <- list(c(0,cyn$b),c(1,cyn$b))
          gradient <- list(c(0,blu$c),c(1,blu$c))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Mean3DDensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=pmeans,name="<i>p</i>(<i>t,G</i>)",mode="lines",line=pmeanline,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans") %>%
            add_trace(.,type="mesh3d",x=pmesh$xvertex,y=pmesh$yvertex,z=pmesh$zvertex,i=pmesh$ivertex,j=pmesh$jvertex,k=pmesh$kvertex,intensity=pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientpmeans,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="pmeans",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly")
          while(i < m)
          {
            i <- i+dt
            lineopacity <- lineopacity-0.07
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly",showlegend=FALSE)
         }
          fig <- add_trace(fig,type="surface",x=y,y=t,z=densities,name="<i>p</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        # OUP_A_Mean3DProbability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
          Pmeans <- vector("double",m)
          coordinatemeans <- vector("character",m)
          coordinatePmeans <- vector("character",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            Pmeans[i] <- private$OUPProbability(s,x,t[i],means[i],rho,mu,sigma,psi)
            coordinatePmeans[i] <- paste(sep="","<i>P</i>(<i>t,G</i>)=",format(Pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>y</i>=",y[j]) }
          }
          Pmesh <- MeshCurtainSmooth(means,t,Pmeans,rep(0,m))
          ygap <- 0.03*(y[n]-y[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          Pmeanline <- list(color=cyn$d,width=8)
          probabilityline <- list(color=grn$d,width=8)
          gradientPmeans <- list(c(0,cyn$b),c(1,cyn$b))
          gradient <- list(c(0,grn$c),c(1,grn$c))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Mean3DProbability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=Pmeans,name="<i>P</i>(<i>t,G</i>)",mode="lines",line=Pmeanline,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans") %>%
            add_trace(.,type="mesh3d",x=Pmesh$xvertex,y=Pmesh$yvertex,z=Pmesh$zvertex,i=Pmesh$ivertex,j=Pmesh$jvertex,k=Pmesh$kvertex,intensity=Pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientPmeans,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Pmeans",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly")
          while(i < m)
          {
            i <- i+dt
            lineopacity <- lineopacity-0.07
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="surface",x=y,y=t,z=probabilities,name="<i>P</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot variances
    #' @param type  = -1,...,2 or 'n','p','d' for next, previous, default
    #' @param pmax    maximum scale for vertical axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotVariance = function(type=NULL,pmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,3)[[1]]
      self$set_plot_args(pmax,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      eps <- private$y_stoch_args[[6]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      blu <- private$plot_colors$blu
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      means <- private$G #no plot or copy
      if(is.null(means)) { means <- self$Mean(who="A")[[1]] }
      variances <- private$H2 #protect against recursive call
      if(is.null(variances)) { variances <- self$Variance(who="A")[[1]] }
      timeeps <- private$H2teps #protect against recursive call
      if(is.null(timeeps)) { timeeps <- self$Variance(who="A")[[2]] }
      densities <- private$p #no plot or copy
      if(is.null(densities)) { densities <- self$Density(who="A")[[1]] }
      probabilities <- private$P #no plot or copy
      if(is.null(probabilities)) { probabilities <- self$Probability(who="A")[[1]] }
      sqrts <- variances^0.5
      asymvar <- sigma^2/(2*rho)
      varianceteps <- Inf
      meansplus <- means+sqrts
      meansminus <- means-sqrts
      if(rho > 0) { varianceteps <- sigma^2/(2*rho)*(1-exp(-2*rho*(timeeps-s))) }
      m <- length(t)
      n <- length(y)
      if(n > 1) { dy <- (y[n]-y[1])/(n-1) }
      else { dy <- 1 }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        meansplus <- meansplus[Ixbeg:Ixend]
        meansminus <- meansminus[Ixbeg:Ixend]
        variances <- variances[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,,drop=FALSE]
        probabilities <- probabilities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend,drop=FALSE]
        probabilities <- probabilities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",""),c("Variance",""),c("s",s),c("rho",rho),c("sigma",sigma),c("eps",eps),c("H2teps",timeeps),c("t","H\u00B2"),cbind(t,variances))
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>e</i>=",esym,format(eps,digits=4),")",esml)
      if(labels == TRUE)
      {
        if(is.null(title)) { title <- "Variance" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # 2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>H</i><sup>2</sup>(<i>t</i>|<i>s</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        # OUP_A_VarianceToConverge2D
        if(type < -0.5)
        {
          varianceline <- list(color=mgn$d,width=4)
          variancetepsline <- list(dash="dot",color=mgn$d,width=4)
          tepsline <- list(dash="dot",color=mgn$d,width=4)
          horzline <- list(color=gry$d,width=1)
          tepsH2 <- list(x=timeeps,y=varianceteps,text=paste(sep="","&nbsp;<i>t</i><sub>",eps,"</sub>",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>H</i><sup>2</sup>(<i>t</i><sub>",eps,"</sub>)",bsym,"=",esym,format(varianceteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="top",align="left",showarrow=FALSE)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_VarianceToConverge2D")
          fig <- plot_ly()  %>%
            add_trace(.,type="scatter",x=t,y=variances,mode="lines",line=varianceline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(asymvar,asymvar),mode="lines",line=horzline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],timeeps),y=c(varianceteps,varianceteps),mode="lines",line=variancetepsline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(timeeps,timeeps),y=c(0,asymvar),mode="lines",line=tepsline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,annotations=tepsH2,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_Variance2D
        else
        {
          varianceline <- list(color=mgn$d,width=4)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Variance2D")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter",x=t,y=variances,mode="lines",line=varianceline) %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        plusline <- list(color=mgn$e,width=8)
        meanline <- list(color=cyn$e,width=8)
        minusline <- list(color=mgn$e,width=8)
        if(x < mu) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
        # OUP_A_Variance3DDensity
        if(type < 1.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
          pmeansplus <- vector("double",m)
          pmeans <- vector("double",m)
          pmeansminus <- vector("double",m)
          coordinatemeansplus <- vector("character",m)
          coordinatemeans <- vector("character",m)
          coordinatemeansminus <- vector("character",m)
          coordinatepmeansplus <- vector("character",m)
          coordinatepmeans <- vector("character",m)
          coordinatepmeansminus <- vector("character",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeansplus[i] <- paste(sep="","<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)=",format(meansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeansminus[i] <- paste(sep="","<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)=",format(meansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            pmeansplus[i] <- private$OUPDensity(s,x,t[i],meansplus[i],rho,mu,sigma,dy)
            coordinatepmeansplus[i] <- paste(sep="","<i>p</i>(<i>t,G</i>+<i>H</i>)=",format(pmeansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G+H</i>=",format(meansplus[i],digits=4))
            pmeans[i] <- private$OUPDensity(s,x,t[i],means[i],rho,mu,sigma,dy)
            coordinatepmeans[i] <- paste(sep="","<i>p</i>(<i>t,G</i>)=",format(pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            pmeansminus[i] <- private$OUPDensity(s,x,t[i],meansminus[i],rho,mu,sigma,dy)
            coordinatepmeansminus[i] <- paste(sep="","<i>p</i>(<i>t,G</i>-<i>H</i>)=",format(pmeansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G-H</i>=",format(meansminus[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>p</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>y</i>=",format(y[j],digits=4)) }
          }
          pmeshplus <- MeshCurtainSmooth(meansplus,t,pmeansplus,rep(0,m))
          pmesh <- MeshCurtainSmooth(means,t,pmeans,rep(0,m))
          pmeshminus <- MeshCurtainSmooth(meansminus,t,pmeansminus,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,range=c(-0.03*pmax,1.03*pmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          pplusline <- list(color=mgn$d,width=8)
          pmeanline <- list(color=cyn$d,width=8)
          pminusline <- list(color=mgn$d,width=8)
          densityline <- list(color=blu$d,width=8)
          gradientplus <- list(c(0,mgn$b),c(1,mgn$b))
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientminus <- list(c(0,mgn$b),c(1,mgn$b))
          gradient <- list(c(0,blu$c),c(1,blu$c))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Variance3DDensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)",mode="lines",line=plusline,hoverinfo="text",text=coordinatemeansplus) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)",mode="lines",line=minusline,hoverinfo="text",text=coordinatemeansminus) %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=pmeansplus,name="<i>p</i>(<i>G</i>+<i>H</i>)",mode="lines",line=pplusline,hoverinfo="text",text=coordinatepmeansplus,legendgroup="pmeansplus") %>%
            add_trace(.,type="mesh3d",x=pmeshplus$xvertex,y=pmeshplus$yvertex,z=pmeshplus$zvertex,i=pmeshplus$ivertex,j=pmeshplus$jvertex,k=pmeshplus$kvertex,intensity=pmeshplus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientplus,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="pmeansplus",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=pmeans,name="<i>p</i>(<i>t,G</i>)",mode="lines",line=pmeanline,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans") %>%
            add_trace(.,type="mesh3d",x=pmesh$xvertex,y=pmesh$yvertex,z=pmesh$zvertex,i=pmesh$ivertex,j=pmesh$jvertex,k=pmesh$kvertex,intensity=pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="pmeans",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=pmeansminus,name="<i>p</i>(<i>G</i>-<i>H</i>)",mode="lines",line=pminusline,hoverinfo="text",text=coordinatepmeansminus,legendgroup="pmeansminus") %>%
            add_trace(.,type="mesh3d",x=pmeshminus$xvertex,y=pmeshminus$yvertex,z=pmeshminus$zvertex,i=pmeshminus$ivertex,j=pmeshminus$jvertex,k=pmeshminus$kvertex,intensity=pmeshminus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientminus,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="pmeansminus",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly")
          while(i < m)
          {
            i <- i+dt
            lineopacity <- lineopacity-0.07
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="surface",x=y,y=t,z=densities,name="<i>p</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        # OUP_A_Variance3DProbability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
          Pmeans <- rep(0.5,m)
          if(psi > 0)
          {
            Pmeansplus <- rep(0.1586553,m)
            Pmeansminus <- rep(0.8413447,m)
          }
          else
          {
            Pmeansplus <- rep(0.8413447,m)
            Pmeansminus <- rep(0.1586553,m)
          }
          coordinatemeansplus <- vector("character",m)
          coordinatemeans <- vector("character",m)
          coordinatemeansminus <- vector("character",m)
          coordinatePmeansplus <- vector("character",m)
          coordinatePmeans <- vector("character",m)
          coordinatePmeansminus <- vector("character",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeansplus[i] <- paste(sep="","<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)=",format(meansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeansminus[i] <- paste(sep="","<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)=",format(meansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatePmeansplus[i] <- paste(sep="","<i>P</i>(<i>t,G+<i>H</i>)=",format(Pmeansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G+H</i>=",format(meansplus[i],digits=4))
            coordinatePmeans[i] <- paste(sep="","<i>P</i>(<i>t,G</i>)=",format(Pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            coordinatePmeansminus[i] <- paste(sep="","<i>P</i>(<i>t,G-<i>H</i>)=",format(Pmeansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G-H</i>=",format(meansminus[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(probabilities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
          }
          Pmeshplus <- MeshCurtainSmooth(meansplus,t,Pmeansplus,rep(0,m))
          Pmesh <- MeshCurtainSmooth(means,t,Pmeans,rep(0,m))
          Pmeshminus <- MeshCurtainSmooth(meansminus,t,Pmeansminus,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(1.03*y[1]-0.03*y[n],1.03*y[n]-0.03*y[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          Pplusline <- list(color=mgn$d,width=8)
          Pmeanline <- list(color=cyn$d,width=8)
          Pminusline <- list(color=mgn$d,width=8)
          probabilityline <- list(color=grn$d,width=8)
          gradientplus <- list(c(0,mgn$b),c(1,mgn$b))
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientminus <- list(c(0,mgn$b),c(1,mgn$b))
          gradient <- list(c(0,grn$c),c(1,grn$c))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Variance3DProbability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)",mode="lines",line=plusline,hoverinfo="text",text=coordinatemeansplus) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)",mode="lines",line=minusline,hoverinfo="text",text=coordinatemeansminus) %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=Pmeansplus,name="<i>P</i>(<i>t,G</i>+<i>H</i>)",mode="lines",line=Pplusline,hoverinfo="text",text=coordinatePmeansplus,legendgroup="Pmeansplus") %>%
            add_trace(.,type="mesh3d",x=Pmeshplus$xvertex,y=Pmeshplus$yvertex,z=Pmeshplus$zvertex,i=Pmeshplus$ivertex,j=Pmeshplus$jvertex,k=Pmeshplus$kvertex,intensity=Pmeshplus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientplus,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Pmeansplus",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=Pmeans,name="<i>P</i>(<i>t,G</i>)",mode="lines",line=Pmeanline,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans") %>%
            add_trace(.,type="mesh3d",x=Pmesh$xvertex,y=Pmesh$yvertex,z=Pmesh$zvertex,i=Pmesh$ivertex,j=Pmesh$jvertex,k=Pmesh$kvertex,intensity=Pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Pmeans",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=Pmeansminus,name="<i>P</i>(<i>t,G</i>-<i>H</i>)",mode="lines",line=Pminusline,hoverinfo="text",text=coordinatePmeansminus,legendgroup="Pmeansminus") %>%
            add_trace(.,type="mesh3d",x=Pmeshminus$xvertex,y=Pmeshminus$yvertex,z=Pmeshminus$zvertex,i=Pmeshminus$ivertex,j=Pmeshminus$jvertex,k=Pmeshminus$kvertex,intensity=Pmeshminus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientminus,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Pmeansminus",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly")
          while(i < m)
          {
            i <- i+dt
            lineopacity <- lineopacity-0.07
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="surface",x=y,y=t,z=probabilities,name="<i>P</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot transition densities
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param pmax    maximum scale for vertical axis
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
      type <- self$set_plot_type(type,4)[[1]]
      self$set_plot_args(pmax,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      pmax <- private$plot_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      blu <- private$plot_colors$blu
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      densities <- private$p #protect against recursive call
      if(is.null(densities)) { densities <- self$Density(who="A")[[1]] }
      m <- length(t)
      n <- length(y)
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
        clip <- rbind(c("Analytical",rep("",n)),c("Transition Densities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("p(t,y)",y),cbind(t,densities))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Transition Densities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Density2D
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
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        densityline <- list(color=blu$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Density2D")
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=y,y=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
          lineopacity <- lineopacity-0.07
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_Density3D
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
        densityline <- list(color=blu$e,width=8)
        gradient <- list(c(0,blu$c),c(1,blu$c))
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Density3D")
        fig <- plot_ly()
        tt <- vector("double",n)
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        i <- 1
        for(j in 1:n) { tt[j] <- t[i] }
        fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p")
        while(i < m)
        {
          i <- i+dt
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=densities,name="<i>p</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly", showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot transition probabilities
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotProbability = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(psi <= 0) { probabilities <- private$Pneg } #protect against recursive call
      else { probabilities <- private$Ppos }
      if(is.null(probabilities)) { probabilities <- self$Probability(who="A")[[1]] }
      m <- length(t)
      n <- length(y)
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        probabilities <- probabilities[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        probabilities <- probabilities[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n)),c("Transition Probabilities",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("P(t,y)",y),cbind(t,probabilities))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>y</i>=",esym,format(psi,digits=4),")",esml)
        if(is.null(title)) { title <- "Transition Probabilities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Probability2D
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
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        probabilityline <- list(color=grn$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Probability2D")
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=y,y=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
          lineopacity <- lineopacity-0.07
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_Probability3D
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
        probabilityline <- list(color=grn$e,width=8)
        gradient <- list(c(0,grn$c),c(1,grn$c))
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Probability3D")
        fig <- plot_ly()
        tt <- vector("double",n)
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        i <- 1
        for(j in 1:n) { tt[j] <- t[i] }
        fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P")
        while(i < m)
        {
          i <- i+dt
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=probabilities,name="<i>P</i>(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly", showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot double integrals of transition densities
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param ybeg    begin value for state axis
    #' @param yend    end value for state axis
    #' @return plot
    PlotDoubleIntegral = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,ybeg=NULL,yend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(psi <= 0 ) { doubleintegrals <- private$PPneg } #protect against recursive call
      else { doubleintegrals <- private$PPpos }
      if(is.null(doubleintegrals)) { doubleintegrals <- self$DoubleIntegral(who="A")[[1]] }
      m <- length(t)
      n <- length(y)
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[,Ixbeg:Ixend,drop=FALSE]
        n <- length(y)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n)),c("Double Integrals",rep("",n)),c("s",s,rep("",n-1)),c("x",x,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("psi",psi,rep("",n-1)),c("\u2119(t,y)",y),cbind(t,doubleintegrals))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>y</i>=",esym,format(psi,digits=4),")",esml)
        if(is.null(title)) { title <- "Double Integrals" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_DoubleIntegral2D
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
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        doubleintegralline <- list(color=red$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_DoubleIntegral2D")
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=y,y=doubleintegrals[i,],mode="lines",line=doubleintegralline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
          lineopacity <- lineopacity-0.07
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_DoubleIntegral3D
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
        doubleintegralline <- list(color=red$e,width=8)
        gradient <- list(c(0,red$c),c(1,red$c))
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_DoubleIntegral3D")
        fig <- plot_ly()
        tt <- vector("double",n)
        dt <- as.integer((m-1)/10)
        if(dt < 1) { dt <- 1 }
        i <- 1
        for(j in 1:n) { tt[j] <- t[i] }
        fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=doubleintegrals[i,],name="\u2119(<i>t,y</i>)",mode="lines",line=doubleintegralline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="PP")
        while(i < m)
        {
          i <- i+dt
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=doubleintegrals[i,],mode="lines",line=doubleintegralline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="PP",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=y,y=t,z=doubleintegrals,name="\u2119(<i>t,y</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly", showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot options
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param sbeg    begin value for time axis
    #' @param send    end value for time axis
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotOption = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,sbeg=NULL,send=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(phi <= 0) { options <- private$OOneg } #protect against recursive call
      else { options <- private$OOpos }
      if(is.null(options)) { options <- self$Option(who="A")[[1]] }
      m <- length(s)
      n <- length(x)
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m || Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        options <- options[Ixend:Ixbeg,,drop=FALSE]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Options",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("c",c,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options)) }
        else { clip <- rbind(c("Analytical",rep("",n)),c("Options",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("\uD835\uDD46(s,x)",x),cbind(s,options)) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        if(phi > 0) { syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>b</i>",bsym,"=",esym,format(b,digits=4),")",esml) }
        else { syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>c</i>",bsym,"=",esym,format(c,digits=4),")",esml) }
        if(is.null(title)) { title <- "Options" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Option2D
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
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        optionline <- list(color=red$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Option2D")
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=x,y=options[i,],mode="lines",line=optionline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+ds
          lineopacity <- lineopacity-0.07
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_Option3D
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
        optionline <- list(color=red$d,width=8)
        gradient <- list(c(0,red$c),c(1,red$c))
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Option3D")
        fig <- plot_ly()
        ss <- vector("double",n)
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        i <- 1
        for(j in 1:n) { ss[j] <- s[i] }
        fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],name="\uD835\uDD46(<i>s,x</i>)",mode="lines",line=optionline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="OO")
        while(i < m)
        {
          i <- i+ds
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { ss[j] <- s[i] }
          fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],mode="lines",line=optionline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="OO",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly", showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the option envelope
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param sbeg    begin value for time axis
    #' @param send    end value for time axis
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotOptionEnvelope = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,sbeg=NULL,send=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      ylw <- private$plot_colors$ylw
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(phi <= 0)
      {
        OOhat <- private$OOhatneg #protect against recursive call
        if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="A")[[1]] }
        shat <- private$shatneg
        options <- private$OOneg #no plot or copy
        if(is.null(options)) { options <- self$Option(who="A")[[1]] }
        dOOdsconvex <- private$dOOdsconvexneg
        if(is.null(dOOdsconvex)) { dOOdsconvex <- private$OUPdOOdsZero()[[1]] }
        dOOdsconcave <- private$dOOdsconcaveneg
        dOOdspatch <- private$dOOdspatchneg
      }
      else
      {
        OOhat <- private$OOhatpos
        if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="A")[[1]] }
        shat <- private$shatpos
        options <- private$OOpos
        if(is.null(options)) { options <- self$Option(who="A")[[1]] }
        dOOdsconvex <- private$dOOdsconvexpos
        if(is.null(dOOdsconvex)) { dOOdsconvex <- private$OUPdOOdsZero()[[1]] }
        dOOdsconcave <- private$dOOdsconcavepos
        dOOdspatch <- private$dOOdspatchpos
      }
      m <- length(s)
      n <- length(x)
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m || Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        options <- options[Ixend:Ixbeg,,drop=FALSE]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        OOhat <- OOhat[Ixbeg:Ixend]
        shat <- shat[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend,drop=FALSE]
        dOOdsconvex <- dOOdsconvex[,Ixbeg:Ixend,drop=FALSE]
        dOOdsconcave <- dOOdsconcave[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Option Envelope",rep("",n)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("c",c,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat)) }
        else { clip <- rbind(c("Analytical",rep("",n)),c("Option Envelope",rep("",n)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat)) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        if(phi > 0) { syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>b</i>",bsym,"=",esym,format(b,digits=4),")",esml) }
        else { syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>c</i>",bsym,"=",esym,format(c,digits=4),")",esml) }
        if(is.null(title)) { title="Option Envelope"  }
      }
      else if(is.null(title)) { title="" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Envelope2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        if(is.null(yaxis)) { yaxis <- "\u00D4(<i>x</i>|<i>y</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        terminal <- vector("double",n)
        for(j in 1:n)
        {
          if(phi > 0)
          {
            if(x[j] > y) { terminal[j] <- x[j]-y+b }
            else { terminal[j] <- b}
          }
          else
          {
            if(y > x[j]) { terminal[j] <- y-x[j]-c }
            else { terminal[j] <- -c }
          }
        }
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(phi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        OOhatline <- list(color=red$d,width=4)
        terminalline <- list(color=gry$c,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Envelope2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=terminal,mode="lines",line=terminalline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=OOhat,mode="lines",line=OOhatline,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_OptionEnvelope3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis)) { zaxis <- "\u00D4(<i>x</i>|<i>y</i>)" }
        OOhold <- vector("double",n)
        OOexercise <- vector("double",n)
        coordinateshat <- vector("character",n)
        coordinatesconvex <- vector("character",n)
        coordinatesconcave <- vector("character",n)
        coordinates <- matrix("",m,n)
        finite <- TRUE
        for(j in 1:n)
        {
          if(is.finite(OOhat[j]))
          {
            if(shat[j] == t)
            {
              OOhold[j] <- NA
              OOexercise[j] <- OOhat[j]
            }
            else
            {
              OOhold[j] <- OOhat[j]
              OOexercise[j] <- NA
            }
            coordinateshat[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(OOhat[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
            coordinatesconvex[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(dOOdsconvex[1,j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
            coordinatesconcave[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(dOOdsconcave[1,j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          }
          else { finite <- FALSE }
          for(i in 1:m) { coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
        }
        if(phi > 0) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        OOmax <- max(options[1,1],options[1,n])
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*x[1]-0.03*x[n],1.03*x[n]-0.03*x[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*s[m]-0.03*s[1],1.03*s[1]-0.03*s[m]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,range=c(-0.03*OOmax,1.03*OOmax),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview, zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_OptionEnvelope3D")
        gradient <- list(c(0,red$c),c(1,red$c))
        fig <- plot_ly()
        if(finite == TRUE)
        {
          OOholdmesh <- MeshCurtainSmooth(x,shat,OOhold,rep(0,n))
          OOexercisemesh <- MeshCurtainSmooth(x,shat,OOexercise,rep(0,n))
          dzero1mesh <- MeshCurtainSmooth(x,dOOdsconvex[2,],dOOdsconvex[1,],rep(0,n))
          dzero2mesh <- MeshCurtainSmooth(x,dOOdsconcave[2,],dOOdsconcave[1,],rep(0,n))
          dzero3mesh <- MeshCurtainSmooth(dOOdspatch[3,],dOOdspatch[2,],dOOdspatch[1,],rep(0,3))
          OOholdline <- list(color=ylw$e,width=10)
          OOexerciseline <- list(color=ylw$e,width=10)
          OOhatline <- list(dash="dash",color=ylw$d,width=10)
          dOOdszeroline <- list(color=red$e,width=8)
          gradienthat <- list(c(0,ylw$d),c(1,ylw$d))
          gradientzero <- list(c(0,red$d),c(1,red$d))
          fig <- add_trace(fig,type="scatter3d",x=x,y=shat,z=OOhat,name="\u00D4(<i>x</i>)",mode="lines",line=OOhatline,hoverinfo="text",text=coordinateshat,legendgroup="OOhat") %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=OOhold,mode="lines",line=OOholdline,hoverinfo="text",text=coordinateshat,legendgroup="OOhat",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=OOexercise,mode="lines",line=OOexerciseline,hoverinfo="text",text=coordinateshat,legendgroup="OOhat",showlegend=FALSE) %>%
            add_trace(.,type="mesh3d",x=OOholdmesh$xvertex,y=OOholdmesh$yvertex,z=OOholdmesh$zvertex,i=OOholdmesh$ivertex,j=OOholdmesh$jvertex,k=OOholdmesh$kvertex,intensity=OOholdmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthat,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="OOhat",showlegend=FALSE) %>%
            add_trace(.,type="mesh3d",x=OOexercisemesh$xvertex,y=OOexercisemesh$yvertex,z=OOexercisemesh$zvertex,i=OOexercisemesh$ivertex,j=OOexercisemesh$jvertex,k=OOexercisemesh$kvertex,intensity=OOexercisemesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthat,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="OOhat",showlegend=FALSE) %>%
            add_trace(.,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
            add_trace(.,type="scatter3d",x=x,y=dOOdsconvex[2,],z=dOOdsconvex[1,],name="d\uD835\uDD46/ds=0",mode="lines",line=dOOdszeroline,hoverinfo="text",text=coordinatesconvex,legendgroup="dOOds",visible="legendonly") %>%
            add_trace(.,type="scatter3d",x=x,y=dOOdsconcave[2,],z=dOOdsconcave[1,],mode="lines",line=dOOdszeroline,hoverinfo="text",text=coordinatesconcave,legendgroup="dOOds",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=dOOdspatch[3,],y=dOOdspatch[2,],z=dOOdspatch[1,],mode="lines",line=dOOdszeroline,hoverinfo="skip",legendgroup="dOOds",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="mesh3d",x=dzero1mesh$xvertex,y=dzero1mesh$yvertex,z=dzero1mesh$zvertex,i=dzero1mesh$ivertex,j=dzero1mesh$jvertex,k=dzero1mesh$kvertex,intensity=dzero1mesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientzero,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="dOOds",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="mesh3d",x=dzero2mesh$xvertex,y=dzero2mesh$yvertex,z=dzero2mesh$zvertex,i=dzero2mesh$ivertex,j=dzero2mesh$jvertex,k=dzero2mesh$kvertex,intensity=dzero2mesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientzero,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="dOOds",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="mesh3d",x=dzero3mesh$xvertex,y=dzero3mesh$yvertex,z=dzero3mesh$zvertex,i=dzero3mesh$ivertex,j=dzero3mesh$jvertex,k=dzero3mesh$kvertex,intensity=dzero3mesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientzero,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="dOOds",visible="legendonly",showlegend=FALSE)
        }
        else { fig <- add_trace(fig,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,showlegend=TRUE) }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the decision threshold
    #' @param type  = 0
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDecisionThreshold = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      if(phi <= 0) { decision <- private$KOOneg } #protect against recursive call
      else { decision <- private$KOOpos }
      if(is.null(decision))
      {
        decision <- self$DecisionThreshold(who="A")
        k <- decision[[1]]
        OO <- decision[[2]]
      }
      else
      {
        k <- decision[1]
        OO <- decision[2]
      }
      if(phi <= 0) { OOhat <- private$OOhatneg } #no plot or copy
      else { OOhat <- private$OOhatpos }
      if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="A")[[1]] }
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        OOhat <- OOhat[Ixbeg:Ixend]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(phi <= 0) { clip <- rbind(c("Analytical",""),c("Decision Threshold",""),c("y",y),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("c",c),c("k",k),c("\u00D4",OO)) }
        else { clip <- rbind(c("Analytical",""),c("Decision Threshold",""),c("y",y),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("b",b),c("k",k),c("\u00D4",OO)) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_A_Decision2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        if(phi > 0) { syms <- paste(sep="",bsml,"(<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>b</i>",bsym,"=",esym,format(b,digits=4),")",esml) }
        else { syms <- paste(sep="",bsml,"(<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>c</i>",bsym,"=",esym,format(c,digits=4),")",esml) }
        if(is.null(title)) { title <- "Decision Threshold" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "\u00D4(<i>x</i>|<i>y</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      terminal <- vector("double",n)
      for(j in 1:n)
      {
        if(phi > 0)
        {
          if(x[j] > y) { terminal[j] <- x[j]-y+b }
          else { terminal[j] <- b}
        }
        else
        {
          if(y > x[j]) { terminal[j] <- y-x[j]-c }
          else { terminal[j] <- -c }
        }
      }
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      if(phi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
      else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
      OOhatline <- list(color=red$d,width=4)
      terminalline <- list(color=gry$c,width=4)
      OOline <- list(dash="dot",color=red$d,width=4)
      kline <- list(dash="dot",color=red$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Decision2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=terminal,mode="lines",line=terminalline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=x,y=OOhat,mode="lines",line=OOhatline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(k,k),y=c(0,OO),mode="lines",line=kline,hoverinfo="x+y")
      if(phi > 0)
      {
        fig <- add_trace(fig,type="scatter",x=c(x[n],k),y=c(OO,OO),mode="lines",line=OOline,hoverinfo="x+y")
        KOO <- list(x=k,y=OO,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(OO,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
      }
      else
      {
        fig <- add_trace(fig,type="scatter",x=c(x[1],k),y=c(OO,OO),mode="lines",line=OOline,hoverinfo="x+y")
        KOO <- list(x=k,y=OO,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(OO,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
      }
      fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,annotations=KOO,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot obligations
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param sbeg    begin value for time axis
    #' @param send    end value for time axis
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotObligation = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,sbeg=NULL,send=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,4)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      ylw <- private$plot_colors$ylw
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(phi <= 0) { obligations <- private$BCneg } #protect against recursive call
      else { obligations <- private$BCpos }
      if(is.null(obligations)) { obligations <- self$Obligation(who="A")[[1]] }
      if(phi <= 0)
      {
        OOhat <- private$OOhatneg #protect against recursive call
        if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="A")[[1]] }
        shat <- private$shatneg
      }
      else
      {
        OOhat <- private$OOhatpos
        if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="A")[[1]] }
        shat <- private$shatpos
      }
      m <- length(s)
      n <- length(x)
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m || Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        obligations <- obligations[Ixend:Ixbeg,,drop=FALSE]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        OOhat <- OOhat[Ixbeg:Ixend]
        shat <- shat[Ixbeg:Ixend]
        obligations <- obligations[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(phi <= 0) { clip <- rbind(c("Analytical",rep("",n)),c("Obligation",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("c",c,rep("",n-1)),c("\uD835\uDD39(s,x)",x),cbind(s,obligations)) }
        else { clip <- rbind(c("Analytical",rep("",n)),c("Prohibition",rep("",n)),c("t",t,rep("",n-1)),c("y",y,rep("",n-1)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("r",r,rep("",n-1)),c("phi",phi,rep("",n-1)),c("b",b,rep("",n-1)),c("c",c,rep("",n-1)),c("\u2102(s,x)",x),cbind(s,obligations)) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>t</i>",bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),",<i>b</i>",bsym,"=",esym,format(b,digits=4),",<i>c</i>",bsym,"=",esym,format(c,digits=4),")",esml)
        if(is.null(title))
        {
          if(phi > 0) { title <- "Prohibition" }
          else { title <- "Obligation" }
        }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Obligation2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        if(is.null(yaxis))
        {
          if(phi > 0) { yaxis <- "\u2102(<i>s,x</i>|<i>t,y</i>)" }
          else { yaxis <- "\uD835\uDD39(<i>s,x</i>|<i>t,y</i>)" }
        }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(phi > 0) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=TRUE) }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=TRUE,side="right") }
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        obligationline <- list(color=ylw$d,width=4)
        lineopacity <- 1
        fig <- plot_ly()
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Obligation2D")
        i <- 1
        if(phi > 0)
        {
          while(i <= m)
          {
            fig <- add_trace(fig,type="scatter",x=x,y=obligations[i,],mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="x+y")
            i <- i+ds
            lineopacity <- lineopacity-0.07
          }
        }
        else
        {
          while(i <= m)
          {
            fig <- add_trace(fig,type="scatter",x=x,y=obligations[i,],mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="x+y")
            i <- i+ds
            lineopacity <- lineopacity-0.07
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_Obligation3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis))
        {
          if(phi > 0) { zaxis <- "\u2102(<i>s,x</i>|<i>t,y</i>)" }
          else { zaxis <- "\uD835\uDD39(<i>s,x</i>|<i>t,y</i>)" }
        }
        origin <- matrix(0.0,m,n)
        originx <- rep(0.0,n)
        originy <- rep(0.0,m)
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            if(phi > 0) { coordinates[i,j] <- paste(sep="","\u2102(<i>s,x</i>)=",format(obligations[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
            else { coordinates[i,j] <- paste(sep="","\uD835\uDD39(<i>s,x</i>)=",format(obligations[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
          }
        }
        if(phi > 0) { tracename <- "\u2102(<i>s,x</i>)" }
        else { tracename <- "\uD835\uDD39(<i>s,x</i>)" }
        if(phi > 0) { spy <- list(x=0.6,y=-2.3,z=0.3) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.3) }
        else { spy <- list(x=-0.6,y=-2.3,z=0.3) }
        BCmax = max(obligations[1,1],obligations[1,n])
        BCmin = min(obligations[1,1],obligations[1,n])
        xview <- list(title=xaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$a,showbackground=walls,range=c(1.03*x[1]-0.03*x[n],1.03*x[n]-0.03*x[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$a,showbackground=walls,range=c(1.03*s[m]-0.03*s[1],1.03*s[1]-0.03*s[m]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$b,showbackground=floor,range=c(1.03*BCmin-0.03*BCmax,1.03*BCmax-0.03*BCmin),tickmode="auto",nticks=5,zeroline=TRUE,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        obligationline <- list(color=ylw$d,width=8)
        obgradient <- list(c(0,ylw$c),c(1,ylw$b))
        zgradient <- list(c(0,ylw$b),c(1,ylw$b))
        originframe <- list(color=ylw$c,width=2)
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Obligation3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=x,y=s,z=origin,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=zgradient,reversescale=reverse,opacity=0.5,hoverinfo="skip",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=rep(s[1],n),z=originx,mode="lines",line=originframe,hoverinfo="skip",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=rep(s[m],n),z=originx,mode="lines",line=originframe,hoverinfo="skip",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=rep(x[1],m),y=s,z=originy,mode="lines",line=originframe,hoverinfo="skip",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=rep(x[n],m),y=s,z=originy,mode="lines",line=originframe,hoverinfo="skip",showlegend=FALSE)
        # # scatter
        ss <- vector("double",n)
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        i <- 1
        for(j in 1:n) { ss[j] <- s[i] }
        fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=obligations[i,],name=tracename,mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="BC")
        while(i < m)
        {
          i <- i+ds
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { ss[j] <- s[i] }
          fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=obligations[i,],mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="BC",showlegend=FALSE)
        }
        # # surface
        fig <- add_trace(fig,type="surface",x=x,y=s,z=obligations,name=tracename,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly", showlegend=TRUE)
        # envelope
        obhat <- vector("double",n)
        obhold <- vector("double",n)
        obexercise <- vector("double",n)
        obinterest <- vector("double",n)
        ophat <- vector("double",n)
        ophold <- vector("double",n)
        opexercise <- vector("double",n)
        opinterest <- vector("double",n)
        zero <- rep(0,n)
        coordinatesobhat <- vector("character",n)
        coordinatesophat <- vector("character",n)
        coordinateszero <- vector("character",n)
        for(j in 1:n)
        {
          if(shat[j] < s[m])
          {
            obhat[j] <- NA
            ophat[j] <- NA
            zero[j] <- NA
          }
          else
          {
            G <- mu+(x[j]-mu)*exp(-rho*(t-shat[j]))
            if(phi > 0) { obhat[j] <- exp(-r*(t-shat[j]))*(y-G-b-c) }
            else { obhat[j] <- exp(-r*(t-shat[j]))*(G-y+b+c) }
            ophat[j] <- OOhat[j]+obhat[j]
          }
          coordinatesobhat[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(obhat[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          coordinatesophat[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(ophat[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          coordinateszero[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(zero[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
        }
        if(phi > 0)
        {
          j <- 0
          while(j < n)
          {
            j <- j+1
            if(shat[j] == t)
            {
              obhold[j] <- NA
              obexercise[j] <- obhat[j]
              ophold[j] <- NA
              opexercise[j] <- ophat[j]
              obinterest[j] <- NA
              opinterest[j] <- NA
              if(j > 1)
              {
                obhat[j] <- obhat[j-1]
                ophat[j] <- ophat[j-1]
                opinterest[j] <- obhat[j]-ophat[j]
              }
            }
            else
            {
              obhold[j] <- obhat[j]
              obexercise[j] <- NA
              ophold[j] <- ophat[j]
              opexercise[j] <- NA
              obinterest[j] <- NA
              opinterest[j] <- NA
            }
          }
        }
        else
        {
          j <- n
          while(j > 0)
          {
            if(shat[j] == t)
            {
              obhold[j] <- NA
              obexercise[j] <- obhat[j]
              ophold[j] <- NA
              opexercise[j] <- ophat[j]
              obinterest[j] <- NA
              opinterest[j] <- NA
              if(j < n)
              {
                obhat[j] <- obhat[j+1]
                ophat[j] <- ophat[j+1]
                obinterest[j] <- obhat[j]-ophat[j]
              }
            }
            else
            {
              obhold[j] <- obhat[j]
              obexercise[j] <- NA
              ophold[j] <- ophat[j]
              opexercise[j] <- NA
              obinterest[j] <- NA
              opinterest[j] <- NA
            }
            j <- j-1
          }
        }
        opholdmesh <- MeshCurtainSmooth(x,shat,ophold,rep(0,n))
        obholdmesh <- MeshCurtainSmooth(x,shat,obhold,rep(0,n))
        exercisemesh <- MeshCurtainSmooth(x,shat,opexercise,obexercise)
        obholdline <- list(color=red$c,width=10)
        obexerciseline <- list(color=red$b,width=8)
        obhatline <- list(dash="dash",color=ylw$c,width=8)
        obinterestline <- list(color=ylw$b,width=8)
        opholdline <- list(color=red$c,width=10)
        opexerciseline <- list(color=red$c,width=8)
        ophatline <- list(dash="dash",color=ylw$c,width=8)
        opinterestline <- list(color=ylw$b,width=8)
        zeroline <- list(color=red$c,width=8)
        opgradient <- list(c(0,red$c),c(1,red$c))
        obgradient <- list(c(0,red$b),c(1,red$b))
        fig <- add_trace(fig,type="scatter3d",x=x,y=shat,z=zero,name="\u00D4(<i>x</i>)",mode="lines",line=zeroline,hoverinfo="text",text=coordinateszero,legendgroup="OOhat",visible="legendonly") %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=obhat,mode="lines",line=obhatline,hoverinfo="text",text=coordinatesobhat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=obhold,mode="lines",line=obholdline,hoverinfo="text",text=coordinatesobhat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=obexercise,mode="lines",line=obexerciseline,hoverinfo="text",text=coordinatesobhat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=obinterest,mode="lines",line=obinterestline,hoverinfo="text",text=coordinatesobhat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=ophat,mode="lines",line=ophatline,hoverinfo="text",text=coordinatesophat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=ophold,mode="lines",line=opholdline,hoverinfo="text",text=coordinatesophat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=opexercise,mode="lines",line=opexerciseline,hoverinfo="text",text=coordinatesophat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=opinterest,mode="lines",line=opinterestline,hoverinfo="text",text=coordinatesobhat,legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="mesh3d",x=opholdmesh$xvertex,y=opholdmesh$yvertex,z=opholdmesh$zvertex,i=opholdmesh$ivertex,j=opholdmesh$jvertex,k=opholdmesh$kvertex,intensity=opholdmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=opgradient,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="mesh3d",x=obholdmesh$xvertex,y=obholdmesh$yvertex,z=obholdmesh$zvertex,i=obholdmesh$ivertex,j=obholdmesh$jvertex,k=obholdmesh$kvertex,intensity=obholdmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          add_trace(.,type="mesh3d",x=exercisemesh$xvertex,y=exercisemesh$yvertex,z=exercisemesh$zvertex,i=exercisemesh$ivertex,j=exercisemesh$jvertex,k=exercisemesh$kvertex,intensity=exercisemesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,opacity=0.9,hoverinfo="skip",legendgroup="OOhat",visible="legendonly",showlegend=FALSE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time modes, medians and means
    #' @param type  = -3,...,2 or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum scale for vertical axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotPassageTimeModeMedianMean = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_plot_args(NULL,ptmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      m <- length(t)
      n <- length(z)
      tmodemedianmean <- private$tmodemedianmean
      tmodesmediansmeans <- private$tmodesmediansmeans
      if(is.null(tmodemedianmean) || is.null(tmodesmediansmeans))
      {
        modesmediansmeans <- self$PassageTimeModeMedianMean(who="A")
        tmodemedianmean <- modesmediansmeans[[1]]
        tmodesmediansmeans <- modesmediansmeans[[2]]
      }
      tmode <- tmodemedianmean[[1]]
      tmedian <- tmodemedianmean[[2]]
      tmean <- tmodemedianmean[[3]]
      tmodes <- tmodesmediansmeans[[1]]
      tmedians <- tmodesmediansmeans[[2]]
      tmeans <- tmodesmediansmeans[[3]]
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) || is.null(pt))
      {
        densities <- self$PassageTimeDensity(who="A")
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) || is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(who="A")
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,,drop=FALSE]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        tmodes <- tmodes[Ixbeg:Ixend]
        tmedians <- tmedians[Ixbeg:Ixend]
        tmeans <- tmeans[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend,drop=FALSE]
        Pt <- Pt[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Mode Median Mean",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("tmode(x)"),tmodemedianmean[[1]],paste0("tmode(z)"),tmodesmediansmeans[[1]]),c(paste0("tmedian(x)"),tmodemedianmean[[2]],paste0("tmedian(z)"),tmodesmediansmeans[[2]]),c(paste0("tmean(x)"),tmodemedianmean[[3]],paste0("tmean(z)"),tmodesmediansmeans[[3]]))
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>w</i>=",esym,format(omega,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title))
        {
          if(rho > 0) { title <- "Passage Time Mode, Median and Mean" }
          else { title <- "Passage Time Mode and Median" }
        }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # 2D
      if(type < -1.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        densityline <- list(color=blu$d,width=4,shape="spline",smoothing=1.3)
        probabilityline <- list(color=grn$d,width=4,shape="spline",smoothing=1.3)
        meandashline <- list(color=cyn$d,dash="longdash",width=3)
        meandotline <- list(color=cyn$d,dash="dot",width=1)
        mediandashline <- list(color=grn$d,dash="dash",width=3)
        mediandotline <- list(color=grn$d,dash="dot",width=1)
        modedashline <- list(color=blu$d,dash="dot",width=3)
        modedotline <- list(color=blu$d,dash="dot",width=1)
        fig <- plot_ly()
        # OUP_A_PassageTimeModeMedianMean2Dpt
        if(type < -2.5)
        {
          if(is.null(yaxis)) { yaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
          else
          {
            mindensity <- 0
            for(i in 1:m) { if(ptx[i] < mindensity) { mindensity <- ptx[i] } }
            vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(mindensity,ptmax))
          }
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2Dpt")
          legendpos <- list(x=1.0,y=1.0,xanchor="right",tracegroupgap=0)
          ptxmmm <- RcppOUPAPassageTimeDensity(c(tmode,tmedian,tmean),k,s,x,omega,rho,mu,sigma)
          ptxmode <- ptxmmm[1]
          ptxmedian <- ptxmmm[2]
          ptxmean <- ptxmmm[3]
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter",x=c(tmean,tmean),y=c(0,ptxmean),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,hoverinfo="x+y",legendgroup="mean") %>%
              add_trace(.,type="scatter",x=c(s,tmean),y=c(ptxmean,ptxmean),mode="lines",line=meandotline,hoverinfo="x+y",legendgroup="mean",showlegend=FALSE)
          }
          fig <-  add_trace(fig,type="scatter",x=c(tmedian,tmedian),y=c(0,ptxmedian),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,hoverinfo="x+y",legendgroup="median") %>%
            add_trace(.,type="scatter",x=c(s,tmedian),y=c(ptxmedian,ptxmedian),mode="lines",line=mediandotline,hoverinfo="x+y",legendgroup="median",showlegend=FALSE) %>%
            add_trace(.,type="scatter",x=c(tmode,tmode),y=c(0,ptxmode),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,hoverinfo="x+y",legendgroup="mode") %>%
            add_trace(.,type="scatter",x=c(s,tmode),y=c(ptxmode,ptxmode),mode="lines",line=modedotline,hoverinfo="x+y",legendgroup="mode",showlegend=FALSE) %>%
            add_trace(.,type="scatter",x=t,y=ptx,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=densityline,hoverinfo="x+y")
        }
        # OUP_A_PassageTimeModeMedianMean2DPt
        else
        {
          if(is.null(yaxis)) { yaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2DPt")
          legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
          Ptxmmm <- RcppOUPAPassageTimeProbability(c(tmode,tmedian,tmean),k,s,x,omega,rho,mu,sigma)
          Ptxmode <- Ptxmmm[1]
          Ptxmedian <- Ptxmmm[2]
          Ptxmean <- Ptxmmm[3]
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter",x=c(tmean,tmean),y=c(0,Ptxmean),name="<i>t</i><sub>mean</sub>",mode="lines",line=meandashline,hoverinfo="x+y",legendgroup="mean") %>%
              add_trace(.,type="scatter",x=c(s,tmean),y=c(Ptxmean,Ptxmean),mode="lines",line=meandotline,hoverinfo="x+y",legendgroup="mean",showlegend=FALSE)
          }
          fig <-  add_trace(fig,type="scatter",x=c(tmedian,tmedian),y=c(0,Ptxmedian),name="<i>t</i><sub>median</sub>",mode="lines",line=mediandashline,hoverinfo="x+y",legendgroup="median") %>%
            add_trace(.,type="scatter",x=c(s,tmedian),y=c(Ptxmedian,Ptxmedian),mode="lines",line=mediandotline,hoverinfo="x+y",legendgroup="median",showlegend=FALSE) %>%
            add_trace(.,type="scatter",x=c(tmode,tmode),y=c(0,Ptxmode),name="<i>t</i><sub>mode</sub>",mode="lines",line=modedashline,hoverinfo="x+y",legendgroup="mode") %>%
            add_trace(.,type="scatter",x=c(s,tmode),y=c(Ptxmode,Ptxmode),mode="lines",line=modedotline,hoverinfo="x+y",legendgroup="mode",showlegend=FALSE) %>%
            add_trace(.,type="scatter",x=t,y=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=probabilityline,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 2D continued
      else if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>z</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>z</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        minm <- s
        maxm <- s+1
        for(j in 1:n)
        {
          if(is.finite(tmeans[j]))
          {
            if(tmeans[j] > maxm) { maxm <- tmeans[j] }
            if(tmeans[j] < minm) { minm <- tmeans[j] }
          }
          if(is.finite(tmedians[j]))
          {
            if(tmedians[j] > maxm) { maxm <- tmedians[j] }
            if(tmedians[j] < minm) { minm <- tmedians[j] }
          }
          if(is.finite(tmodes[j]))
          {
            if(tmodes[j] > maxm) { maxm <- tmodes[j] }
            if(tmodes[j] < minm) { minm <- tmodes[j] }
          }
        }
        maxm <- 1.2*(maxm-minm)
        maxmscale <- 1
        while(maxm > maxmscale) { maxmscale <- 10*maxmscale }
        maxm <- round(maxm/maxmscale,2)*maxmscale+s
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(z[1],z[n]),zeroline=FALSE)
        if(x > k) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(minm,maxm),zeroline=FALSE,side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(minm,maxm),zeroline=FALSE) }
        meanline <- list(color=cyn$d,width=4,shape="spline",smoothing=1.3)
        medianline <- list(color=grn$d,width=4,shape="spline",smoothing=1.3)
        modeline <- list(color=blu$d,width=4,shape="spline",smoothing=1.3)
        horzline <- list(color=gry$d,width=1)
        fig <- plot_ly()
        # OUP_A_PassageTimeModeMedianMean2D
        if(type < -0.5)
        {
          legendx <- (k-z[1])/(z[n]-z[1])
          if(legendx < 0.15) { legendx <- 0.15 }
          else if(legendx > 0.85) { legendx <- 0.85 }
          legendpos <- list(x=legendx,y=1.0,xanchor="center",yanchor="top")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2D")
          if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=z,y=tmeans,name="<i>t</i><sub>mean</sub>(<i>z</i>)",mode="lines",line=meanline,hoverinfo="x+y") }
          fig <-  add_trace(fig,type="scatter",x=z,y=tmedians,name="<i>t</i><sub>median</sub>(<i>z</i>)",mode="lines",line=medianline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tmodes,name="<i>t</i><sub>mode</sub>(<i>z</i>)",mode="lines",line=modeline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(z[1],z[n]),y=c(s,s),mode="lines",line=horzline,hoverinfo="x+y",showlegend=FALSE) %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimeModeMedianMean2Dx
        else
        {
          meandashline <- list(color=cyn$d,dash="longdash",width=3)
          meandotline <- list(color=cyn$d,dash="dot",width=1)
          mediandashline <- list(color=grn$d,dash="dash",width=3)
          mediandotline <- list(color=grn$d,dash="dot",width=1)
          modedashline <- list(color=blu$d,dash="dot",width=3)
          modedotline <- list(color=blu$d,dash="dot",width=1)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2Dx")
          legendy <- minm
          if(is.finite(tmean)) { if(tmean > legendy) { legendy <- tmean } }
          if(is.finite(tmedian)) { if(tmedian > legendy) { legendy <- tmedian } }
          if(is.finite(tmode)) { if(tmode > legendy) { legendy <- tmode } }
          if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=c(x,x),y=c(tmedian,tmean),name="mean",mode="lines",line=meandotline,hoverinfo="x+y") }
          fig <- add_trace(fig,type="scatter",x=c(x,x),y=c(tmode,tmedian),name="median",mode="lines",line=mediandotline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(x,x),y=c(s,tmode),name="mode",mode="lines",line=modedotline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(z[1],z[n]),y=c(s,s),mode="lines",line=horzline,hoverinfo="x+y")
          if(x > k)
          {
            if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=c(x,z[n]),y=c(tmean,tmean),name="mean",mode="lines",line=meandashline,hoverinfo="x+y") }
            fig <- add_trace(fig,type="scatter",x=c(x,z[n]),y=c(tmedian,tmedian),name="median",mode="lines",line=mediandashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(x,z[n]),y=c(tmode,tmode),name="mode",mode="lines",line=modedashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>x</i>",bsym,"=",esym,format(x,digits=4),"<br><i>t</i><sub>mean</sub>(<i>x</i>)",bsym,"=",esym,format(tmean,digits=4),"<br><i>t</i><sub>median</sub>(<i>x</i>)",bsym,"=",esym,format(tmedian,digits=4),"<br><i>t</i><sub>mode</sub>(<i>x</i>)",bsym,"=",esym,format(tmode,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
          }
          else
          {
            if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tmean,tmean),name="mean",mode="lines",line=meandashline,hoverinfo="x+y") }
            fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tmedian,tmedian),name="median",mode="lines",line=mediandashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(tmode,tmode),name="mode",mode="lines",line=modedashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>x</i>",bsym,"=",esym,format(x,digits=4),"<br><i>t</i><sub>mean</sub>(<i>x</i>)",bsym,"=",esym,format(tmean,digits=4),"<br><i>t</i><sub>median</sub>(<i>x</i>)",bsym,"=",esym,format(tmedian,digits=4),"<br><i>t</i><sub>mode</sub>(<i>x</i>)",bsym,"=",esym,format(tmode,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
          }
          if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=z,y=tmeans,name="mean(<i>z</i>)",mode="lines",line=meanline,hoverinfo="x+y") }
          fig <- add_trace(fig,type="scatter",x=z,y=tmedians,name="median(<i>z</i>)",mode="lines",line=medianline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tmodes,name="mode(<i>z</i>)",mode="lines",line=modeline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,annotations=xmmm,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
        # OUP_A_PassageTimeModeMedianMean3DDensity
        if(type < 1.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          spy <- list(x=2.35,y=-0.85,z=0.5)
          mindensity <- 0
          ptmeans <- vector("double",n)
          ptmedians <- vector("double",n)
          ptmodes <- vector("double",n)
          coordinatemeans <- vector("character",n)
          coordinatemedians <- vector("character",n)
          coordinatemodes <- vector("character",n)
          coordinateptmeans <- vector("character",n)
          coordinateptmedians <- vector("character",n)
          coordinateptmodes <- vector("character",n)
          coordinatex <- vector("character",m)
          coordinates <- matrix("",m,n)
          xx <- vector("double",m)
          ptxmmm <- RcppOUPAPassageTimeDensity(c(tmode,tmedian,tmean),k,s,x,omega,rho,mu,sigma)
          ptxmode <- ptxmmm[1]
          ptxmedian <- ptxmmm[2]
          ptxmean <- ptxmmm[3]
          for(j in 1:n)
          {
            ptmmm <- RcppOUPAPassageTimeDensity(c(tmodes[j],tmedians[j],tmeans[j]),k,s,z[j],omega,rho,mu,sigma)
            ptmodes[j] <- ptmmm[1]
            ptmedians[j] <- ptmmm[2]
            ptmeans[j] <- ptmmm[3]
            if(rho > 0)
            {
              coordinatemeans[j] <- paste(sep="","mean=",format(tmeans[j],digits=4),"<br><i>x</i>=",z[j])
              coordinateptmeans[j] <- paste(sep="","<i>p<sub>t</sub></i>(mean)=",format(ptmeans[j],digits=4),"<br><i>t</i>=",tmeans[j],"<br><i>x</i>=",z[j])
              if(ptmeans[j] < mindensity) { mindensity <- ptmeans[j] }
            }
            coordinatemedians[j] <- paste(sep="","median=",format(tmedians[j],digits=4),"<br><i>x</i>=",z[j])
            coordinateptmedians[j] <- paste(sep="","<i>p<sub>t</sub></i>(median)=",format(ptmedians[j],digits=4),"<br><i>t</i>=",tmedians[j],"<br><i>x</i>=",z[j])
            if(ptmedians[j] < mindensity) { mindensity <- ptmedians[j] }
            coordinatemodes[j] <- paste(sep="","mode=",format(tmodes[j],digits=4),"<br><i>x</i>=",z[j])
            coordinateptmodes[j] <- paste(sep="","<i>p<sub>t</sub></i>(mode)=",format(ptmodes[j],digits=4),"<br><i>t</i>=",tmodes[j],"<br><i>x</i>=",z[j])
            if(ptmodes[j] < mindensity) { mindensity <- ptmodes[j] }
          }
          for(i in 1:m)
          {
            if(ptx[i] < mindensity) { mindensity <- ptx[i] }
            coordinatex[i] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
            for(j in 1:n)
            {
              if(pt[i,j] < mindensity) { mindensity <- pt[i,j] }
              coordinates[i,j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j])
            }
            xx[i] <- x
          }
          ptxmesh <- MeshCurtainSmooth(xx,t,ptx,rep(0,m))
          ptmeshmeans <- MeshCurtainSmooth(z,tmeans,ptmeans,rep(0,m))
          ptmeshmedians <- MeshCurtainSmooth(z,tmedians,ptmedians,rep(0,m))
          ptmeshmodes <- MeshCurtainSmooth(z,tmodes,ptmodes,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(1.03*mindensity-0.03*ptmax,1.03*ptmax-0.03*mindensity),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          meanline <- list(color=cyn$e,width=8)
          medianline <- list(color=grn$e,width=8)
          modeline <- list(color=blu$e,width=8)
          meandashline <- list(color=blu$e,dash="longdash",width=8)
          mediandashline <- list(color=blu$e,dash="longdash",width=8)
          modedashline <- list(color=blu$e,dash="longdash",width=8)
          ptmeanline <- list(color=cyn$d,width=8)
          ptmedianline <- list(color=grn$d,width=8)
          ptmodeline <- list(color=blu$d,width=8)
          ptxline <- list(color=gry$e,width=8)
          ptline <- list(color=blu$e,width=8)
          gradientptx <- list(c(0,blu$a),c(1,blu$a))
          gradient <- list(c(0,blu$c),c(1,blu$c))
          gradientmean <- list(c(0,cyn$d),c(1,cyn$d))
          gradientmedian <- list(c(0,grn$c),c(1,grn$c))
          gradientmode <- list(c(0,blu$b),c(1,blu$b))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeModeMedianMean3DDensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,legendgroup="ptx") %>%
            add_trace(.,type="mesh3d",x=ptxmesh$xvertex,y=ptxmesh$yvertex,z=ptxmesh$zvertex,i=ptxmesh$ivertex,j=ptxmesh$jvertex,k=ptxmesh$kvertex,intensity=ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientptx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE)
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=c(x,x),y=c(tmean,tmean),z=c(0,ptxmean),mode="lines",line=meandashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) }
          fig <- add_trace(fig,type="scatter3d",x=c(x,x),y=c(tmedian,tmedian),z=c(0,ptxmedian),mode="lines",line=mediandashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tmode,tmode),z=c(0,ptxmode),mode="lines",line=modedashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE)
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=rep(0,n),name="<i>t</i><sub>mean</sub>(<i>z</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=rep(0,n),name="<i>t</i><sub>medn</sub>(<i>z</i>)",mode="lines",line=medianline,hoverinfo="text",text=coordinatemedians) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=rep(0,n),name="<i>t</i><sub>mode</sub>(<i>z</i>)",mode="lines",line=modeline,hoverinfo="text",text=coordinatemodes)
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=ptmeans,name="<i>p<sub>t</sub></i>(mean)",mode="lines",line=ptmeanline,hoverinfo="text",text=coordinateptmeans,legendgroup="ptmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmeans$xvertex,y=ptmeshmeans$yvertex,z=ptmeshmeans$zvertex,i=ptmeshmeans$ivertex,j=ptmeshmeans$jvertex,k=ptmeshmeans$kvertex,intensity=ptmeshmeans$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="ptmeans",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=ptmedians,name="<i>p<sub>t</sub></i>(medn)",mode="lines",line=ptmedianline,hoverinfo="text",text=coordinateptmedians,legendgroup="ptmedians",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmedians$xvertex,y=ptmeshmedians$yvertex,z=ptmeshmedians$zvertex,i=ptmeshmedians$ivertex,j=ptmeshmedians$jvertex,k=ptmeshmedians$kvertex,intensity=ptmeshmedians$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmedian,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="ptmedians",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=ptmodes,name="<i>p<sub>t</sub></i>(mode)",mode="lines",line=ptmodeline,hoverinfo="text",text=coordinateptmodes,legendgroup="ptmodes",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmodes$xvertex,y=ptmeshmodes$yvertex,z=ptmeshmodes$zvertex,i=ptmeshmodes$ivertex,j=ptmeshmodes$jvertex,k=ptmeshmodes$kvertex,intensity=ptmeshmodes$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmode,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="ptmodes",visible="legendonly",showlegend=FALSE)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 1
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly")
          while(j < n)
          {
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity-0.07 }
            else { lineopacity <- lineopacity+0.07 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="surface",x=z,y=t,z=pt,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        # OUP_A_PassageTimeModeMedianMean3DProbability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          spy <- list(x=2.35,y=-0.85,z=0.5)
          Ptmeans <- vector("double",n)
          Ptmedians <- vector("double",n)
          Ptmodes <- vector("double",n)
          coordinatemeans <- vector("character",n)
          coordinatemedians <- vector("character",n)
          coordinatemodes <- vector("character",n)
          coordinatePtmeans <- vector("character",n)
          coordinatePtmedians <- vector("character",n)
          coordinatePtmodes <- vector("character",n)
          coordinatex <- vector("character",m)
          coordinates <- matrix("",m,n)
          xx <- vector("double",m)
          kindex <- 0
          Ptxmmm <- RcppOUPAPassageTimeProbability(c(tmode,tmedian,tmean),k,s,x,omega,rho,mu,sigma)
          Ptxmode <- Ptxmmm[1]
          Ptxmedian <- Ptxmmm[2]
          Ptxmean <- Ptxmmm[3]
          for(j in 1:n)
          {
            Ptmmm <- RcppOUPAPassageTimeProbability(c(tmodes[j],tmedians[j],tmeans[j]),k,s,z[j],omega,rho,mu,sigma)
            Ptmodes[j] <- Ptmmm[1]
            Ptmedians[j] <- Ptmmm[2]
            Ptmeans[j] <- Ptmmm[3]
            if(rho > 0)
            {
              coordinatemeans[j] <- paste(sep="","mean=",format(tmeans[j],digits=4),"<br><i>x</i>=",z[j])
              coordinatePtmeans[j] <- paste(sep="","<i>P<sub>t</sub></i>(mean)=",format(Ptmeans[j],digits=4),"<br><i>t</i>=",tmeans[j],"<br><i>x</i>=",z[j])
            }
            coordinatemedians[j] <- paste(sep="","median=",format(tmedians[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatePtmedians[j] <- paste(sep="","<i>P<sub>t</sub></i>(median)=",format(Ptmedians[j],digits=4),"<br><i>t</i>=",tmedians[j],"<br><i>x</i>=",z[j])
            coordinatemodes[j] <- paste(sep="","mode=",format(tmodes[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatePtmodes[j] <- paste(sep="","<i>P<sub>t</sub></i>(mode)=",format(Ptmodes[j],digits=4),"<br><i>t</i>=",tmodes[j],"<br><i>x</i>=",z[j])
            if(z[j] == k) { kindex <- j }
          }
          for(i in 1:m)
          {
            coordinatex[i] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(Ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
            for(j in 1:n)
            {
              coordinates[i,j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(Pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j])
            }
            xx[i] <- x
          }
          if(kindex > 1 && kindex < n && sigma > 0)
          {
            Ptmeans[kindex] <- 0.5*(Ptmeans[kindex-1]+Ptmeans[kindex+1])
            Ptmedians[kindex] <- 0.5*(Ptmedians[kindex-1]+Ptmedians[kindex+1])
            Ptmodes[kindex] <- 0.5*(Ptmodes[kindex-1]+Ptmodes[kindex+1])
          }
          Ptxmesh <- MeshCurtainSmooth(xx,t,Ptx,rep(0,m))
          Ptmeshmeans <- MeshCurtainSmooth(z,tmeans,Ptmeans,rep(0,m))
          Ptmeshmedians <- MeshCurtainSmooth(z,tmedians,Ptmedians,rep(0,m))
          Ptmeshmodes <- MeshCurtainSmooth(z,tmodes,Ptmodes,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          meanline <- list(color=cyn$e,width=8)
          medianline <- list(color=grn$e,width=8)
          modeline <- list(color=blu$e,width=8)
          meandashline <- list(color=cyn$e,dash="longdash",width=8)
          mediandashline <- list(color=grn$e,dash="longdash",width=8)
          modedashline <- list(color=blu$e,dash="longdash",width=8)
          Ptmeanline <- list(color=cyn$d,width=8)
          Ptmedianline <- list(color=grn$d,width=8)
          Ptmodeline <- list(color=blu$d,width=8)
          Ptxline <- list(color=gry$e,width=8)
          Ptline <- list(color=grn$e,width=8)
          gradientPtx <- list(c(0,grn$a),c(1,grn$a))
          gradient <- list(c(0,grn$c),c(1,grn$c))
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientmedian <- list(c(0,grn$b),c(1,grn$b))
          gradientmode <- list(c(0,blu$b),c(1,blu$b))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeModeMedianMean3DProbability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,legendgroup="Ptx") %>%
            add_trace(.,type="mesh3d",x=Ptxmesh$xvertex,y=Ptxmesh$yvertex,z=Ptxmesh$zvertex,i=Ptxmesh$ivertex,j=Ptxmesh$jvertex,k=Ptxmesh$kvertex,intensity=Ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientPtx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE)
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=c(x,x),y=c(tmean,tmean),z=c(0,Ptxmean),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=meandashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) }
          fig <- add_trace(fig,type="scatter3d",x=c(x,x),y=c(tmedian,tmedian),z=c(0,Ptxmedian),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=mediandashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tmode,tmode),z=c(0,Ptxmode),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=modedashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE)
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=rep(0,n),name="<i>t</i><sub>mean</sub>(<i>z</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=rep(0,n),name="<i>t</i><sub>medn</sub>(<i>z</i>)",mode="lines",line=medianline,hoverinfo="text",text=coordinatemedians) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=rep(0,n),name="<i>t</i><sub>mode</sub>(<i>z</i>)",mode="lines",line=modeline,hoverinfo="text",text=coordinatemodes)
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=Ptmeans,name="<i>P<sub>t</sub></i>(mean)",mode="lines",line=Ptmeanline,hoverinfo="text",text=coordinatePtmeans,legendgroup="Ptmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmeans$xvertex,y=Ptmeshmeans$yvertex,z=Ptmeshmeans$zvertex,i=Ptmeshmeans$ivertex,j=Ptmeshmeans$jvertex,k=Ptmeshmeans$kvertex,intensity=Ptmeshmeans$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Ptmeans",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=Ptmedians,name="<i>P<sub>t</sub></i>(medn)",mode="lines",line=Ptmedianline,hoverinfo="text",text=coordinatePtmedians,legendgroup="Ptmedians",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmedians$xvertex,y=Ptmeshmedians$yvertex,z=Ptmeshmedians$zvertex,i=Ptmeshmedians$ivertex,j=Ptmeshmedians$jvertex,k=Ptmeshmedians$kvertex,intensity=Ptmeshmedians$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmedian,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Ptmedians",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=Ptmodes,name="<i>P<sub>t</sub></i>(mode)",mode="lines",line=Ptmodeline,hoverinfo="text",text=coordinatePtmodes,legendgroup="Ptmodes",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmodes$xvertex,y=Ptmeshmodes$yvertex,z=Ptmeshmodes$zvertex,i=Ptmeshmodes$ivertex,j=Ptmeshmodes$jvertex,k=Ptmeshmodes$kvertex,intensity=Ptmeshmodes$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmode,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="Ptmodes",visible="legendonly",showlegend=FALSE)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 1
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly")
          while(j < n)
          {
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity-0.07 }
            else { lineopacity <- lineopacity+0.07}
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly",showlegend=FALSE)
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Ptz",visible="legendonly",showlegend=FALSE) }
          fig <-add_trace(fig,type="surface",x=z,y=t,z=Pt,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time lower, median and upper percentiles
    #' @param type  = -3,...,2 or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum scale for vertical axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotPassageTimePercentiles = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,5)[[1]]
      self$set_plot_args(NULL,ptmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      Ppct <- private$t_stoch_args[[7]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      if(Ppct > 0.5)
      {
        Pupper <- Ppct
        Phalf <- 0.5
        Plower <- 1-Ppct
      }
      else
      {
        Pupper <- 1-Ppct
        Phalf <- 0.5
        Plower <- Ppct
      }
      m <- length(t)
      n <- length(z)
      tpercentile <- private$tpercentile
      tpercentiles <- private$tpercentiles
      if(is.null(tpercentile) || is.null(tpercentiles))
      {
        percentiles <- self$PassageTimePercentiles(who="A")
        tpercentile <- percentiles[[1]]
        tpercentiles <- percentiles[[2]]
      }
      tlower <- tpercentile[[1]]
      thalf <- tpercentile[[2]]
      tupper <- tpercentile[[3]]
      tlowers <- tpercentiles[[1]]
      thalfs <- tpercentiles[[2]]
      tuppers <- tpercentiles[[3]]
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) || is.null(pt))
      {
        densities <- self$PassageTimeDensity(who="A")
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) || is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(who="A")
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,,drop=FALSE]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        tlowers <- tlowers[Ixbeg:Ixend]
        thalfs <- thalfs[Ixbeg:Ixend]
        tuppers <- tuppers[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend,drop=FALSE]
        Pt <- Pt[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        if(Ppct > 0.5) { clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Percentiles",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("Ppct",Ppct,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("t",1-Ppct,"(x)"),tpercentile[[1]],paste0("t",1-Ppct,"(z)"),tpercentiles[[1]]),c(paste0("t0.5(x)"),tpercentile[[2]],paste0("t0.5(z)"),tpercentiles[[2]]),c(paste0("t",Ppct,"(x)"),tpercentile[[3]],paste0("t",Ppct,"(z)"),tpercentiles[[3]])) }
        else { clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Percentiles",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("Ppct",Ppct,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("x",x,"z",z),c(paste0("t",Ppct,"(x)"),tpercentile[[1]],paste0("t",Ppct,"(z)"),tpercentiles[[1]]),c(paste0("t0.5(x)"),tpercentile[[2]],paste0("t0.5(z)"),tpercentiles[[2]]),c(paste0("t",1-Ppct,"(x)"),tpercentile[[3]],paste0("t",1-Ppct,"(z)"),tpercentiles[[3]])) }
        private$CopyToClipboard(clip)
      }
      # plot ----
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>w</i>=",esym,format(omega,digits=4),",<i>P</i><sub>%</sub>",bsym,"=",esym,format(Ppct,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Passage Time Percentiles"  }
      }
      else if(is.null(title)) { title <- ""  }
      lookup <- list(text=title,yref="container",y=0.95)
      # 2D
      if(type < -1.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        densityline <- list(color=blu$d,width=4,shape="spline",smoothing=1.3)
        probabilityline <- list(color=grn$d,width=4,shape="spline",smoothing=1.3)
        upperdashline <- list(color=grn$d,dash="longdash",width=3)
        upperdotline <- list(color=grn$d,dash="dot",width=1)
        halfdashline <- list(color=grn$d,dash="dash",width=3)
        halfdotline <- list(color=grn$d,dash="dot",width=1)
        lowerdashline <- list(color=grn$d,dash="dot",width=3)
        lowerdotline <- list(color=grn$d,dash="dot",width=1)
        fig <- plot_ly()
        # OUP_A_PassageTimePercentiles2Dpt
        if(type < -2.5)
        {
          if(is.null(yaxis)) { yaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
          else
          {
            mindensity <- 0
            for(i in 1:m) { if(ptx[i] < mindensity) { mindensity <- ptx[i] } }
            vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(mindensity,ptmax))
          }
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimePercentiles2Dpt")
          legendpos <- list(x=1.0,y=1.0,xanchor="right",tracegroupgap=0)
          ptxlhu <- RcppOUPAPassageTimeDensity(c(tlower,thalf,tupper),k,s,x,omega,rho,mu,sigma)
          ptxlower <- ptxlhu[1]
          ptxhalf <- ptxlhu[2]
          ptxupper <- ptxlhu[3]
          fig <- add_trace(fig,type="scatter",x=c(tupper,tupper),y=c(0,ptxupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tupper),y=c(ptxupper,ptxupper),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(thalf,thalf),y=c(0,ptxhalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,thalf),y=c(ptxhalf,ptxhalf),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tlower,tlower),y=c(0,ptxlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tlower),y=c(ptxlower,ptxlower),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=ptx,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=densityline,hoverinfo="x+y")
        }
        # OUP_A_PassageTimePercentiles2DPt
        else
        {
          if(is.null(yaxis)) { yaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimePercentiles2DPt")
          legendpos <- list(x=1.0,y=0.2,xanchor="right",tracegroupgap=0)
          Ptxlhu <- RcppOUPAPassageTimeProbability(c(tlower,thalf,tupper),k,s,x,omega,rho,mu,sigma)
          Ptxlower <- Ptxlhu[1]
          Ptxhalf <- Ptxlhu[2]
          Ptxupper <- Ptxlhu[3]
          fig <- add_trace(fig,type="scatter",x=c(tupper,tupper),y=c(0,Ptxupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tupper),y=c(Ptxupper,Ptxupper),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(thalf,thalf),y=c(0,Ptxhalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,thalf),y=c(Ptxhalf,Ptxhalf),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tlower,tlower),y=c(0,Ptxlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tlower),y=c(Ptxlower,Ptxlower),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=probabilityline,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 2D continued
      else if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>z</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        if(is.null(xaxis)) { xaxis <- "<i>z</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        maxpct <- s+1
        for(j in 1:n) { if(tuppers[j] > maxpct) { maxpct <- tuppers[j] } }
        maxpct <- 1.2*(maxpct-s)
        maxpctscale <- 1
        while(maxpct > maxpctscale) { maxpctscale <- 10*maxpctscale }
        maxpct <- round(maxpct/maxpctscale,2)*maxpctscale+s
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(z[1],z[n]),zeroline=FALSE)
        if(x > k) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(s,maxpct),zeroline=FALSE,side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(s,maxpct),zeroline=FALSE) }
        upperline <- list(color=grn$c,width=4,shape="spline",smoothing=1.3)
        halfline <- list(color=grn$d,width=4,shape="spline",smoothing=1.3)
        lowerline <- list(color=grn$c,width=4,shape="spline",smoothing=1.3)
        horzline <- list(color=gry$d,width=1)
        fig <- plot_ly()
        # OUP_A_PassageTimePercentiles2D
        if(type < -0.5)
        {
          legendx <- (k-z[1])/(z[n]-z[1])
          if(legendx < 0.15) { legendx <- 0.15 }
          else if(legendx > 0.85) { legendx <- 0.85 }
          legendpos <- list(x=legendx,y=1.0,xanchor="center",yanchor="top")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimePercentiles2D")
          fig <- add_trace(fig,type="scatter",x=z,y=tuppers,name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=upperline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=thalfs,name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=halfline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tlowers,name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=lowerline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(z[1],z[n]),y=c(s,s),mode="lines",line=horzline,showlegend=FALSE) %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimePercentiles2Dx
        else
        {
          upperdashline <- list(color=grn$d,dash="longdash",width=3)
          upperdotline <- list(color=grn$d,dash="dot",width=1)
          halfdashline <- list(color=grn$d,dash="dash",width=3)
          halfdotline <- list(color=grn$d,dash="dot",width=1)
          lowerdashline <- list(color=grn$d,dash="dot",width=3)
          lowerdotline <- list(color=grn$d,dash="dot",width=1)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimePercentiles2Dx")
          legendy <- s
          if(tupper > legendy) { legendy <- tupper }
          fig <- add_trace(fig,type="scatter",x=c(x,x),y=c(thalf,tupper),name="<i>x</i>",mode="lines",line=upperdotline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(x,x),y=c(tlower,thalf),name="<i>x</i>",mode="lines",line=halfdotline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(x,x),y=c(s,tlower),name="<i>x</i>",mode="lines",line=lowerdotline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(z[1],z[n]),y=c(s,s),mode="lines",line=horzline,hoverinfo="x+y")
          if(x > k)
          {
            fig <- add_trace(fig,type="scatter",x=c(x,z[n]),y=c(tupper,tupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(x,z[n]),y=c(thalf,thalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(x,z[n]),y=c(tlower,tlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>x</i>",bsym,"=",esym,format(x,digits=4),"<br><i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(tupper,digits=4),"<br><i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(thalf,digits=4),"<br><i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(tlower,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
          }
          else
          {
            fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tupper,tupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(thalf,thalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(tlower,tlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>x</i>",bsym,"=",esym,format(x,digits=4),"<br><i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(tupper,digits=4),"<br><i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(thalf,digits=4),"<br><i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>x</i>)",bsym,"=",esym,format(tlower,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
          }
          fig <- add_trace(fig,type="scatter",x=z,y=tuppers,name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=thalfs,name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tlowers,name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
            layout(.,title=lookup,annotations=xmmm,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
        # OUP_A_PassageTimePercentiles3DDensity
        if(type < 1.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          spy <- list(x=2.35,y=-0.85,z=0.5)
          mindensity <- 0
          ptuppers <- vector("double",n)
          pthalfs <- vector("double",n)
          ptlowers <- vector("double",n)
          coordinateuppers <- vector("character",n)
          coordinatehalfs <- vector("character",n)
          coordinatelowers <- vector("character",n)
          coordinateptuppers <- vector("character",n)
          coordinatepthalfs <- vector("character",n)
          coordinateptlowers <- vector("character",n)
          coordinatex <- vector("character",m)
          coordinates <- matrix("",m,n)
          xx <- vector("double",m)
          ptxlhu <- RcppOUPAPassageTimeDensity(c(tlower,thalf,tupper),k,s,x,omega,rho,mu,sigma)
          ptxlower <- ptxlhu[1]
          ptxhalf <- ptxlhu[2]
          ptxupper <- ptxlhu[3]
          for(j in 1:n)
          {
            ptlhu <- RcppOUPAPassageTimeDensity(c(tlowers[j],thalfs[j],tuppers[j]),k,s,z[j],omega,rho,mu,sigma)
            ptlowers[j] <- ptlhu[1]
            pthalfs[j] <- ptlhu[2]
            ptuppers[j] <- ptlhu[3]
            coordinateuppers[j] <- paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>=",format(tuppers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatehalfs[j] <- paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>=",format(thalfs[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatelowers[j] <- paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>=",format(tlowers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinateptuppers[j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)=",format(ptuppers[j],digits=4),"<br><i>t</i>=",format(tuppers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatepthalfs[j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)=",format(pthalfs[j],digits=4),"<br><i>t</i>=",format(thalfs[j],digits=4),"<br><i>x</i>=",z[j])
            coordinateptlowers[j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)=",format(ptlowers[j],digits=4),"<br><i>t</i>=",format(tlowers[j],digits=4),"<br><i>x</i>=",z[j])
          }
          for(i in 1:m)
          {
            if(ptx[i] < mindensity) { mindensity <- ptx[i] }
            coordinatex[i] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
            for(j in 1:n)
            {
              if(pt[i,j] < mindensity) { mindensity <- pt[i,j] }
              coordinates[i,j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j])
            }
            xx[i] <- x
          }
          ptxmesh <- MeshCurtainSmooth(xx,t,ptx,rep(0,m))
          ptmeshuppers <- MeshCurtainSmooth(z,tuppers,ptuppers,rep(0,m))
          ptmeshhalfs <- MeshCurtainSmooth(z,thalfs,pthalfs,rep(0,m))
          ptmeshlowers <- MeshCurtainSmooth(z,tlowers,ptlowers,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(1.03*mindensity-0.03*ptmax,1.03*ptmax-0.03*mindensity),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          upperline <- list(color=blu$e,width=8)
          halfline <- list(color=blu$e,width=8)
          lowerline <- list(color=blu$e,width=8)
          upperdashline <- list(color=blu$e,dash="longdash",width=8)
          halfdashline <- list(color=blu$e,dash="longdash",width=8)
          lowerdashline <- list(color=blu$e,dash="longdash",width=8)
          ptupperline <- list(color=blu$d,width=8)
          pthalfline <- list(color=blu$d,width=8)
          ptlowerline <- list(color=blu$d,width=8)
          ptxline <- list(color=gry$e,width=8)
          ptline <- list(color=blu$e,width=8)
          gradientptx <- list(c(0,blu$a),c(1,blu$a))
          gradient <- list(c(0,blu$c),c(1,blu$c))
          gradientupper <- list(c(0,blu$d),c(1,blu$d))
          gradienthalf <- list(c(0,blu$c),c(1,blu$c))
          gradientlower <- list(c(0,blu$b),c(1,blu$b))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimePercentiles3DDensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,legendgroup="ptx") %>%
            add_trace(.,type="mesh3d",x=ptxmesh$xvertex,y=ptxmesh$yvertex,z=ptxmesh$zvertex,i=ptxmesh$ivertex,j=ptxmesh$jvertex,k=ptxmesh$kvertex,intensity=ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientptx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tupper,tupper),z=c(0,ptxupper),name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=upperdashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(thalf,thalf),z=c(0,ptxhalf),name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=halfdashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tlower,tlower),z=c(0,ptxlower),name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=lowerdashline,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=upperline,hoverinfo="text",text=coordinateuppers) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=halfline,hoverinfo="text",text=coordinatehalfs) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=lowerline,hoverinfo="text",text=coordinatelowers) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=ptuppers,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)"),mode="lines",line=ptupperline,hoverinfo="text",text=coordinateptuppers,legendgroup="ptuppers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshuppers$xvertex,y=ptmeshuppers$yvertex,z=ptmeshuppers$zvertex,i=ptmeshuppers$ivertex,j=ptmeshuppers$jvertex,k=ptmeshuppers$kvertex,intensity=ptmeshuppers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientupper,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="ptuppers",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=pthalfs,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)"),mode="lines",line=pthalfline,hoverinfo="text",text=coordinatepthalfs,legendgroup="pthalfs",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshhalfs$xvertex,y=ptmeshhalfs$yvertex,z=ptmeshhalfs$zvertex,i=ptmeshhalfs$ivertex,j=ptmeshhalfs$jvertex,k=ptmeshhalfs$kvertex,intensity=ptmeshhalfs$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthalf,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="pthalfs",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=ptlowers,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)"),mode="lines",line=ptlowerline,hoverinfo="text",text=coordinateptlowers,legendgroup="ptlowers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshlowers$xvertex,y=ptmeshlowers$yvertex,z=ptmeshlowers$zvertex,i=ptmeshlowers$ivertex,j=ptmeshlowers$jvertex,k=ptmeshlowers$kvertex,intensity=ptmeshlowers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientlower,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="ptlowers",visible="legendonly",showlegend=FALSE)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 1
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly")
          while(j < n)
          {
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity-0.07 }
            else { lineopacity <- lineopacity+0.07 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="surface",x=z,y=t,z=pt,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        # OUP_A_PassageTimePercentiles3DProbability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          spy <- list(x=2.35,y=-0.85,z=0.5)
          Ptuppers <- vector("double",n)
          Pthalfs <- vector("double",n)
          Ptlowers <- vector("double",n)
          coordinateuppers <- vector("character",n)
          coordinatehalfs <- vector("character",n)
          coordinatelowers <- vector("character",n)
          coordinatePtuppers <- vector("character",n)
          coordinatePthalfs <- vector("character",n)
          coordinatePtlowers <- vector("character",n)
          coordinatex <- vector("character",m)
          coordinates <- matrix("",m,n)
          xx <- vector("double",m)
          kindex <- 0
          Ptxlhu <- RcppOUPAPassageTimeProbability(c(tlower,thalf,tupper),k,s,x,omega,rho,mu,sigma)
          Ptxlower <- Ptxlhu[1]
          Ptxhalf <- Ptxlhu[2]
          Ptxupper <- Ptxlhu[3]
          for(j in 1:n)
          {
            PInf <- private$OUPPassageTimeProbabilityInf(z[j],k,omega,rho,mu,sigma)
            Ptuppers[j] <- Pupper*PInf
            Pthalfs[j] <- Phalf*PInf
            Ptlowers[j] <- Plower*PInf
            coordinateuppers[j] <- paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>=",format(tuppers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatehalfs[j] <- paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>=",format(thalfs[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatelowers[j] <- paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>=",format(tlowers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatePtuppers[j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)=",format(Ptuppers[j],digits=4),"<br><i>t</i>=",format(tuppers[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatePthalfs[j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)=",format(Pthalfs[j],digits=4),"<br><i>t</i>=",format(thalfs[j],digits=4),"<br><i>x</i>=",z[j])
            coordinatePtlowers[j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)=",format(Ptlowers[j],digits=4),"<br><i>t</i>=",format(tlowers[j],digits=4),"<br><i>x</i>=",z[j])
            if(z[j] == k) { kindex <- j }
          }
          for(i in 1:m)
          {
            coordinatex[i] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(Ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
            for(j in 1:n)
            {
              coordinates[i,j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(Pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j])
            }
            xx[i] <- x
          }
          Ptxmesh <- MeshCurtainSmooth(xx,t,Ptx,rep(0,m))
          Ptmeshuppers <- MeshCurtainSmooth(z,tuppers,Ptuppers,rep(0,m))
          Ptmeshhalfs <- MeshCurtainSmooth(z,thalfs,Pthalfs,rep(0,m))
          Ptmeshlowers <- MeshCurtainSmooth(z,tlowers,Ptlowers,rep(0,m))
          xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          upperline <- list(color=grn$e,width=8)
          halfline <- list(color=grn$e,width=8)
          lowerline <- list(color=grn$e,width=8)
          upperdashline <- list(color=grn$e,dash="longdash",width=8)
          halfdashline <- list(color=grn$e,dash="longdash",width=8)
          lowerdashline <- list(color=grn$e,dash="longdash",width=8)
          Ptupperline <- list(color=grn$d,width=8)
          Pthalfline <- list(color=grn$d,width=8)
          Ptlowerline <- list(color=grn$d,width=8)
          Ptxline <- list(color=gry$e,width=8)
          Ptline <- list(color=grn$e,width=8)
          gradientPtx <- list(c(0,grn$a),c(1,grn$a))
          gradient <- list(c(0,grn$c),c(1,grn$c))
          gradientupper <- list(c(0,grn$b),c(1,grn$b))
          gradienthalf <- list(c(0,grn$c),c(1,grn$c))
          gradientlower <- list(c(0,grn$d),c(1,grn$d))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimePercentiles3DProbability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,legendgroup="Ptx") %>%
            add_trace(.,type="mesh3d",x=Ptxmesh$xvertex,y=Ptxmesh$yvertex,z=Ptxmesh$zvertex,i=Ptxmesh$ivertex,j=Ptxmesh$jvertex,k=Ptxmesh$kvertex,intensity=Ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientPtx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tupper,tupper),z=c(0,Ptxupper),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=upperdashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(thalf,thalf),z=c(0,Ptxhalf),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=halfdashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=c(x,x),y=c(tlower,tlower),z=c(0,Ptxlower),name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=lowerdashline,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=upperline,hoverinfo="text",text=coordinateuppers) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=halfline,hoverinfo="text",text=coordinatehalfs) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=lowerline,hoverinfo="text",text=coordinatelowers) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=Ptuppers,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)"),mode="lines",line=Ptupperline,hoverinfo="text",text=coordinatePtuppers,legendgroup="Ptuppers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshuppers$xvertex,y=Ptmeshuppers$yvertex,z=Ptmeshuppers$zvertex,i=Ptmeshuppers$ivertex,j=Ptmeshuppers$jvertex,k=Ptmeshuppers$kvertex,intensity=Ptmeshuppers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientupper,reversescale=reverse,opacity=0.3,hoverinfo="skip",legendgroup="Ptuppers",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=Pthalfs,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)"),mode="lines",line=Pthalfline,hoverinfo="text",text=coordinatePthalfs,legendgroup="Pthalfs",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshhalfs$xvertex,y=Ptmeshhalfs$yvertex,z=Ptmeshhalfs$zvertex,i=Ptmeshhalfs$ivertex,j=Ptmeshhalfs$jvertex,k=Ptmeshhalfs$kvertex,intensity=Ptmeshhalfs$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthalf,reversescale=reverse,opacity=0.3,hoverinfo="skip",legendgroup="Pthalfs",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=Ptlowers,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)"),mode="lines",line=Ptlowerline,hoverinfo="text",text=coordinatePtlowers,legendgroup="Ptlowers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshlowers$xvertex,y=Ptmeshlowers$yvertex,z=Ptmeshlowers$zvertex,i=Ptmeshlowers$ivertex,j=Ptmeshlowers$jvertex,k=Ptmeshlowers$kvertex,intensity=Ptmeshlowers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientlower,reversescale=reverse,opacity=0.3,hoverinfo="skip",legendgroup="Ptlowers",visible="legendonly",showlegend=FALSE)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 1
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly")
          while(j < n)
          {
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity-0.07 }
            else { lineopacity <- lineopacity+0.07 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly",showlegend=FALSE)
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Ptz",visible="legendonly",showlegend=FALSE) }
          fig <-add_trace(fig,type="surface",x=z,y=t,z=Pt,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time densities
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param ptmax   maximum scale for vertical axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotPassageTimeDensity = function(type=NULL,ptmax=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,6)[[1]]
      self$set_plot_args(NULL,ptmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      ptmax <- private$plot_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      m <- length(t)
      n <- length(z)
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) || is.null(pt))
      {
        densities <- self$PassageTimeDensity(who="A")
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Densities",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("t","ptx","pt(t,z)",z),cbind(t,ptx,t,pt))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>w</i>=",esym,format(omega,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Passage Time Densities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_PassageTimeDensity2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,rangemode="tozero",ticks="outside") }
        else
        {
          mindensity <- 0
          for(i in 1:m)
          {
            if(ptx[i] < mindensity) { mindensity <- ptx[i] }
            for(j in 1:n)
            {
              if(pt[i,j] < mindensity) { mindensity <- pt[i,j] }
            }
          }
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,range=c(mindensity,ptmax),ticks="outside")
        }
        ptxline <- list(color=blu$e,width=4)
        ptline <- list(color=blu$c,width=4)
        legendpos <- list(orientation="h",x=1.0,y=1.0,xanchor="right")
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeDensity2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=t,y=ptx,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=ptxline,hoverinfo="x+y")
        dx <- as.integer((n-1)/10)
        if(dx < 1) { dx <- 1 }
        lineopacity <- 1
        j <- 1
        q <- 1
        fig <- add_trace(fig,type="scatter",x=t,y=pt[,j],name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)"),mode="lines",line=ptline,opacity=lineopacity,hoverinfo="x+y",legendgroup="pt",visible="legendonly")
        while(j < n)
        {
          j <- j+dx
          q <- q+1
          if(q < 7) { lineopacity <- lineopacity-0.07 }
          else { lineopacity <- lineopacity+0.07 }
          fig <- add_trace(fig,type="scatter",x=t,y=pt[,j],name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|",z[j],")"),mode="lines",line=ptline,opacity=lineopacity,hoverinfo="x+y",legendgroup="pt",visible="legendonly",showlegend=FALSE)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_PassageTimeDensity3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        mindensity <- 0
        coordinatex <- vector("character",m)
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          coordinatex[i] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
          if(ptx[i] < mindensity) { mindensity <- ptx[i] }
          for(j in 1:n)
          {
            coordinates[i,j] <- paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j])
            if(pt[i,j] < mindensity) { mindensity <- pt[i,j] }
          }
        }
        spy <- list(x=2.35,y=-0.85,z=0.5)
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(1.03*mindensity-0.03*ptmax,1.03*ptmax-0.03*mindensity),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        xx <- vector("double",m)
        for(i in 1:m) { xx[i] <- x }
        ptxmesh <- MeshCurtainSmooth(xx,t,ptx,rep(0,m))
        ptxline <- list(color=font$color,width=8)
        ptline <- list(color=blu$d,width=8)
        gradientptx <- list(c(0,blu$a),c(1,blu$a))
        gradient <- list(c(0,blu$c),c(1,blu$c))
        if(k > mu) { rise <- list(x=10,y=100,z=0) }
        else if(k == mu) { rise <- list(x=0,y=100,z=0) }
        else { rise <- list(x=-10,y=100,z=0) }
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeDensity3D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,legendgroup="ptx") %>%
          add_trace(.,type="mesh3d",x=ptxmesh$xvertex,y=ptxmesh$yvertex,z=ptxmesh$zvertex,i=ptxmesh$ivertex,j=ptxmesh$jvertex,k=ptxmesh$kvertex,intensity=ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientptx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="ptx",showlegend=FALSE)
        dx <- as.integer((n-1)/10)
        if(dx < 1) { dx <- 1 }
        lineopacity <- 1
        j <- 1
        q <- 1
        for(i in 1:m) { xx[i] <- z[j] }
        fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly")
        while(j < n)
        {
          j <- j+dx
          q <- q+1
          if(q < 7) { lineopacity <- lineopacity-0.07 }
          else { lineopacity <- lineopacity+0.07 }
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="ptz",visible="legendonly",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=z,y=t,z=pt,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time probabilities
    #' @param type  = 0, 1 or 'n','p','d' for next, previous, default
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zaxis   text for z-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotPassageTimeProbability = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,6)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      m <- length(t)
      n <- length(z)
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) || is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(who="A")
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,,drop=FALSE]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        Pt <- Pt[,Ixbeg:Ixend,drop=FALSE]
        n <- length(z)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Analytical",rep("",n+2)),c("Passage Time Probabilities",rep("",n+2)),c("k",k,rep("",n+1)),c("s",s,rep("",n+1)),c("x",x,rep("",n+1)),c("k",k,rep("",n+1)),c("omega",omega,rep("",n+1)),c("rho",rho,rep("",n+1)),c("mu",mu,rep("",n+1)),c("sigma",sigma,rep("",n+1)),c("t","Ptx","Pt(t,z)",z),cbind(t,Ptx,t,Pt))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>k</i>",bsym,"=",esym,format(k,digits=4),",<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>w</i>=",esym,format(omega,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Passage Time Probabilities" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_PassageTimeProbability2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        if(is.null(yaxis)) { yaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(0,1),side="right")
        Ptxline <- list(color=grn$e,width=4)
        Ptline <- list(color=grn$c,width=4)
        legendy <- 1
        if(Pt[m,1] > 0.8 || Pt[m,n] > 0.8)
        {
          if((Pt[m,1] < 0.5 || Pt[m,n] < 0.5) && sigma/(2*rho)^0.5 < 15) { legendy <- 0.9 }
          else { legendy <- 0.15 }
        }
        kindex <- 0
        for(j in 1:n)  { if(z[j] == k) { kindex <- j } }
        legendpos <- list(orientation="h",x=1.0,y=legendy,xanchor="right")
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeProbability2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=t,y=Ptx,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=Ptxline,hoverinfo="x+y")
        dx <- as.integer((n-1)/10)
        if(dx < 1) { dx <- 1 }
        lineopacity <- 1
        j <- 1
        q <- 1
        fig <- add_trace(fig,type="scatter",x=t,y=Pt[,j],name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)"),mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="x+y",legendgroup="Pt",visible="legendonly",showlegend=TRUE)
        while(j < n)
        {
          j <- j+dx
          q <- q+1
          if(q < 7) { lineopacity <- lineopacity-0.07 }
          else { lineopacity <- lineopacity+0.07 }
          fig <- add_trace(fig,type="scatter",x=t,y=Pt[,j],name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|",z[j],")"),mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="x+y",legendgroup="Pt",visible="legendonly",showlegend=FALSE)
        }
        if(kindex > 0) { fig <- add_trace(fig,type="scatter",x=c(t[1],t[1]),y=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="x+y",legendgroup="Pt",visible="legendonly",showlegend=FALSE) }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_PassageTimeProbability3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        kindex <- 0
        coordinatex <- vector("character",m)
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          coordinatex[i] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(Ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(Pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j]) }
        }
        for(j in 1:n) { if(z[j] == k) { kindex <- j } }
        spy <- list(x=2.35,y=-0.85,z=0.5)
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*z[1]-0.03*z[n],1.03*z[n]-0.03*z[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(1.03*t[1]-0.03*t[m],1.03*t[m]-0.03*t[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(-0.03,1.03),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        xx <- vector("double",m)
        for(i in 1:m) { xx[i] <- x }
        Ptxmesh <- MeshCurtainSmooth(xx,t,Ptx,rep(0,m))
        Ptxline <- list(color=gry$d,width=8)
        Ptline <- list(color=grn$e,width=8)
        gradientPtx <- list(c(0,grn$a),c(1,grn$a))
        gradient <- list(c(0,grn$c),c(1,grn$c))
        if(k > mu) { rise <- list(x=10,y=100,z=0) }
        else if(k == mu) { rise <- list(x=0,y=100,z=0) }
        else { rise <- list(x=-10,y=100,z=0) }
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeProbability3D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,legendgroup="Ptx") %>%
          add_trace(.,type="mesh3d",x=Ptxmesh$xvertex,y=Ptxmesh$yvertex,z=Ptxmesh$zvertex,i=Ptxmesh$ivertex,j=Ptxmesh$jvertex,k=Ptxmesh$kvertex,intensity=Ptxmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientPtx,reversescale=reverse,opacity=0.7,hoverinfo="skip",legendgroup="Ptx",showlegend=FALSE)
        dx <- as.integer((n-1)/10)
        if(dx < 1) { dx <- 1 }
        lineopacity <- 1
        j <- 1
        q <- 1
        for(i in 1:m) { xx[i] <- z[j] }
        fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly")
        while(j < n)
        {
          j <- j+dx
          q <- q+1
          if(q < 7) { lineopacity <- lineopacity-0.07 }
          else { lineopacity <- lineopacity+0.07 }
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Ptz",visible="legendonly",showlegend=FALSE)
        }
        if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Ptz",visible="legendonly",showlegend=FALSE) }
        fig <-add_trace(fig,type="surface",x=z,y=t,z=Pt,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
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
    z_stoch_args = NULL,
    y_stoch_args = NULL,
    x_stoch_args = NULL,
    t_stoch_args = NULL,
    plot_args = NULL,
    undo_args = NULL,
    plot_info = NULL,
    flags = NULL,
    # private intermediate fields ----
    dOOdsconvexneg = NULL,
    dOOdsconvexpos = NULL,
    dOOdsconcaveneg = NULL,
    dOOdsconcavepos = NULL,
    dOOdspatchneg = NULL,
    dOOdspatchpos = NULL,
    # private output fields ----
    g = NULL,
    h2 = NULL,
    G = NULL,
    Gteps = NULL,
    H2 = NULL,
    H2teps = NULL,
    p = NULL,
    Pneg = NULL,
    Ppos = NULL,
    PPneg = NULL,
    PPpos = NULL,
    OOneg = NULL,
    OOpos = NULL,
    OOhatneg = NULL,
    OOhatpos = NULL,
    shatneg = NULL,
    shatpos = NULL,
    KOOneg = NULL,
    KOOpos = NULL,
    BCneg = NULL,
    BCpos = NULL,
    tmodemedianmean = NULL,
    tmodesmediansmeans = NULL,
    tpercentile = NULL,
    tpercentiles = NULL,
    ptx = NULL,
    pt = NULL,
    Ptx = NULL,
    Pt = NULL,
    # private globals ----
    psiphi = NULL,
    undoIx = NULL,
    modebar_2D = NULL,
    modebar_3D = NULL,
    plot_types = NULL,
    # private colors ----
    plot_colors = NULL,
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
    # private calculate methods ----
    OUPdOOdsZero = function()
    {
      # get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
      # calculate ----
      n <- length(x)
      dOOdszero <- RcppOUPAdOOdsZero(s,x,t,y,rho,mu,sigma,r,phi,b,c)
      dOOdsconvex <- dOOdszero[1:2,1:n,drop=FALSE]
      dOOdsconcave <- dOOdszero[3:4,1:n,drop=FALSE]
      dOOdspatch <- dOOdszero[1:3,(n+1):(n+3),drop=FALSE]
      if(phi > 0 )
      {
        private$dOOdsconvexpos <- dOOdsconvex
        private$dOOdsconcavepos <- dOOdsconcave
        private$dOOdspatchpos <- dOOdspatch
      }
      else
      {
        private$dOOdsconvexneg <- dOOdsconvex
        private$dOOdsconcaveneg <- dOOdsconcave
        private$dOOdspatchneg <- dOOdspatch
      }
      return(list(dOOdsconvex=dOOdsconvex,dOOdsconcave=dOOdsconcave,dOOdspatch=dOOdspatch))
    },
    # private plot methods ----
    MeshCurtainSmooth = function(x,y,ztop,zbottom)
    {
      m <- min(length(x),length(y),length(ztop),length(zbottom))
      xv <- vector("double",3*m-1)
      yv <- vector("double",3*m-1)
      zv <- vector("double",3*m-1)
      iv <- vector("double",4*(m-1))
      jv <- vector("double",4*(m-1))
      kv <- vector("double",4*(m-1))
      i <- 0
      while(i < m)
      {
        i <- i+1
        xv[i] <- x[i]
        yv[i] <- y[i]
        zv[i] <- zbottom[i]
        xv[i+m] <- x[i]
        yv[i+m] <- y[i]
        zv[i+m] <- ztop[i]
      }
      i <- 0
      while(i < m-1)
      {
        i <- i+1
        xv[i+2*m] <- 0.5*(x[i]+x[i+1])
        yv[i+2*m] <- 0.5*(y[i]+y[i+1])
        zv[i+2*m] <- 0.25*(ztop[i]+zbottom[i]+ztop[i+1]+zbottom[i+1])
      }
      i <- 0
      while(i < m-1)
      {
        i <- i+1
        # right-hand side indexes are zero-based
        iv[4*i-3] <- i-1+2*m
        jv[4*i-3] <- i-1
        kv[4*i-3] <- i
        iv[4*i-2] <- i-1+2*m
        jv[4*i-2] <- i
        kv[4*i-2] <- i+m
        iv[4*i-1] <- i-1+2*m
        jv[4*i-1] <- i+m
        kv[4*i-1] <- i-1+m
        iv[4*i] <- i-1+2*m
        jv[4*i] <- i-1+m
        kv[4*i] <- i-1
      }
      return(list(xvertex=xv,yvertex=yv,zvertex=zv,ivertex=iv,jvertex=jv,kvertex=kv))
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
    OUPDensity = function(s,x,t,y,rho,mu,sigma,dy)
    {
      mean <- mu+(x-mu)*exp(-rho*(t-s))
      if(rho < 0.0000000001) { variance <- sigma^2*(t-s) }
      else { variance <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
      if(variance^0.5 < dy)
      {
        lomean <- mean-0.5*dy
        himean <- mean+0.5*dy
        if(y >= lomean && y < himean)
        {
          if(dy > 0) { density <- 1/dy }
          else { density <- 10 }
        }
        else { density <- 0 }
      }
      else
      {
        u <- 1/variance^0.5
        v <- 0.5*(y-mean)^2/variance
        density <- 1/(2^0.5*1.77245385090552)*u*exp(-v)
      }
      return(density)
    },
    OUPDnsty = function(s,x,t,y,rho,mu,sigma)
    {
      # Do not call this function with t=s or sigma=0
      # This routine is for calculating options where the checks are done
      mean <- mu+(x-mu)*exp(-rho*(t-s))
      if(rho < 0.0000000001) { variance <- sigma^2*(t-s) }
      else { variance <- sigma^2/(2*rho)*(1-exp(-2*rho*(t-s))) }
      if(variance < 0.0000000001) { density <- 0 }
      else
      {
        u <- 1/variance^0.5
        v <- 0.5*(y-mean)^2/variance
        density <- 1/(2^0.5*1.77245385090552)*u*exp(-v)
      }
      return(density)
    },
    OUPProbability = function(s,x,t,y,rho,mu,sigma,psi)
    {
      epsilon <- 0.000000000001
      if( t <= s)
      {
        if(y < x-epsilon) { probability <- 0 }
        else if(y <= x+epsilon) { probability <- 0.5 }
        else { probability <- 1 }
      }
      else
      {
        mean <- mu+(x-mu)*exp(-rho*(t-s))
        if(rho < 0.0000000001) { variance <- sigma^2*(t-s) }
        else { variance <- sigma^2/(2*rho)*(1-exp(-2*rho*(t-s))) }
        if(variance < 0.0000000001)
        {
          if(y <= mean-epsilon) { probability <- 0 }
          else if(y <= mean+epsilon) { probability <- 0.5 }
          else { probability <- 1 }
        }
        else
        {
          v2 <- 0.5*(y-mean)^2/variance
          if(y < mean) { probability <- 0.5*private$GammaBigOneHalf(v2)/1.77245385090552 }
          else if(y == mean) { probability <- 0.5 }
          else { probability <- 0.5*(1+private$GammaSmallOneHalf(v2)/1.77245385090552) }
        }
      }
      if(psi > 0) { probability <- 1-probability }
      return(probability)
    },
    OUPPrbblty = function(s,x,t,y,rho,mu,sigma,psi)
    {
      # This routine has fewer checks for use in options where the checks are done
      epsilon <- 0.000000000001
      mean <- mu+(x-mu)*exp(-rho*(t-s))
      if(rho < 0.0000000001) { variance <- sigma^2*(t-s) }
      else { variance <- sigma^2/(2*rho)*(1-exp(-2*rho*(t-s))) }
      if(variance < 0.0000000001)
      {
        if(y <= mean-epsilon) { probability <- 0 }
        else if(y <= mean+epsilon) { probability <- 0.5 }
        else { probability <- 1 }
      }
      else
      {
        v2 <- 0.5*(y-mean)^2/variance
        if(y < mean) { probability <- 0.5*private$GammaBigOneHalf(v2)/1.77245385090552 }
        else if(y == mean) { probability <- 0.5 }
        else { probability <- 0.5*(1+private$GammaSmallOneHalf(v2)/1.77245385090552) }
      }
      if(psi > 0) { probability <- 1-probability }
      return(probability)
    },
    OUPOption = function(s,x,t,y,rho,mu,sigma,r,phi,b,c)
    {
      ltgt <- 1
      bc <- -c
      if(phi > 0)
      {
        ltgt = -1
        bc <- b
      }
      if(t <= s) { option <- ltgt*(y-x)*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,phi) }
      else if(sigma == 0) { option <- ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,phi) }
      else if(rho < 0.0000000001) { option <- ltgt*(y-x)*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,phi)+sigma^2*(t-s)*private$OUPDnsty(s,x,t,y,rho,mu,sigma) }
      else { option <- ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,phi)+sigma^2/(2*rho)*(1-exp(-2*rho*(t-s)))*private$OUPDnsty(s,x,t,y,rho,mu,sigma) }
      option <- exp(-r*(t-s))*(option+bc)
      return(option)
    },
    OUPPassageTimeDensity = function(s,x,t,k,omega,rho,mu,sigma,dt)
    {
      if(k == Inf || k == -Inf) { density <- 0 }
      else
      {
        if(rho < 0.0000000001) { varx <- sigma^2*(t-s) }
        else { varx <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
        if(varx < 0.0000000001)
        {
          if(x < mu)
          {
            lomeanx <- mu+(x-mu)*exp(-rho*(t-0.5*dt-s))
            himeanx <- mu+(x-mu)*exp(-rho*(t+0.5*dt-s))
          }
          else
          {
            lomeanx <- mu+(x-mu)*exp(-rho*(t+0.5*dt-s))
            himeanx <- mu+(x-mu)*exp(-rho*(t-0.5*dt-s))
          }
          if(k > lomeanx && k <= himeanx) { density <- 9 }
          else { density <- 0 }
        }
        else if(x == k)
        {
          u <- abs(k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))*sigma^2/varx^1.5
          v2 <- 0.5*((k-mu)*(1-exp(-rho*(t-s))))^2/varx
          density <- (1-omega)*u*exp(-v2)/(2*2.50662827431001)
        }
        else if(k == mu)
        {
          u <- abs(x-k)*exp(-rho*(t-s))*sigma^2/varx^1.5
          v2 <- 0.5*((k-x)*exp(-rho*(t-s)))^2/varx
          density <- (1+omega)*u*exp(-v2)/(2*2.50662827431001)
        }
        else
        {
          meanx <- mu+(x-mu)*exp(-rho*(t-s))
          boundx <- mu+(k-mu)*exp(-rho*(t-s))-(x-k)*exp(-rho*(t-s))
          vm2 <- 0.5*(k-meanx)^2/varx
          vb2 <- 0.5*(k-boundx)^2/varx
          smallgammab <- private$GammaSmallOneHalf(vb2)
          biggammab <- private$GammaBigOneHalf(vb2)
          lnlambda <- 2*(k-mu)*(1-exp(-rho*(t-s)))*(x-k)*exp(-rho*(t-s))/varx
          dlnlambda <- -(k-mu)*(x-k)*sigma^2*(exp(-rho*(t-s))-2*exp(-2*rho*(t-s))+exp(-3*rho*(t-s)))/varx^2
          if(x > k)
          {
            u <- ((x-mu)*exp(-rho*(t-s))-(k-mu)*exp(-2*rho*(t-s))+omega*((x-k)*exp(-rho*(t-s))-(k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))))*sigma^2/varx^1.5
            if(k < mu && k < boundx) { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552) }
            else
            {
              if(lnlambda > 709.782712893384) { density <- u*exp(-vm2)/(2*2.50662827431001) }
              else { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*biggammab*dlnlambda*exp(lnlambda)/(2*1.77245385090552) }
            }
          }
          else
          {
            u <- ((k-mu)*exp(-2*rho*(t-s))-(x-mu)*exp(-rho*(t-s))+omega*((k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))-(x-k)*exp(-rho*(t-s))))*sigma^2/varx^1.5
            if(k > mu && k > boundx) { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552) }
            else
            {
              if(lnlambda > 709.782712893384) { density <- u*exp(-vm2)/(2*2.50662827431001) }
              else { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*biggammab*dlnlambda*exp(lnlambda)/(2*1.77245385090552) }
            }
          }
        }
      }
      return(density)
    },
    OUPPassageTimeProbabilityInf = function(x,k,omega,rho,mu,sigma)
    {
      if(k == Inf || k == -Inf) { PInf <- 0 }
      else
      {
        if(omega == 1) { PInf <- 1 }
        else if(k == mu) { PInf <- 0.5*(1+omega) }
        else if(sigma^2 < 0.0000000001)
        {
          if(x == k) { PInf <- 1 }
          else if(x > k)
          {
            if(k > mu) { PInf <- 1 }
            else { PInf <- omega }
          }
          else
          {
            if(k < mu) { PInf <- 1 }
            else { PInf <- omega }
          }
        }
        else
        {
          v2 <- rho*((k-mu)/sigma)^2
          smallgamma <- private$GammaSmallOneHalf(v2)
          biggamma <- private$GammaBigOneHalf(v2)
          if(x == k) { PInf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
          else if(x > k)
          {
            if(k > mu) { PInf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
            else { PInf <- (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552) }
          }
          else if(x < k)
          {
            if(k < mu) { PInf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
            else { PInf <- (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552) }
          }
        }
      }
      return(PInf)
    },
    GammaSmallOneHalf = function(x)
    {
      sigdig <- 13
      n <- 2000
      if(x <= 0) { gamma <- 0 }
      else if(x >= 33.8551300352457) { gamma <- 1.77245385090552 }
      else
      {
        dgm <- 1
        gm <- dgm
        i <- 0
        while(10^sigdig*dgm >= gm && i < 2000)
        {
          i <- i+1
          dgm <- dgm*x/(1.5+i)
          gm <- gm+dgm
        }
        if(i < n)
        {
          gm <- x^(1.5)/(1.5)*exp(-x)*gm
          gm <- ((x^0.5)*exp(-x)+gm)/0.5
          gamma <- gm
        }
        else { gamma <- NaN }
      }
      return(gamma)
    },
    GammaBigOneHalf = function(x)
    {
      sigdig <- 13
      n <- 2000
      if(x <= 0) { gamma <- 1.77245385090552 }
      else if(x > 708.396418532264) { gamma <- 0.0 }
      else if(x<27)
      {
        dgm <- 1
        gm <- dgm
        i <- 0
        while(10^sigdig*dgm >= gm && i < n)
        {
          i <- i+1
          dgm <- dgm*x/(1.5+i)
          gm <- gm+dgm
        }
        if(i < n)
        {
          gm <- x^(1.5)/(1.5)*exp(-x)*gm
          gm <- ((x^0.5)*exp(-x)+gm)/0.5
          gamma <- 1.77245385090552-gm
        }
        else { gamma <- NaN }
      }
      else
      {
        n <- round(x,0)
        dgm <- 1
        gm <- 1
        i <- 0
        while(i<n)
        {
          i <- i+1
          dgm <- dgm*(0.5-i)/x
          gm <- gm+dgm
        }
        gamma <- x^(-0.5)*exp(-x)*gm
      }
      return(gamma)
    },
    # private clipboard methods ----
    CopyToClipboard = function(clip)
    {
      if(!is.null(private$OUP)) { OUP$CopyToClipboard(clip) }
      else if(interactive() && clipr_available()) { write_clip(clip,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    }
  )
)
