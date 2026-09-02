library(R6)
library(plotly)
library(stringr)
library(clipr)

# roxygen ----
#' @title R6 class implementing a Finite Difference method for the Ornstein-Uhlenbeck Process.
#'
#' @description
#' A Finite Difference method--without boundary conditions on the state but with
#'  arbitrary terminal values--for the pricing of index insurance, perpetual options
#'  and sequential options. From the option prices, an option envelope and a decision
#'  threshold can be calculated.  Drift, diffusion and several possible terminal values
#'  are pre-programmed.  User-defined terminal values can also be entered.
#'
#' @details # Formulas and Methods:
#'     x stochastic
#'       Drift
#'       Diffusion
#'       TerminalValue_Linear
#'       TerminalValue_Degenerate
#'       TerminalValue_Stepped
#'       TerminalValue_Kinked
#'       TerminalValue_Butterfly
#'       TerminalValue_Mitscherlich
#'       TerminalValue_Gompertz
#'       TerminalValue_Logistic
#'       TerminalValue_Transcendental
#'       TerminalValue_YieldIndex
#'       TerminalValue
#'       Option
#'       OptionEnvelope
#'       DecisionThreshold
#'
#' @details # Plots
#'       PlotDrift
#'       PlotDiffusion
#'       PlotTerminalValue
#'       PlotOption
#'       PlotOptionEnvelope
#'       PlotDecisionThreshold
#'
#' @details # Arguments of functions:
#'       All arguments are optional in all functions.
#'       s:     vector of times
#'       x:     vector of states
#'       V:     vector of terminal values
#'       r:     discount rate
#'       phi:   search direction for exit or entry options
#'       theta: weight of current time in time stepping
#'       skip:  divide the time interval and report every skip result
#'       rho:   rate parameter
#'       mu:    location parameter
#'       sigma: scale parameter
#'       xo:    state at the intercept, spike, step or kink
#'       xi:    state at the inflection point
#'       xm:    state at the maximum or kink
#'       vs:    slope or direction of a step
#'       vr:    rate of change
#'       Vmax:  maximum terminal value
#'       Vmin:  minimum terminal value
#'
#' @details # Usage:
#' The FiniteDifference object must first be instantiated before its methods
#'  are called. There are two ways.  The first way instantiates the OUProcess
#'  object and calls a function to get a pointer:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       FD <- OUP$get_FiniteDifference()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'
#' The FiniteDifference object will coordinate arguments to functions with the
#'  Analytial, MaximumLikelihood and MonteCarlo objects.  The second
#'  way instantiates the FiniteDifference object by itself with no coordination:
#'
#'       FD <- FiniteDifference$new()
#'
#' Once the object is instantiated, its methods can be called, to calculate and
#'  plot a Decision Threshold, for example:
#'
#'       FD$DecisionThreshold()
#'
#' The plot methods can be used to customize the plots, with a title, for example:
#'
#'       FD$PlotDecisionThreshold(title="My Decision")
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
#' The Finite Difference Method can solve problems with no analytical solutions.
#'  As a numerical procedure, it will have errors. The larger the scale parameter,
#'  the larger the errors.  To manage the errors, the domain for x should be as
#'  large as practical. The increment between times is ds.  Parameter skip
#'  subdivides ds into smaller intervals. Only calculations for ds are reported.
#'  The increment between states is dx and ds/skip should be about dx/100.
#'  Increasing skip or dx may reduce errors up to a point. If possible, the
#'  Finite Difference Method should be calibrated against an analytical solution
#'  for a simpler problem before solving more complicated problems.
#'
#' The functions include many different types of terminal values, but suppose
#'  the bottom line is a simple exit option:
#'
#'       FD <- FiniteDifference$new()
#'       FD$DecisionThreshold(mu=10)
#'
#' Or a simple entry option:
#'
#'       FD$TerminalValue_Kinked(xo=5,vs=1)
#'       FD$DecisionThreshold()
#'
#' (That's x oh, not x zero). And finally a sequence of an entry with an option
#'  to exit, calculating the exit decision first:
#'
#'       FD$TerminalValue_Kinked(xo=0,vs=-1)
#'       exit <- FD$OptionEnvelope()$Ohat
#'       V <- exit+FD$TerminalValue_Kinked(xo=5,vs=1)$V
#'       FD$DecisionThreshold(V=V,phi=1)
#'
#' The functions all return named lists. Before being used in calculations, a
#'  list must be stripped of its names.  In the above example, the names are Ohat
#'  for the option envelope and V for the terminal value. The functions try to
#'  help by accepting named lists as arguments.  For example:
#'
#'       V <- FD$TerminalValue_Kinked(xo=5,vs=1)
#'       FD$DecisionThreshold(V=V)
#'
#' But the terminal value functions automatically set the terminal value in the
#'  object, so entering V as as argument does nothing useful in this case.
#'
#' You may wish to calculate your own terminal value in the console. You must
#'  first get rid of the names. The names to get rid of are listed under
#'  'Returns' in the help for each function. You can also use double brackets
#'  to eliminate names if you know the position in the list. Both 'OOhat' and 'V'
#'  are first in their lists:
#'
#'       FD$TerminalValue_Kinked(xo=0,vs=-1)
#'       exit <- FD$OptionEnvelope()[[1]]
#'       V <- exit+FD$TerminalValue_Kinked(xo=5,vs=1)[[1]]
#'       FD$DecisionThreshold(V=V,phi=1)
#'
#' Or you can unlist everything:
#'
#'       V <- unlist(FD$OptionEnvelope())+unlist(FD$TerminalValue_Kinked(xo=5,vs=1))
#'
#' Sorry.

# class ----
FiniteDifference <- R6::R6Class("FiniteDifference",
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
    #' Create a new FiniteDifference object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new FiniteDifference object
    initialize = function(OUP=NULL)
    {
      # pointer to container object ----
      if(!is.null(OUP) && class(OUP)[[1]] == "OUProcess") { private$OUP <- OUP }
      # arguments ----
      s <- seq(from=10,to=0,by=-0.1)
      x <- seq(from=-30,to=30,by=0.6)
      V <- c(seq(from=30,to=0.6,by=-0.6),rep(0,51))
      private$oup_params <- list(rho=0.5,mu=15,sigma=15)
      private$x_stoch_args <- list(s=s,x=x,V=V,r=0.05,phi=0,theta=0.5,skip=10,ds=0.1,dx=0.6)
      private$V_linear_args <- list(xo=0,vs=1)
      private$V_degenerate_args <- list(xo=0,Vmax=0.5,Vmin=0)
      private$V_stepped_args <- list(xo=0,vs=-1,Vmax=1,Vmin=0)
      private$V_kinked_args <- list(xo=0,vs=-1,Vmax=30,Vmin=0)
      private$V_butterfly_args <- list(xo=-5,xm=5,vs=-1,Vmax=25,Vmin=0)
      private$V_mitscherlich_args <- list(xo=0,vr=-0.1,Vmax=10,Vmin=0)
      private$V_gompertz_args <- list(xi=-5,vr=-0.1,Vmax=10,Vmin=0)
      private$V_logistic_args <- list(xi=-5,vr=-0.1,Vmax=10,Vmin=0)
      private$V_transcendental_args <- list(xo=25,xi=10,xm=-5,Vmax=10,Vmin=0)
      private$V_yieldindex_args <- list(xo=25,xi=10,xm=-5,Vmax=8,Vmin=-2)
      private$V_info <- list(Ix=4,name="Kinked",names=list("Linear","Degenerate","Stepped","Kinked","Butterfly","Mitscherlich","Gompertz","Logistic","Transcendental","Yield Index"),text="Linear, Degenerate, Stepped, Kinked, Butterfly, Mitscherlich, Gompertz, Logistic, Transcendental, Yield Index")
      # undo ----
      private$undo_args <- list(list(oup_params=private$oup_params,x_stoch_args=private$x_stoch_args,V_linear_args=private$V_linear_args,V_degenerate_args=private$V_degenerate_args,V_stepped_args=private$V_stepped_args,V_kinked_args=private$V_kinked_args,V_butterfly_args=private$V_butterfly_args,V_mitscherlich_args=private$V_mitscherlich_args,V_gompertz_args=private$V_gompertz_args,V_logistic_args=private$V_logistic_args,V_transcendental_args=private$V_transcendental_args,V_yieldindex_args=private$V_yieldindex_args,V_info=private$V_info))
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
      private$plot_types <- c(rep(0,3))
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
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"FD") }
      if(!is.null(rho))
      {
        sca <- private$extract_scalar(rho)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            message("negative rho set to zero.")
            sca <- 0.0
          }
          if(sca != private$oup_params$rho)
          {
            private$oup_params$rho <- sca
            private$g <- NULL
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
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
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
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
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
        }
        else { message("sigma not set.")}
      }
      return(private$oup_params)
    },
    #' @description
    #' Set x as a stochastic state and its arguments
    #' @param s     vector of times
    #' @param x     vector of states
    #' @param V     vector of terminal values
    #' @param r     discount rate 0<=r<inf, scalar
    #' @param phi   search direction for exit or entry options
    #' @param theta weight of current time in time stepping 0.5<=theta<=1
    #' @param skip  subdivide time intervals but report at times s 1<=skip<=1000
    #' @param who   object id of sender
    #' @return list(s,x,V,r,phi,theta,skip,ds,dx)
    set_x_stoch_args = function(s=NULL,x=NULL,V=NULL,r=NULL,phi=NULL,theta=NULL,skip=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,NULL,r,phi,"FD") }
      if(!is.null(s))
      {
        vec <- private$extract_vector(s,-1)
        if(!is.null(vec))
        {
          m <- length(vec)
          if(m > 1)
          {
            ds <- vec[1]-vec[2]
            if(ds > 0)
            {
              ok <- TRUE
              i=2
              while(i < m && ok == TRUE)
              {
                i <- i+1
                if(abs(vec[i-1]-vec[i]-ds) > 0.00001)
                {
                  message("time increments are not equal.")
                  message("s not set.")
                  ok <- FALSE
                }
              }
              if(ok == TRUE)
              {
                if(!private$vecs_equal(vec,private$x_stoch_args$s))
                {
                  private$x_stoch_args$s <- vec
                  private$x_stoch_args$ds <- ds
                  private$OO <- NULL
                }
              }
            }
            else
            {
              message("time increment is zero.")
              message("s not set.")
            }
          }
          else
          {
            message("must have at least two times.")
            message("s not set.")
          }
        }
        else { message("s not set.")}
      }
      if(!is.null(x))
      {
        xx <- private$extract_vector(x)
        if(!is.null(xx))
        {
          nx <- length(xx)
          if(nx > 100)
          {
            index <- order(xx)
            xx <- xx[index]
            dx <- xx[2]-xx[1]
            if(dx > 0)
            {
              i <- 2
              cancel <- FALSE
              while(i < nx && cancel == FALSE)
              {
                i <- i+1
                if(abs(xx[i]-xx[i-1]-dx) > 0.00001)
                {
                  message("x increments are not equal.")
                  cancel <- TRUE
                  xx <- NULL
                }
              }
            }
            else
            {
              message("x increment is zero.")
              xx <- NULL
            }
          }
          else
          {
            message("x vector must have at least 101 elements")
            xx <- NULL
          }
        }
        if(!is.null(xx))
        {
          if(!private$vecs_equal(xx,private$x_stoch_args$x))
          {
            private$x_stoch_args$x <- xx
            private$x_stoch_args[3] <- list(V=NULL)
            private$x_stoch_args$dx <- dx
            private$g <- NULL
            private$h2 <- NULL
            private$V_linear <- NULL
            private$V_degenerate <- NULL
            private$V_stepped <- NULL
            private$V_kinked <- NULL
            private$V_butterfly <- NULL
            private$V_mitscherlich <- NULL
            private$V_gompertz <- NULL
            private$V_logistic <- NULL
            private$V_transcendental <- NULL
            private$V_yieldindex <- NULL
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
          if(!is.null(V))
          {
            VV <- private$extract_vector(V)
            if(!is.null(VV))
            {
              nV <- length(VV)
              if(nV < nx)
              {
                message("vector V is shorter than vector x.")
                VV <- NULL
              }
              else if(nV > nx)
              {
                message("vector V is longer than vector x.")
                VV <- NULL
              }
            }
            if(!is.null(VV))
            {
              VVV <- VV[index]
              if(!private$vecs_equal(VVV,private$x_stoch_args$V))
              {
                private$x_stoch_args$V <- VVV
                private$V_info$name <- "Terminal Value"
              }
            }
            else { message("V not set.")}
          }
        }
        else { message("x not set.") }
      }
      else if(!is.null(V))
      {
        VV <- private$extract_vector(V)
        if(!is.null(VV))
        {
          xx <- private$x_stoch_args$x
          nx <- length(xx)
          nV <- length(VV)
          if(nV < nx)
          {
            message("vector V is shorter than vector x.")
            VV <- NULL
          }
          else if(nV > nx)
          {
            message("vector V is longer than vector x.")
            VV <- NULL
          }
        }
        if(!is.null(VV))
        {
          if(!private$vecs_equal(VV,private$x_stoch_args$V))
          {
            private$x_stoch_args$V <- VV
            private$V_info$name <- "Terminal Value"
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
        }
        else { message("V not set.")}
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
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
        }
        else { message("r not set.")}
      }
      if(!is.null(phi))
      {
        sca <- private$extract_scalar(phi)
        if(!is.null(sca))
        {
          if(sca != private$x_stoch_args$phi)
          {
            private$x_stoch_args$phi <- sca
            private$kOOhat <- NULL
          }
        }
        else { message("phi not set.")}
      }
      if(!is.null(theta))
      {
        sca <- private$extract_scalar(theta)
        if(!is.null(sca))
        {
          if(sca < 0.5)
          {
            message("theta has been set to 0.5.")
            sca <- 0.5
          }
          else if(sca > 1)
          {
            message("theta has been set to 1.")
            sca <- 1
          }
          if(sca != private$x_stoch_args$theta)
          {
            private$x_stoch_args$theta <- sca
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
        }
        else { message("theta not set.")}
      }
      if(!is.null(skip))
      {
        sca <- private$extract_scalar(skip)
        if(!is.null(sca))
        {
          sca <- as.integer(sca)
          if(sca < 1)
          {
            message("skip has been set to 1.")
            sca <- 1
          }
          else if(sca > 1000)
          {
            message("skip has been set to 1000.")
            sca <- 1000
          }
          if(sca != private$x_stoch_args$skip)
          {
            private$x_stoch_args$skip <- sca
            private$OO <- NULL
            private$OOhat <- NULL
            private$kOOhat <- NULL
          }
        }
        else { message("skip not set.")}
      }
      return(private$x_stoch_args)
    },
    #' @description
    #' Set V linear arguments
    #' @param xo state at the intercept
    #' @param vs slope
    #' @return list(xo,vs)
    set_V_linear_args = function(xo=NULL,vs=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_linear_args$xo)
          {
            private$V_linear_args$xo <- sca
            private$V_linear <- NULL
            if(private$V_info[[1]] == 1) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(vs))
      {
        sca <- private$extract_scalar(vs)
        if(!is.null(sca))
        {
          if(sca != private$V_linear_args$vs)
          {
            private$V_linear_args$vs <- sca
            private$V_linear <- NULL
            if(private$V_info[[1]] == 1) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vs not set.")}
      }
      return(private$V_linear_args)
    },
    #' @description
    #' Set V degenerate arguments
    #' @param xo   state at the step
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,vs,Vmax,Vmin)
    set_V_degenerate_args = function(xo=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_degenerate_args$xo)
          {
            private$V_degenerate_args$xo <- sca
            private$V_degenerate <- NULL
            if(private$V_info[[1]] == 2) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_degenerate_args$Vmax || mn != private$V_degenerate_args$Vmin)
            {
              private$V_degenerate_args$Vmax <- mx
              private$V_degenerate_args$Vmin <- mn
              private$V_degenerate <- NULL
              if(private$V_info[[1]] == 2) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_degenerate_args[[3]]
          if(mx > mn)
          {
            if(mx != private$V_degenerate_args$Vmax)
            {
              private$V_degenerate_args$Vmax <- mx
              private$V_degenerate <- NULL
              if(private$V_info[[1]] == 2) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_degenerate_args[[2]]
          if(mn < mx)
          {
            if(mn != private$V_degenerate_args$Vmin)
            {
              private$V_degenerate_args$Vmin <- mn
              private$V_degenerate <- NULL
              if(private$V_info[[1]] == 2) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_degenerate_args)
    },
    #' @description
    #' Set V stepped arguments
    #' @param xo   state at the step
    #' @param vs   direction of the step
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,vs,Vmax,Vmin)
    set_V_stepped_args = function(xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_stepped_args$xo)
          {
            private$V_stepped_args$xo <- sca
            private$V_stepped <- NULL
            if(private$V_info[[1]] == 3) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(vs))
      {
        sca <- private$extract_scalar(vs)
        if(!is.null(sca))
        {
          if(sca != private$V_stepped_args$vs)
          {
            private$V_stepped_args$vs <- sca
            private$V_stepped <- NULL
            if(private$V_info[[1]] == 3) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vs not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_stepped_args$Vmax || mn != private$V_stepped_args$Vmin)
            {
              private$V_stepped_args$Vmax <- mx
              private$V_stepped_args$Vmin <- mn
              private$V_stepped <- NULL
              if(private$V_info[[1]] == 3) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_stepped_args[[4]]
          if(mx > mn)
          {
            if(mx != private$V_stepped_args$Vmax)
            {
              private$V_stepped_args$Vmax <- mx
              private$V_stepped <- NULL
              if(private$V_info[[1]] == 3) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_stepped_args[[3]]
          if(mn < mx)
          {
            if(mn != private$V_stepped_args$Vmin)
            {
              private$V_stepped_args$Vmin <- mn
              private$V_stepped <- NULL
              if(private$V_info[[1]] == 3) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_stepped_args)
    },
    #' @description
    #' Set V kinked arguments
    #' @param xo   state at the kink
    #' @param vs   slope
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,vs,Vmax,Vmin)
    set_V_kinked_args = function(xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_kinked_args$xo)
          {
            private$V_kinked_args$xo <- sca
            private$V_kinked <- NULL
            if(private$V_info[[1]] == 4) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(vs))
      {
        sca <- private$extract_scalar(vs)
        if(!is.null(sca))
        {
          if(sca != private$V_kinked_args$vs)
          {
            private$V_kinked_args$vs <- sca
            private$V_kinked <- NULL
            if(private$V_info[[1]] == 4) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vs not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_kinked_args$Vmax || mn != private$V_kinked_args$Vmin)
            {
              private$V_kinked_args$Vmax <- mx
              private$V_kinked_args$Vmin <- mn
              private$V_kinked <- NULL
              if(private$V_info[[1]] == 4) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_kinked_args[[4]]
          if(mx > mn)
          {
            if(mx != private$V_kinked_args$Vmax)
            {
              private$V_kinked_args$Vmax <- mx
              private$V_kinked <- NULL
              if(private$V_info[[1]] == 4) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_kinked_args[[3]]
          if(mn < mx)
          {
            if(mn != private$V_kinked_args$Vmin)
            {
              private$V_kinked_args$Vmin <- mn
              private$V_kinked <- NULL
              if(private$V_info[[1]] == 4) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_kinked_args)
    },
    #' @description
    #' Set V butterfly arguments
    #' @param xo   state at the left wing
    #' @param xm   state at the right wing
    #' @param vs   slope
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,xm,vs,Vmax,Vmin)
    set_V_butterfly_args = function(xo=NULL,xm=NULL,vs=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo) && !is.null(xm))
      {
        scao <- private$extract_scalar(xo)
        scam <- private$extract_scalar(xm)
        if(!is.null(scao) && !is.null(scam))
        {
          if(scao <= scam)
          {
            if(scao != private$V_butterfly_args$xo || scam != private$V_butterfly_args$xm)
            {
              private$V_butterfly_args$xo <- scao
              private$V_butterfly_args$xm <- scam
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("xo > xm.")
            message("x0 and xm not set.")
          }
        }
        else { message("xo and xm not set.")}
      }
      else if(!is.null(xo))
      {
        scao <- private$extract_scalar(xo)
        if(!is.null(scao))
        {
          scam <- private$V_butterfly_args[[2]]
          if(scao <= scam)
          {
            if(scao != private$V_butterfly_args$scao)
            {
              private$V_butterfly_args$xo <- scao
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("xo > mx=",scam))
            message("xo not set.")
          }
        }
        else { message("xo not set.")}
      }
      else if(!is.null(xm))
      {
        scam <- private$extract_scalar(xm)
        if(!is.null(scam))
        {
          scao <- private$V_butterfly_args[[1]]
          if(scam >= scao)
          {
            if(scam != private$V_butterfly_args$xm)
            {
              private$V_butterfly_args$xm <- scam
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("xm < existing xo=",scao))
            message("xm not set.")
          }
        }
        else { message("xm not set.")}
      }
      if(!is.null(vs))
      {
        sca <- private$extract_scalar(vs)
        if(!is.null(sca))
        {
          if(sca != private$V_butterfly_args$vs)
          {
            private$V_butterfly_args$vs <- sca
            private$V_butterfly <- NULL
            if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vs not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_butterfly_args$Vmax || mn != private$V_butterfly_args$Vmin)
            {
              private$V_butterfly_args$Vmax <- mx
              private$V_butterfly_args$Vmin <- mn
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_butterfly_args[[5]]
          if(mx > mn)
          {
            if(mx != private$V_butterfly_args$Vmax)
            {
              private$V_butterfly_args$Vmax <- mx
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_butterfly_args[[4]]
          if(mn < mx)
          {
            if(mn != private$V_butterfly_args$Vmin)
            {
              private$V_butterfly_args$Vmin <- mn
              private$V_butterfly <- NULL
              if(private$V_info[[1]] == 5) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_butterfly_args)
    },
    #' @description
    #' Set V mitscherlich arguments
    #' @param xo   state at the intercept
    #' @param vr   rate of change
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,vr,Vmax,Vmin)
    set_V_mitscherlich_args = function(xo=NULL,vr=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_mitscherlich_args$xo)
          {
            private$V_mitscherlich_args$xo <- sca
            private$V_mitscherlich <- NULL
            if(private$V_info[[1]] == 6) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(vr))
      {
        sca <- private$extract_scalar(vr)
        if(!is.null(sca))
        {
          if(sca != private$V_mitscherlich_args$vr)
          {
            private$V_mitscherlich_args$vr <- sca
            private$V_mitscherlich <- NULL
            if(private$V_info[[1]] == 6) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vr not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_mitscherlich_args$Vmax || mn != private$V_mitscherlich_args$Vmin)
            {
              private$V_mitscherlich_args$Vmax <- mx
              private$V_mitscherlich_args$Vmin <- mn
              private$V_mitscherlich <- NULL
            if(private$V_info[[1]] == 6) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_mitscherlich_args[[4]]
          if(mx > mn)
          {
            if(mx != private$V_mitscherlich_args$Vmax)
            {
              private$V_mitscherlich_args$Vmax <- mx
              private$V_mitscherlich <- NULL
              if(private$V_info[[1]] == 6) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_mitscherlich_args[[3]]
          if(mn < mx)
          {
            if(mn != private$V_mitscherlich_args$Vmin)
            {
              private$V_mitscherlich_args$Vmin <- mn
              private$V_mitscherlich <- NULL
              if(private$V_info[[1]] == 6) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_mitscherlich_args)
    },
    #' @description
    #' Set V gompertz arguments
    #' @param xi   state at the inflection point
    #' @param vr   rate of change
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xi,vr,Vmax,Vmin)
    set_V_gompertz_args = function(xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xi))
      {
        sca <- private$extract_scalar(xi)
        if(!is.null(sca))
        {
          if(sca != private$V_gompertz_args$xi)
          {
            private$V_gompertz_args$xi <- sca
            private$V_gompertz <- NULL
            if(private$V_info[[1]] == 7) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xi not set.")}
      }
      if(!is.null(vr))
      {
        sca <- private$extract_scalar(vr)
        if(!is.null(sca))
        {
          if(sca != private$V_gompertz_args$vr)
          {
            private$V_gompertz_args$vr <- sca
            private$V_gompertz <- NULL
            if(private$V_info[[1]] == 7) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vr not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_gompertz_args$Vmax || mn != private$V_gompertz_args$Vmin)
            {
              private$V_gompertz_args$Vmax <- mx
              private$V_gompertz_args$Vmin <- mn
              private$V_gompertz <- NULL
              if(private$V_info[[1]] == 7) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_gompertz_args[[4]]
          if(mx > mn)
          {
            if(mx != private$V_gompertz_args$Vmax)
            {
              private$V_gompertz_args$Vmax <- mx
              private$V_gompertz <- NULL
              if(private$V_info[[1]] == 7) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_gompertz_args[[3]]
          if(mn < mx)
          {
            if(mn != private$V_gompertz_args$Vmin)
            {
              private$V_gompertz_args$Vmin <- mn
              private$V_gompertz <- NULL
              if(private$V_info[[1]] == 7) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_gompertz_args)
    },
    #' @description
    #' Set V logistic arguments
    #' @param xi   state at the inflection point
    #' @param vr   rate of change
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xi,vr,Vmax,Vmin)
    set_V_logistic_args = function(xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xi))
      {
        sca <- private$extract_scalar(xi)
        if(!is.null(sca))
        {
          if(sca != private$V_logistic_args$xi)
          {
            private$V_logistic_args$xi <- sca
            private$V_logistic <- NULL
            if(private$V_info[[1]] == 8) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xi not set.")}
      }
      if(!is.null(vr))
      {
        sca <- private$extract_scalar(vr)
        if(!is.null(sca))
        {
          if(sca != private$V_logistic_args$vr)
          {
            private$V_logistic_args$vr <- sca
            private$V_logistic <- NULL
            if(private$V_info[[1]] == 8) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("vr not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_logistic_args$Vmax || mn != private$V_logistic_args$Vmin)
            {
              private$V_logistic_args$Vmax <- mx
              private$V_logistic_args$Vmin <- mn
              private$V_logistic <- NULL
              if(private$V_info[[1]] == 8) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_logistic_args[[4]]
          if(mx > mn)
          {
            if(mx != private$V_logistic_args$Vmax)
            {
              private$V_logistic_args$Vmax <- mx
              private$V_logistic <- NULL
              if(private$V_info[[1]] == 8) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_logistic_args[[3]]
          if(mn < mx)
          {
            if(mn != private$V_logistic_args$Vmin)
            {
              private$V_logistic_args$Vmin <- mn
              private$V_logistic <- NULL
              if(private$V_info[[1]] == 8) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_logistic_args)
    },
    #' @description
    #' Set V transcendental arguments
    #' @param xo   state at the intercept
    #' @param xi   state at the inflection point
    #' @param xm   state at the maximum
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,xi,xm,Vmax,Vmin)
    set_V_transcendental_args = function(xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_transcendental_args$xo)
          {
            private$V_transcendental_args$xo <- sca
            private$V_transcendental <- NULL
            if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(xi))
      {
        sca <- private$extract_scalar(xi)
        if(!is.null(sca))
        {
          if(sca != private$V_transcendental_args$xi)
          {
            private$V_transcendental_args$xi <- sca
            private$V_transcendental <- NULL
            if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xi not set.")}
      }
      if(!is.null(xm))
      {
        sca <- private$extract_scalar(xm)
        if(!is.null(sca))
        {
          if(sca != private$V_transcendental_args$xm)
          {
            private$V_transcendental_args$xm <- sca
            private$V_transcendental <- NULL
            if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xm not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_transcendental_args$Vmax || mn != private$V_transcendental_args$Vmin)
            {
              private$V_transcendental_args$Vmax <- mx
              private$V_transcendental_args$Vmin <- mn
              private$V_transcendental <- NULL
              if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_transcendental_args[[5]]
          if(mx > mn)
          {
            if(mx != private$V_transcendental_args$Vmax)
            {
              private$V_transcendental_args$Vmax <- mx
              private$V_transcendental <- NULL
              if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_transcendental_args[[4]]
          if(mn < mx)
          {
            if(mn != private$V_transcendental_args$Vmin)
            {
              private$V_transcendental_args$Vmin <- mn
              private$V_transcendental <- NULL
              if(private$V_info[[1]] == 9) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_transcendental_args)
    },
    #' @description
    #' Set V yield index arguments
    #' @param xo   state at the intercept
    #' @param xi   state at the inflection point
    #' @param xm   state at the maximum
    #' @param Vmax maximum terminal value
    #' @param Vmin minimum terminal value
    #' @return list(xo,xi,xm,Vmax,Vmin)
    set_V_yieldindex_args = function(xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL)
    {
      if(!is.null(xo))
      {
        sca <- private$extract_scalar(xo)
        if(!is.null(sca))
        {
          if(sca != private$V_yieldindex_args$xo)
          {
            private$V_yieldindex_args$xo <- sca
            private$V_yieldindex <- NULL
            if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xo not set.")}
      }
      if(!is.null(xi))
      {
        sca <- private$extract_scalar(xi)
        if(!is.null(sca))
        {
          if(sca != private$V_yieldindex_args$xi)
          {
            private$V_yieldindex_args$xi <- sca
            private$V_yieldindex <- NULL
            if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xi not set.")}
      }
      if(!is.null(xm))
      {
        sca <- private$extract_scalar(xm)
        if(!is.null(sca))
        {
          if(sca != private$V_yieldindex_args$xm)
          {
            private$V_yieldindex_args$xm <- sca
            private$V_yieldindex <- NULL
            if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
          }
        }
        else { message("xm not set.")}
      }
      if(!is.null(Vmax) && !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) && !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_yieldindex_args$Vmax || mn != private$V_yieldindex_args$Vmin)
            {
              private$V_yieldindex_args$Vmax <- mx
              private$V_yieldindex_args$Vmin <- mn
              private$V_yieldindex <- NULL
              if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message("Vmax <= Vmin.")
            message("Vmax and Vmin not set.")
          }
        }
        else { message("Vmax and Vmin not set.")}
      }
      else if(!is.null(Vmax))
      {
        mx <- private$extract_scalar(Vmax)
        if(!is.null(mx))
        {
          mn <- private$V_yieldindex_args[[5]]
          if(mx > mn)
          {
            if(mx != private$V_yieldindex_args$Vmax)
            {
              private$V_yieldindex_args$Vmax <- mx
              private$V_yieldindex <- NULL
              if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmax <= existing Vmin=",mn))
            message("Vmax not set.")
          }
        }
        else { message("Vmax not set.")}
      }
      else if(!is.null(Vmin))
      {
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mn))
        {
          mx <- private$V_yieldindex_args[[4]]
          if(mn < mx)
          {
            if(mn != private$V_yieldindex_args$Vmin)
            {
              private$V_yieldindex_args$Vmin <- mn
              private$V_yieldindex <- NULL
              if(private$V_info[[1]] == 10) { private$x_stoch_args[3] <- list(V=NULL) }
            }
          }
          else
          {
            message(paste("Vmin >= existing Vmax=",mx))
            message("Vmin not set.")
          }
        }
        else { message("Vmin not set.")}
      }
      return(private$V_yieldindex_args)
    },
    #' @description
    #' Set V arguments by position in list
    #' @param Ix   index of terminal value
    #' @param name name of terminal value
    #' @param v1   first parameter
    #' @param v2   second parameter
    #' @param v3   third parameter
    #' @param v4   fourth parameter
    #' @param v5   fifth parameter
    #' @return list(Ix,name,v1,v2,v3,v4,v5)
    set_V_args = function(Ix=NULL,name=NULL,v1=NULL,v2=NULL,v3=NULL,v4=NULL,v5=NULL)
    {
      OK <- FALSE
      if(!is.null(Ix))
      {
        if(is.numeric(Ix) && is.finite(Ix) && !is.na(Ix))
        {
          Ix <- as.integer(private$extract_scalar(Ix))
          if(Ix > 0 && Ix < 11) { OK <- TRUE }
          else { message("Ix must be from 1 to 10.") }
        }
        else if(is.character(Ix))
        {
          chr <- private$extract_character(Ix)
          Ix <- match(chr,private$V_info$names)
          if(!is.na(Ix)) { OK <- TRUE }
          else { message(paste(sep="",chr," is not recognized. Recognized names are: ",private$V_info$text)) }
        }
        else {  message(paste(sep="","name is not recognized. Recognized names are: ",private$V_info$text)) }
      }
      else
      {
        if(!is.null(name))
        {
          chr <- private$extract_character(name)
          if(!is.null(chr))
          {
            Ix <- match(chr,private$V_info$names)
            if(!is.na(Ix)) { OK <- TRUE }
            else {  message(paste(sep="",chr," is not recognized. Recognized names are: ",private$V_info$text)) }
          }
          else {  message(paste(sep="","name is not recognized. Recognized names are: ",private$V_info$text)) }
        }
      }
      if(OK == FALSE) { Ix <- private$V_info$Ix }
      V_args <- NULL
      if(Ix == 1) { V_args <- self$set_V_linear_args(v1,v2) }
      else if(Ix == 2) { V_args <- self$set_V_degenerate_args(v1,v2,v3) }
      else if(Ix == 3) { V_args <- self$set_V_stepped_args(v1,v2,v3,v4) }
      else if(Ix == 4) { V_args <- self$set_V_kinked_args(v1,v2,v3,v4) }
      else if(Ix == 5) { V_args <- self$set_V_butterfly_args(v1,v2,v3,v4,v5) }
      else if(Ix == 6) { V_args <- self$set_V_mitscherlich_args(v1,v2,v3,v4) }
      else if(Ix == 7) { V_args <- self$set_V_gompertz_args(v1,v2,v3,v4) }
      else if(Ix == 8) { V_args <- self$set_V_logistic_args(v1,v2,v3,v4) }
      else if(Ix == 9) { V_args <- self$set_V_transcendental_args(v1,v2,v3,v4,v5) }
      else if(Ix == 10) { V_args <- self$set_V_yieldindex_args(v1,v2,v3,v4,v5) }

      return(V_args)
    },
    #' @description
    #' Set information for terminal value
    #' @param Ix   index of terminal value
    #' @param name name of terminal value
    #' @return list(Ix,name,names,text)
    set_V_info = function(Ix=NULL,name=NULL)
    {
      if(!is.null(Ix))
      {
        if(is.numeric(Ix))
        {
          sca <- private$extract_scalar(Ix)
          if(!is.null(sca))
          {
            Ix <- as.integer(sca)
            if(Ix > 0 && Ix < 11)
            {
              if(Ix != private$V_info$Ix)
              {
                private$V_info$Ix <- Ix
                private$V_info$name <- private$V_info$names[[Ix]]
                private$x_stoch_args[3] <- list(V=NULL)
                private$OO <- NULL
                private$OOhat <- NULL
                private$kOOhat <- NULL
              }
            }
            else { message("Ix must be from 1 to 10.") }
          }
          else { message("Ix must be from 1 to 10.") }
        }
        else if(is.character(Ix))
        {
          chr <- private$extract_character(Ix)
          Ix <- match(chr,private$V_info$names)
          if(!is.na(Ix))
          {
            if(Ix != private$V_info$Ix)
            {
              private$V_info$Ix <- Ix
              private$V_info$name <- chr
              private$x_stoch_args[3] <- list(V=NULL)
              private$OO <- NULL
              private$OOhat <- NULL
              private$kOOhat <- NULL
            }
          }
          else { message(paste(sep="",chr," is not recognized. Recognized names are: ",private$V_info$text,". For a new name, use set_V_info(name='newname').")) }
        }
      }
      else
      {
        if(!is.null(name))
        {
          chr <- private$extract_character(name)
          Ix <- match(chr,private$V_info$names)
          if(!is.na(Ix))
          {
            if(Ix != private$V_info$Ix)
            {
              private$V_info$Ix <- Ix
              private$x_stoch_args[3] <- list(V=NULL)
              private$OO <- NULL
              private$OOhat <- NULL
              private$kOOhat <- NULL
            }
          }
          private$V_info$name <- chr
        }
      }
      return(private$V_info)
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
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,"FD") }
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
    #' @return list(types,group)
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
        else if(group == 3) # Option, Option Envelope
        {
          least <- 0
          most <- 1
        }
        else # Drift, Terminal Value, Decision Threshold
        {
          group <- 1
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
    #' @param who    object id of sender
    #' @return list(plotit,copyit)
    set_flags = function(plotit=NULL,copyit=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_flags(plotit,copyit,"FD") }
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
    #' @return list(oup_params,x_stoch_args,V_linear_args,V_stepped_args,V_kinked_args,V_butterfly_args,V_mitscherlich_args,V_gompertz_args,V_logistic_args,V_transcendental_args,V_yieldindex_args,plot_info)
    get_all = function()
    {
      all <- list(oup_params=private$oup_params,
        x_stoch_args=private$x_stoch_args,
        V_linear_args=private$V_linear_args,
        V_degenerate_args=private$V_degenerate_args,
        V_stepped_args=private$V_stepped_args,
        V_kinked_args=private$V_kinked_args,
        V_butterfly_args=private$V_butterfly_args,
        V_mitscherlich_args=private$V_mitscherlich_args,
        V_gompertz_args=private$V_gompertz_args,
        V_logistic_args=private$V_logistic_args,
        V_transcendental_args=private$V_transcendental_args,
        V_yieldindex_args=private$V_yieldindex_args,
        V_info=private$V_info,
        plot_info=private$plot_info)

      return(all)
    },
    #' @description
    #' Get OUP parameters
    #' @return list(rho,mu,sigma)
    get_oup_params = function() { return(private$oup_params) },
    #' @description
    #' Get x as a stochastic state and its arguments
    #' @return list(s,x,V,r,phi,theta,skip,ds,dx)
    get_x_stoch_args = function() { return(private$x_stoch_args) },
    #' @description
    #' Get V linear arguments
    #' @return list(xo,vs)
    get_V_linear_args = function() { return(private$V_linear_args) },
    #' @description
    #' Get V degenerate arguments
    #' @return list(xo,Vmax,Vmin)
    get_V_degenerate_args = function() { return(private$V_degenerate_args) },
    #' @description
    #' Get V stepped arguments
    #' @return list(xo,vs,Vmax,Vmin)
    get_V_stepped_args = function() { return(private$V_stepped_args) },
    #' @description
    #' Get V kinked arguments
    #' @return list(xo,vs,Vmax,Vmin)
    get_V_kinked_args = function() { return(private$V_kinked_args) },
    #' @description
    #' Get V butterfly arguments
    #' @return list(xo,xm,vs,Vmax,Vmin)
    get_V_butterfly_args = function() { return(private$V_butterfly_args) },
    #' @description
    #' Get V mitscherlich arguments
    #' @return list(xo,vr,Vmax,Vmin)
    get_V_mitscherlich_args = function() { return(private$V_mitscherlich_args) },
    #' @description
    #' Get V gompertz arguments
    #' @return list(xi,vr,Vmax,Vmin)
    get_V_gompertz_args = function() { return(private$V_gompertz_args) },
    #' @description
    #' Get V logistic arguments
    #' @return list(xi,vr,Vmax,Vmin)
    get_V_logistic_args = function() { return(private$V_logistic_args) },
    #' @description
    #' Get V transcendental arguments
    #' @return list(xo,xi,xm,Vmax,Vmin)
    get_V_transcendental_args = function() { return(private$V_transcendental_args) },
    #' @description
    #' Get V yield index arguments
    #' @return list(xo,xi,xm,Vmax,Vmin)
    get_V_yieldindex_args = function() { return(private$V_yieldindex_args) },
    #' @description
    #' Get V arguments by index or name
    #' @param Ix   index of terminal value
    #' @param name name of terminal value
    #' @return list(xo,xi,xm,vs,vr,Vmax,Vmin)
    get_V_args = function(Ix=NULL,name=NULL)
    {
      OK <- FALSE
      if(!is.null(Ix))
      {
        if(is.numeric(Ix) && is.finite(Ix) && !is.na(Ix))
        {
          Ix <- as.integer(private$extract_scalar(Ix))
          if(Ix > 0 && Ix < 11) { OK <- TRUE }
        }
        else if(is.character(Ix))
        {
          chr <- private$extract_character(Ix)
          Ix <- match(chr,private$V_info$names)
          if(!is.na(Ix)) { OK <- TRUE }
        }
      }
      else
      {
        if(!is.null(name))
        {
          chr <- private$extract_character(name)
          if(!is.null(chr))
          {
            Ix <- match(chr,private$V_info$names)
            if(!is.na(Ix)) { OK <- TRUE }
          }
        }
      }
      if(OK == FALSE) { Ix <- private$V_info$Ix }
      if(Ix == 1) { V_args <- self$get_V_linear_args() }
      else if(Ix == 2) { V_args <- self$get_V_degenerate_args() }
      else if(Ix == 3) { V_args <- self$get_V_stepped_args() }
      else if(Ix == 4) { V_args <- self$get_V_kinked_args() }
      else if(Ix == 5) { V_args <- self$get_V_butterfly_args() }
      else if(Ix == 6) { V_args <- self$get_V_mitscherlich_args() }
      else if(Ix == 7) { V_args <- self$get_V_gompertz_args() }
      else if(Ix == 8) { V_args <- self$get_V_logistic_args() }
      else if(Ix == 9) { V_args <- self$get_V_transcendental_args() }
      else if(Ix == 10) { V_args <- self$get_V_yieldindex_args() }
      return(V_args)
    },
    #' @description
    #' Get information for terminal values
    #' @return list(Ix,name,names,text)
    get_V_info = function() { return(private$V_info) },
    #' @description
    #' Get information for plotting options
    #' @return list(font,file,theme,3D,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
    #' @description
    #' Get current types for plot routines
    #' @return list(types,description)
    get_plot_types = function()
    {
      text <- rbind(c("FiniteDifference groups, types and plots (default type is 0):"),c("  group  types  plots"),c("    1        0  Drift TerminalValue DecisionThreshold"),c("    2     -1,0  Diffusion"),c("    3      0,1  Option OptionEnvelope"))
      return(list(types=private$plot_types,description=text))
    },
    #' @description
    #' Get flags for plotting and copying
    #' @return list(plotit,copyit)
    get_flags = function() { return(private$flags) },
    # public axis method ----
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
      x <- private$x_stoch_args[[2]]
      kOOhat <- private$kOOhat
      if(!is.null(kOOhat)) { k <- kOOhat[[1]] }
      else { k <- mu }
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
      # state
      if(rho > 0) { x <- 8*abs(sigma)/(2*rho)^0.5 }
      else { x <- 8*abs(sigma) }
      if(x < 1) { x <- 1}
      xscale <- 1
      while(x > xscale) { xscale <- 10*xscale }
      x <- round(x/xscale,1)*xscale
      xby <- x/50
      xscale <- 1
      while(abs(k) > xscale) { xscale <- 10*xscale }
      x <- round(k/xscale,2)*xscale
      xfrom <- x-50*xby
      xto <- xfrom+100*xby
      xseq <- seq(from=xfrom,to=xto,by=xby)
      self$set_x_stoch_args(sseq,xseq,NULL,NULL,NULL,NULL,NULL)
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
        x_stoch_args=private$x_stoch_args,
        V_linear_args=private$V_linear_args,
        V_degenerate_args=private$V_degenerate_args,
        V_stepped_args=private$V_stepped_args,
        V_kinked_args=private$V_kinked_args,
        V_butterfly_args=private$V_butterfly_args,
        V_mitscherlich_args=private$V_mitscherlich_args,
        V_gompertz_args=private$V_gompertz_args,
        V_logistic_args=private$V_logistic_args,
        V_transcendental_args=private$V_transcendental_args,
        V_yieldindex_args=private$V_yieldindex_args,
        V_info=private$V_info))
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
        if(private$lists_equal(last_undo_args[[2]],private$x_stoch_args))
        {
          if(private$lists_equal(last_undo_args[[3]],private$V_linear_args))
          {
            if(private$lists_equal(last_undo_args[[4]],private$V_degenerate_args))
            {
              if(private$lists_equal(last_undo_args[[5]],private$V_stepped_args))
              {
                if(private$lists_equal(last_undo_args[[6]],private$V_kinked_args))
                {
                  if(private$lists_equal(last_undo_args[[7]],private$V_butterfly_args))
                  {
                    if(private$lists_equal(last_undo_args[[8]],private$V_mitscherlich_args))
                    {
                      if(private$lists_equal(last_undo_args[[9]],private$V_gompertz_args))
                      {
                        if(private$lists_equal(last_undo_args[[10]],private$V_logistic_args))
                        {
                          if(private$lists_equal(last_undo_args[[11]],private$V_transcendental_args))
                          {
                            if(private$lists_equal(last_undo_args[[12]],private$V_yieldindex_args))
                            {
                              if(last_undo_args[[13]][[1]] == private$V_info[[1]])
                              {
                                if(last_undo_args[[13]][[2]] == private$V_info[[2]])
                                {
                                  not_equal <- FALSE
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
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
            x_stoch_args=private$x_stoch_args,
            V_linear_args=private$V_linear_args,
            V_degenerate_args=private$V_degenerate_args,
            V_stepped_args=private$V_stepped_args,
            V_kinked_args=private$V_kinked_args,
            V_butterfly_args=private$V_butterfly_args,
            V_mitscherlich_args=private$V_mitscherlich_args,
            V_gompertz_args=private$V_gompertz_args,
            V_logistic_args=private$V_logistic_args,
            V_transcendental_args=private$V_transcendental_args,
            V_yieldindex_args=private$V_yieldindex_args,
            V_info=private$V_info)))
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
      x_stoch <- these_undo[[2]]
      V_linear <- these_undo[[3]]
      V_degenerate <- these_undo[[4]]
      V_stepped <- these_undo[[5]]
      V_kinked <- these_undo[[6]]
      V_butterfly <- these_undo[[7]]
      V_mitscherlich <- these_undo[[8]]
      V_gompertz <- these_undo[[9]]
      V_logistic <- these_undo[[10]]
      V_transcendental <- these_undo[[11]]
      V_yieldindex <- these_undo[[12]]
      V_info <- these_undo[[13]]
      self$set_oup_params(oup[[1]],oup[[2]],oup[[3]])
      self$set_x_stoch_args(x_stoch[[1]],x_stoch[[2]],x_stoch[[3]],x_stoch[[4]],x_stoch[[5]],x_stoch[[6]],x_stoch[[7]])
      self$set_V_linear_args(V_linear[[1]],V_linear[[2]])
      self$set_V_degenerate_args(V_degenerate[[1]],V_degenerate[[2]],V_degenerate[[3]])
      self$set_V_stepped_args(V_stepped[[1]],V_stepped[[2]],V_stepped[[3]],V_stepped[[4]])
      self$set_V_kinked_args(V_kinked[[1]],V_kinked[[2]],V_kinked[[3]],V_kinked[[4]])
      self$set_V_butterfly_args(V_butterfly[[1]],V_butterfly[[2]],V_butterfly[[3]],V_butterfly[[4]],V_butterfly[[5]])
      self$set_V_mitscherlich_args(V_mitscherlich[[1]],V_mitscherlich[[2]],V_mitscherlich[[3]],V_mitscherlich[[4]])
      self$set_V_gompertz_args(V_gompertz[[1]],V_gompertz[[2]],V_gompertz[[3]],V_gompertz[[4]])
      self$set_V_logistic_args(V_logistic[[1]],V_logistic[[2]],V_logistic[[3]],V_logistic[[4]])
      self$set_V_transcendental_args(V_transcendental[[1]],V_transcendental[[2]],V_transcendental[[3]],V_transcendental[[4]],V_transcendental[[5]])
      self$set_V_yieldindex_args(V_yieldindex[[1]],V_yieldindex[[2]],V_yieldindex[[3]],V_yieldindex[[4]],V_yieldindex[[5]])
      self$set_V_info(V_info[[1]],V_info[[2]])
      private$undoIx <- undoIx

      return(c(undoIx,n))
    },
    # public calculate methods ----
    #' @description
    #' Calculate, plot and return drifts
    #' @param x       vector of n states
    #' @param rho     rate parameter 0<=<rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param who     object id of caller
    #' @return list(g(1xn))
    Drift = function(x=NULL,rho=NULL,mu=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      x <- private$x_stoch_args[[2]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      drift <- private$g
      if(is.null(drift))
      {
        drift <- RcppOUPFDDrift(x,rho,mu)
        private$g <- drift
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDrift()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Drift",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("x",x),c("g",drift))
          private$CopyToClipboard(clip)
        }
      }
      return(list(g=drift))
    },
    #' @description
    #' Calculate, plot and return diffusions
    #' @param x       vector of n states
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(h2(1xn))
    Diffusion = function(x=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(NULL,NULL,sigma)
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      diffusion <- private$h2
      if(is.null(diffusion))
      {
        n <- length(x)
        diffusion <- RcppOUPFDDiffusion(x,sigma)
        private$h2 <- diffusion
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDiffusion()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Diffusion",rep("",n)),c("sigma",sigma,rep("",n-1)),c("x",x),c("h\u00B2",diffusion))
          private$CopyToClipboard(clip)
        }
      }
      return(list(h2=diffusion))
    },
    #' @description
    #' Create and plot linear terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param vs      slope
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Linear = function(x=NULL,xo=NULL,vs=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_linear_args(xo,vs)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_linear_args[[1]]
      vs <- private$V_linear_args[[2]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_linear
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Linear(x,xo,vs)
        private$V_linear <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 1
      private$V_info$name <- "Linear"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Linear",rep("",n)),c("xo",xo,rep("",n-1)),c("vs",vs,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot degenerate terminal values
    #' @param x       vector of n states
    #' @param xo      state at the spike
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Degenerate = function(x=NULL,xo=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_degenerate_args(xo,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_degenerate_args[[1]]
      Vmax <- private$V_degenerate_args[[2]]
      Vmin <- private$V_degenerate_args[[3]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_degenerate
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Degenerate(x,xo,Vmax,Vmin)
        private$V_degenerate <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 2
      private$V_info$name <- "Degenerate"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Degenerate",rep("",n)),c("xo",xo,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot stepped terminal values
    #' @param x       vector of n states
    #' @param xo      state at the step
    #' @param vs      direction of the step
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Stepped = function(x=NULL,xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_stepped_args(xo,vs,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_stepped_args[[1]]
      vs <- private$V_stepped_args[[2]]
      Vmax <- private$V_stepped_args[[3]]
      Vmin <- private$V_stepped_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_stepped
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Stepped(x,xo,vs,Vmax,Vmin)
        private$V_stepped <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 3
      private$V_info$name <- "Stepped"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Stepped",rep("",n)),c("xo",xo,rep("",n-1)),c("vs",vs,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot kinked terminal values
    #' @param x       vector of n states
    #' @param xo      state at the kink
    #' @param vs      slope
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Kinked = function(x=NULL,xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_kinked_args(xo,vs,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_kinked_args[[1]]
      vs <- private$V_kinked_args[[2]]
      Vmax <- private$V_kinked_args[[3]]
      Vmin <- private$V_kinked_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_kinked
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Kinked(x,xo,vs,Vmax,Vmin)
        private$V_kinked <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 4
      private$V_info$name <- "Kinked"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Kinked",rep("",n)),c("xo",xo,rep("",n-1)),c("vs",vs,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot butterfly terminal values
    #' @param x       vector of n states
    #' @param xo      state at the left wing
    #' @param xm      state at the right wing
    #' @param vs      slope
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Butterfly = function(x=NULL,xo=NULL,xm=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_butterfly_args(xo,xm,vs,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_butterfly_args[[1]]
      xm <- private$V_butterfly_args[[2]]
      vs <- private$V_butterfly_args[[3]]
      Vmax <- private$V_butterfly_args[[4]]
      Vmin <- private$V_butterfly_args[[5]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_butterfly
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Butterfly(x,xo,xm,vs,Vmax,Vmin)
        private$V_butterfly <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 5
      private$V_info$name <- "Butterfly"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Butterfly",rep("",n)),c("xo",xo,rep("",n-1)),c("xm",xm,rep("",n-1)),c("vs",vs,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Mitscherlich terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Mitscherlich = function(x=NULL,xo=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_mitscherlich_args(xo,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_mitscherlich_args[[1]]
      vr <- private$V_mitscherlich_args[[2]]
      Vmax <- private$V_mitscherlich_args[[3]]
      Vmin <- private$V_mitscherlich_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_mitscherlich
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Mitscherlich(x,xo,vr,Vmax,Vmin)
        private$V_mitscherlich <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 6
      private$V_info$name <- "Mitscherlich"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Mitscherlich",rep("",n)),c("xo",xo,rep("",n-1)),c("vr",vr,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Gompertz terminal values
    #' @param x       vector of n states
    #' @param xi      state at the inflection point
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Gompertz = function(x=NULL,xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_gompertz_args(xi,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xi <- private$V_gompertz_args[[1]]
      vr <- private$V_gompertz_args[[2]]
      Vmax <- private$V_gompertz_args[[3]]
      Vmin <- private$V_gompertz_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_gompertz
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Gompertz(x,xi,vr,Vmax,Vmin)
        private$V_gompertz <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 7
      private$V_info$name <- "Gompertz"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Gompertz",rep("",n)),c("xi",xi,rep("",n-1)),c("vr",vr,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Logistic terminal values
    #' @param x       vector of n states
    #' @param xi      state at the inflection point
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Logistic = function(x=NULL,xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_logistic_args(xi,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xi <- private$V_logistic_args[[1]]
      vr <- private$V_logistic_args[[2]]
      Vmax <- private$V_logistic_args[[3]]
      Vmin <- private$V_logistic_args[[4]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_logistic
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Logistic(x,xi,vr,Vmax,Vmin)
        private$V_logistic <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 8
      private$V_info$name <- "Logistic"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Logistic",rep("",n)),c("xi",xi,rep("",n-1)),c("vr",vr,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Transcendental terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param xi      state at the inflection point
    #' @param xm      state at the maximum
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_Transcendental = function(x=NULL,xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_transcendental_args(xo,xi,xm,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_transcendental_args[[1]]
      xi <- private$V_transcendental_args[[2]]
      xm <- private$V_transcendental_args[[3]]
      Vmax <- private$V_transcendental_args[[4]]
      Vmin <- private$V_transcendental_args[[5]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_transcendental
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_Transcendental(x,xo,xi,xm,Vmax,Vmin)
        private$V_transcendental <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 9
      private$V_info$name <- "Transcendental"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Transcendental",rep("",n)),c("xo",xo,rep("",n-1)),c("xi",xi,rep("",n-1)),c("xm",xm,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Yield Index terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param xi      state at the inflection point
    #' @param xm      state at the maximum
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue_YieldIndex = function(x=NULL,xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL,who=NULL)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_yieldindex_args(xo,xi,xm,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_yieldindex_args[[1]]
      xi <- private$V_yieldindex_args[[2]]
      xm <- private$V_yieldindex_args[[3]]
      Vmax <- private$V_yieldindex_args[[4]]
      Vmin <- private$V_yieldindex_args[[5]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      terminalvalues <- private$V_yieldindex
      if(is.null(terminalvalues))
      {
        terminalvalues <- RcppOUPFDTerminalValue_YieldIndex(x,xo,xi,xm,Vmax,Vmin)
        private$V_yieldindex <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 10
      private$V_info$name <- "Yield Index"
      private$OO <- NULL
      private$OOhat <- NULL
      private$kOOhat <- NULL
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Terminal Value",rep("",n)),c("Yield Index",rep("",n)),c("xo",xo,rep("",n-1)),c("xi",xi,rep("",n-1)),c("xm",xm,rep("",n-1)),c("Vmax",Vmax,rep("",n-1)),c("Vmin",Vmin,rep("",n-1)),c("x",x),c("V",terminalvalues))
          private$CopyToClipboard(clip)
        }
      }
      return(list(V=terminalvalues))
    },
    #' @description
    #' Retrieves and plots terminal values
    #' @param Ix      index number or name of terminal values
    #' @param name    name of terminal values
    #' @param who     object id of caller
    #' @return list(V(1xn))
    TerminalValue = function(Ix=NULL,name=NULL,who=NULL)
    {
      # set/get ----
      self$set_V_info(Ix,name)
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      Ix <- private$V_info[[1]]
      Vname <- private$V_info[[2]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # retrieve ----
      if(is.null(V))
      {
        if(Ix == 1) { V <- self$TerminalValue_Linear(who=who)[[1]] }
        else if(Ix == 2) { V <- self$TerminalValue_Degenerate(who=who)[[1]] }
        else if(Ix == 3) { V <- self$TerminalValue_Stepped(who=who)[[1]] }
        else if(Ix == 4) { V <- self$TerminalValue_Kinked(who=who)[[1]] }
        else if(Ix == 5) { V <- self$TerminalValue_Butterfly(who=who)[[1]] }
        else if(Ix == 6) { V <- self$TerminalValue_Mitscherlich(who=who)[[1]] }
        else if(Ix == 7) { V <- self$TerminalValue_Gompertz(who=who)[[1]] }
        else if(Ix == 8) { V <- self$TerminalValue_Logistic(who=who)[[1]] }
        else if(Ix == 9) { V <- self$TerminalValue_Transcendental(who=who)[[1]] }
        else if(Ix == 10) { V <- self$TerminalValue_YieldIndex(who=who)[[1]] }
      }
      # plot or copy ----
      else
      {
        if(is.null(who))
        {
          if(plotit == TRUE) { print(self$PlotTerminalValue()) }
          else if(copyit == TRUE)
          {
            args <- self$get_V_args(Ix)
            arg_names <- names(args)
            args <- unname(args)
            n <- length(x)
            q <- length(args)
            clip <- matrix("",q+5,n+1)
            clip[1,1] <- "Finite Difference"
            clip[2,1] <- "Terminal Value"
            clip[3,1] <- Vname
            for(i in 1:q)
            {
              clip[i+3,1] <- arg_names[[i]]
              clip[i+3,2] <- args[[i]]
            }
            clip[q+4,1] <- "x"
            clip[q+4,2:(n+1)] <- x
            clip[q+5,1] <- "V"
            clip[q+5,2:(n+1)] <- V
            private$CopyToClipboard(clip)
          }
        }
      }
      return(list(V=V))
    },
    #' @description
    #' Calculate and plot option prices
    #' @param s       vector of m times
    #' @param x       vector of n states
    #' @param V       vector of n terminal values
    #' @param r       discount rate 0<=r<inf
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    subdivide time intervals but report at times s 1<=skip<=1000
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(OO(mxn))
    Option = function(s=NULL,x=NULL,V=NULL,r=NULL,theta=NULL,skip=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,x,V,r,NULL,theta,skip)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      c <- private$OO
      if(is.null(c))
      {
        if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
        c <- RcppOUPFDOption(s,x,V,r,theta,skip,rho,mu,sigma)
        private$OO <- c
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotOption()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Options",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("theta",theta,rep("",n-1)),c("ds",ds,rep("",n-1)),c("skip",skip,rep("",n-1)),c("dx",dx,rep("",n-1)),c("\uD835\uDD46(t,x)",x),cbind(s,c))
          private$CopyToClipboard(clip)
        }
      }
      return(list(OO=c))
    },
    #' @description
    #' Calculate and plot the envelope of option prices
    #' @param s       vector of m times
    #' @param x       vector of n states
    #' @param V       vector of n terminal values
    #' @param r       discount rate 0<=r<inf
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    subdivide time intervals but report at times s 1<=skip<=1000
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(OOhat(1xn))
    OptionEnvelope = function(s=NULL,x=NULL,V=NULL,r=NULL,theta=NULL,skip=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,x,V,r,NULL,theta,skip)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      OOhat <- private$OOhat
      shat <- private$shat
      if(is.null(OOhat) || is.null(shat))
      {
        if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
        env <- RcppOUPFDOptionEnvelope(s,x,V,r,theta,skip,rho,mu,sigma)
        OOhat <- env[1,]
        shat <- s[1]-env[2,]
        private$OOhat <- OOhat
        private$shat <- shat
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotOptionEnvelope()) }
        else if(copyit == TRUE)
        {
          n <- length(x)
          clip <- rbind(c("Finite Difference",rep("",n)),c("Option Envelope",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("theta",theta,rep("",n-1)),c("ds",ds,rep("",n-1)),c("skip",skip,rep("",n-1)),c("dx",dx,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat))
          private$CopyToClipboard(clip)
        }
      }
      return(list(OOhat=OOhat,shat=shat))
    },
    #' @description
    #' Calculate and plot the decision threshold
    #' @param s       vector of m times
    #' @param x       vector of n states
    #' @param V       vector of n terminal values
    #' @param r       discount rate 0<=r<inf
    #' @param phi     search direction for exit or entry options
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    subdivide time intervals but report at times s 1<=skip<=1000
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param who     object id of caller
    #' @return list(k,OOhat)
    DecisionThreshold = function(s=NULL,x=NULL,V=NULL,r=NULL,phi=NULL,theta=NULL,skip=NULL,rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(s,x,V,r,phi,theta,skip)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      decision <- private$kOOhat
      if(is.null(decision))
      {
        if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
        OOhat <- private$OOhat #no plot or copy
        if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="FD")[[1]] }
        decision <- RcppOUPFDDecisionThreshold(x,V,OOhat,phi)
        private$kOOhat <- decision
      }
      # plot or copy ----
      if(is.null(who))
      {
        if(plotit == TRUE) { print(self$PlotDecisionThreshold()) }
        else if(copyit == TRUE)
        {
          clip <- rbind(c("Finite Difference",""),c("Decision Threshold",""),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("theta",theta),c("ds",ds),c("skip",skip),c("dx",dx),c("k",decision[1]),c("\u00D4",decision[2]))
          private$CopyToClipboard(clip)
        }
      }
      return(list(k=decision[1],OOhat=decision[2]))
    },
    # public plot methods ----
    #' @description
    #' Plot drifts
    #' @param type  = 0
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDrift = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      x <- private$x_stoch_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      drift <- private$g #protect against recursive call
      if(is.null(drift)) { drift <- self$Drift(who="FD")[[1]] }
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Finite Difference",rep("",n)),c("Drift",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("x",x),c("g",drift))
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_FD_Drift2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),")",esml)
        if(is.null(title)) { title <- "Drift" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "<i>g</i>(<i>x</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside")
      driftline <- list(color=cyn$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Drift2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=drift,name="<i>g</i>(<i>z</i>)",mode="lines",line=driftline) %>%
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
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDiffusion = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      type <- self$set_plot_type(type,2)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      drift <- private$g #no plot or copy
      if(is.null(drift)) { drift <- self$Drift(who="FD")[[1]] }
      diffusion <- private$h2 #protect against recursive call
      if(is.null(diffusion)) { diffusion <- self$Diffusion(who="FD")[[1]]  }
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        diffusion <- diffusion[Ixbeg:Ixend]
        n <- length(x)
      }
      sqrt <- diffusion^0.5
      driftplus <- drift+sqrt
      driftminus <- drift-sqrt
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Finite Difference",rep("",n)),c("Diffusion",rep("",n)),c("sigma",sigma,rep("",n-1)),c("x",x),c("h\u00B2",diffusion))
        private$CopyToClipboard(clip)
      }
      # plot ----
      #2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Diffusion" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
      }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      #OUP_FD_Diffusion2Dg
      if(type < -0.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>g</i>(<i>x</i>)&plusmn;<i>h</i>" }
        lookleft <- list(text=yaxis)
        horz=list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert=list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        diffusionline <- list(color=mgn$d,width=4)
        driftline <- list(color=cyn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Diffusion2Dg")
        legendpos <- list(orientation="h",x=0.5,y=1.05,xanchor="center")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=drift,name="<i>g</i>(<i>x</i>)",mode="lines",line=driftline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=driftplus,name="<i>g</i>(<i>x</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=driftminus,name="<i>g</i>(<i>x</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",showlegend=FALSE,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_FD_Diffusion2D
      else
      {
        if(is.null(yaxis)) { yaxis <- "<i>h</i><sup>2</sup>" }
        lookleft <- list(text=yaxis)
        horz=list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert=list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        diffusionline <- list(color=mgn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Diffusion2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=diffusion,name="<i>h</i><sup>2</sup>",mode="lines",line=diffusionline) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot terminal values
    #' @param type  = 0
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotTerminalValue = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # set/get ----
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
      Ix <- private$V_info[[1]]
      Vname <- private$V_info[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      args <- self$get_V_args(Ix)
      arg_names <- names(args)
      args <- unname(args)
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        q <- length(args)
        clip <- matrix("",q+5,n+1)
        clip[1,1] <- "Finite Difference"
        clip[2,1] <- "Terminal Value"
        clip[3,1] <- Vname
        for(i in 1:q)
        {
          clip[i+3,1] <- arg_names[[i]]
          clip[i+3,2] <- args[[i]]
        }
        clip[q+4,1] <- "x"
        clip[q+4,2:(n+1)] <- x
        clip[q+5,1] <- "V"
        clip[q+5,2:(n+1)] <- V
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_FD_Terminal2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        if(Ix == 1) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>v</i><sub>s</sub>=",format(args[2],digits=4),")",esml) }
        else if(Ix == 2) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>V</i><sub>max</sub>=",format(args[2],digits=4),",<i>V</i><sub>min</sub>=",format(args[3],digits=4),")",esml) }
        else if(Ix == 3) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>v</i><sub>s</sub>=",format(args[2],digits=4),",<i>V</i><sub>max</sub>=",format(args[3],digits=4),",<i>V</i><sub>min</sub>=",format(args[4],digits=4),")",esml) }
        else if(Ix == 4) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>v</i><sub>s</sub>=",format(args[2],digits=4),",<i>V</i><sub>max</sub>=",format(args[3],digits=4),",<i>V</i><sub>min</sub>=",format(args[4],digits=4),")",esml) }
        else if(Ix == 5) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>x</i><sub>m</sub>=",format(args[2],digits=4),",<i>v</i><sub>s</sub>=",format(args[3],digits=4),",<i>V</i><sub>max</sub>=",format(args[4],digits=4),",<i>V</i><sub>min</sub>=",format(args[5],digits=4),")",esml) }
        else if(Ix == 6) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>v</i><sub>r</sub>=",format(args[2],digits=4),",<i>V</i><sub>max</sub>=",format(args[3],digits=4),",<i>V</i><sub>min</sub>=",format(args[4],digits=4),")",esml) }
        else if(Ix == 7) { syms <- paste(sep="",bsml,"(<i>x</i><sub>i</sub>=",format(args[1],digits=4),",<i>v</i><sub>r</sub>=",format(args[2],digits=4),",<i>V</i><sub>max</sub>=",format(args[3],digits=4),",<i>V</i><sub>min</sub>=",format(args[4],digits=4),")",esml) }
        else if(Ix == 8) { syms <- paste(sep="",bsml,"(<i>x</i><sub>i</sub>=",format(args[1],digits=4),",<i>v</i><sub>r</sub>=",format(args[2],digits=4),",<i>V</i><sub>max</sub>=",format(args[3],digits=4),",<i>V</i><sub>min</sub>=",format(args[4],digits=4),")",esml) }
        else if(Ix == 9) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>x</i><sub>i</sub>=",format(args[2],digits=4),",<i>x</i><sub>m</sub>=",format(args[3],digits=4),",<i>V</i><sub>max</sub>=",format(args[4],digits=4),",<i>V</i><sub>min</sub>=",format(args[5],digits=4),")",esml) }
        else if(Ix == 10) { syms <- paste(sep="",bsml,"(<i>x</i><sub>o</sub>=",format(args[1],digits=4),",<i>x</i><sub>i</sub>=",format(args[2],digits=4),",<i>x</i><sub>m</sub>=",format(args[3],digits=4),",<i>V</i><sub>max</sub>=",format(args[4],digits=4),",<i>V</i><sub>min</sub>=",format(args[5],digits=4),")",esml) }
        else { syms <- "" }
        if(is.null(title)) { title <- Vname }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "<i>V</i>(<i>x</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      if(V[n] > V[1]) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
      else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
      Vline <- list(color=gry$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Terminal2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=V,name="<i>V</i>(<i>x</i>)",mode="lines",line=Vline) %>%
        config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot options
    #' @param type = 0, 1, or 'n','p','d' for next, previous, default
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
      type <- self$set_plot_type(type,3)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
      r <- private$x_stoch_args[[4]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      m <- length(s)
      n <- length(x)
      options <- private$OO #protect against recursive call
      if(is.null(options)) { options <- self$Option(who="FD")[[1]] }
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
        V <- V[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Finite Difference",rep("",n)),c("Options",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("theta",theta,rep("",n-1)),c("ds",ds,rep("",n-1)),c("skip",skip,rep("",n-1)),c("dx",dx,rep("",n-1)),c("\uD835\uDD46(t,x)",x),cbind(s,options))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family:Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),bsym,",<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
        if(is.null(title)) { title <- "Options"}
      }
      else if(is.null(title)) { title <- ""}
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_FD_Option2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        if(is.null(yaxis)) { yaxis <- "\uD835\uDD46(<i>s,x</i>|<i>g,h</i><sup>2</sup><i>,V</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(V[n] > V[1]) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        optionline <- list(color=red$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Option2D")
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=x,y=options[i,],name=paste(sep="","\uD835\uDD46(",s[i],"<i>,x</i>)"),mode="lines",line=optionline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+ds
          lineopacity <- lineopacity-0.07
        }
        fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_FD_Option3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis)) { zaxis <- "\uD835\uDD46(<i>s,x</i>|<i>g,h</i><sup>2</sup><i>,V</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            coordinates[i,j] <- paste(sep="","&#x1D546;(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
          }
        }
        if(V[n] > V[1]) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(V[n] == V[1]) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        Vmax <- max(V)
        Vmin <- min(V)
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*x[1]-0.03*x[n],1.03*x[n]-0.03*x[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*s[m]-0.03*s[1],1.03*s[1]-0.03*s[m]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,range=c(1.03*Vmin-0.03*Vmax,1.03*Vmax-0.03*Vmin),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        optionline <- list(color=red$e,width=8)
        gradient <- list(c(0,red$c),c(1,red$c))
        lineopacity <- 1
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions=list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_Option3D")
        fig <- plot_ly()
        ss <- vector("double",n)
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        i <- 1
        for(j in 1:n) { ss[j] <- s[i] }
        fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],name="\uD835\uDD46(<i>s,x</i>)",mode="lines",line=optionline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="O",showlegend=TRUE)
        while(i < m)
        {
          i <- i+ds
          lineopacity <- lineopacity-0.07
          for(j in 1:n) { ss[j] <- s[i] }
          fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],mode="lines",line=optionline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="O",showlegend=FALSE)
        }
        fig <- add_trace(fig,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the option envelope
    #' @param type  = 0, 1, or 'n','p','d' for next, previous, default
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
      type <- self$set_plot_type(type,3)[[1]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
      r <- private$x_stoch_args[[4]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      ylw <- private$plot_colors$ylw
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      copyit <- private$flags[[2]]
      m <- length(s)
      n <- length(x)
      OOhat <- private$OOhat #protect against recursive call
      if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="FD")[[1]] }
      shat <- private$shat
      options <- private$OO #no plot or copy
      if(is.null(options)) { options <- self$Option(who="FD")[[1]] }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        OOhat <- OOhat[Ixbeg:Ixend]
        shat <- shat[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend,drop=FALSE]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Finite Difference",rep("",n)),c("Option Envelope",rep("",n)),c("rho",rho,rep("",n-1)),c("mu",mu,rep("",n-1)),c("sigma",sigma,rep("",n-1)),c("r",r,rep("",n-1)),c("theta",theta,rep("",n-1)),c("ds",ds,rep("",n-1)),c("skip",skip,rep("",n-1)),c("dx",dx,rep("",n-1)),c("x",x),c("\u00D4",OOhat),c("\u015D",shat))
        private$CopyToClipboard(clip)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),bsym,",<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,r,",",bsym,"<i>q</i>=",esym,theta,",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
        if(is.null(title)) { title <- "Option Envelope"}
      }
      else if(is.null(title)) { title <- ""}
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_FD_Envelope2D
      if(type < 0.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
        if(is.null(yaxis)) { yaxis <- "\u00D4(<i>x</i>|<i>g,h</i><sup>2</sup>,<i>V</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(V[n] > V[1]) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
        else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
        OOhatline <- list(color=red$d,width=4)
        terminalline <- list(color=gry$c,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Envelope2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=V,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=OOhat,name="\u00D4(<i>x</i>)",mode="lines",line=OOhatline,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_FD_Envelope3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis)) { zaxis <- "\u00D4(<i>x</i>|<i>g,h</i><sup>2</sup>,<i>V</i>)" }
        OOhold <- vector("double",n)
        OOexercise <- vector("double",n)
        coordinatesenv <- vector("character",n)
        coordinates <- matrix("",m,n)
        for(j in 1:n)
        {
          if(shat[j] == s[1])
          {
            OOhold[j] <- NA
            OOexercise[j] <- OOhat[j]
          }
          else
          {
            OOhold[j] <- OOhat[j]
            OOexercise[j] <- NA
          }
          coordinatesenv[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(OOhat[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          for(i in 1:m) { coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
        }
        OOholdmesh <- MeshCurtainSmooth(x,shat,OOhold,rep(0,n))
        OOexercisemesh <- MeshCurtainSmooth(x,shat,OOexercise,rep(0,n))
        if(V[n] > V[1]) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(V[n] == V[1]) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        Vmax <- max(V)
        Vmin <- min(V)
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*x[1]-0.03*x[n],1.03*x[n]-0.03*x[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(1.03*s[m]-0.03*s[1],1.03*s[1]-0.03*s[m]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,range=c(1.03*Vmin-0.03*Vmax,1.03*Vmax-0.03*Vmin),tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        OOholdline <- list(color=ylw$e,width=10)
        OOexerciseline <- list(color=ylw$e,width=8)
        OOhatline <- list(dash="dash",color=ylw$e,width=6)
        gradientOO <- list(c(0,ylw$d),c(1,ylw$d))
        gradient <- list(c(0,red$c),c(1,red$c))
        legendpos <- list(orientation="h",x=0.5,y=0.92,xanchor="center")
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_OptionEnvelope3D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=OOhat,name="\u00D4(<i>x</i>)",mode="lines",line=OOhatline,hoverinfo="text",text=coordinatesenv,legendgroup="OOhat",showlegend=TRUE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=OOhold,mode="lines",line=OOholdline,hoverinfo="text",text=coordinatesenv,legendgroup="OOhat",showlegend=FALSE) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=OOexercise,mode="lines",line=OOexerciseline,hoverinfo="text",text=coordinatesenv,legendgroup="OOhat",showlegend=FALSE) %>%
          add_trace(.,type="mesh3d",x=OOholdmesh$xvertex,y=OOholdmesh$yvertex,z=OOholdmesh$zvertex,i=OOholdmesh$ivertex,j=OOholdmesh$jvertex,k=OOholdmesh$kvertex,intensity=OOholdmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientOO,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="OOhat",showlegend=FALSE) %>%
          add_trace(.,type="mesh3d",x=OOexercisemesh$xvertex,y=OOexercisemesh$yvertex,z=OOexercisemesh$zvertex,i=OOexercisemesh$ivertex,j=OOexercisemesh$jvertex,k=OOexercisemesh$kvertex,intensity=OOexercisemesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientOO,reversescale=reverse,opacity=0.5,hoverinfo="skip",legendgroup="OOhat",showlegend=FALSE) %>%
          add_trace(.,type="surface",x=x,y=s,z=options,name="\uD835\uDD46(<i>s,x</i>)",showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,visible="legendonly",showlegend=TRUE) %>%
          config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_3D,displaylogo=FALSE) %>%
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
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(who="FD")[[1]] }
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      m <- length(s)
      n <- length(x)
      OOhat <- private$OOhat #no plot or copy
      if(is.null(OOhat)) { OOhat <- self$OptionEnvelope(who="FD")[[1]] }
      decision <- private$kOOhat #protect against recursive call
      if(is.null(decision))
      {
        decision <- self$DecisionThreshold(who="FD")
        k <- decision[[1]]
        OOhat <- decision[[2]]
      }
      else
      {
        k <- decision[1]
        OOhat <- decision[2]
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        OOhat <- OOhat[Ixbeg:Ixend]
        n <- length(x)
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Finite Difference",""),c("Decision Threshold",""),c("rho",rho),c("mu",mu),c("sigma",sigma),c("r",r),c("phi",phi),c("theta",theta),c("ds",ds),c("skip",skip),c("dx",dx),c("k",k),c("\u00D4",OOhat))
        private$CopyToClipboard(clip)
      }
      # plot ----
      # OUP_FD_Decision2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),bsym,",<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,r,",",bsym,"<i>f</i>=",esym,phi,",",bsym,"<i>q</i>=",esym,theta,",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
        if(is.null(title)) { title <- "Decision Threshold"}
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>x</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- ""}
        if(is.null(xaxis)) { xaxis <- "<i>x</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "\u00D4(<i>x</i>|<i>g,h</i><sup>2</sup>,<i>V</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      if(phi > 0 || (phi == 0 && V[n] > V[1])) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
      else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
      OOhatline <- list(color=red$d,width=4)
      terminalline <- list(color=gry$c,width=4)
      OOline <- list(dash="dot",color=red$d,width=4)
      kline <- list(dash="dot",color=red$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Decision2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=V,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=x,y=OOhat,name="\u00D4(<i>x</i>)",mode="lines",line=OOhatline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(k,k),y=c(0,OOhat),name="<i>k</i>",mode="lines",line=kline,hoverinfo="x+y")
      if(phi > 0 || (phi == 0 && V[n] > V[1]))
      {
        fig <- add_trace(fig,type="scatter",x=c(x[n],k),y=c(OOhat,OOhat),mode="lines",line=OOline,hoverinfo="x+y")
        kOOhat <- list(x=k,y=OOhat,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(OOhat,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",showarrow=FALSE)
      }
      else
      {
        fig <- add_trace(fig,type="scatter",x=c(x[1],k),y=c(OOhat,OOhat),mode="lines",line=OOline,hoverinfo="x+y")
        kOOhat <- list(x=k,y=OOhat,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(OOhat,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",showarrow=FALSE)
      }
      fig <- config(fig,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,annotations=kOOhat,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    OUP = NULL,
    # private attributes ----
    oup_params = NULL,
    x_stoch_args = NULL,
    V_linear_args = NULL,
    V_degenerate_args = NULL,
    V_stepped_args = NULL,
    V_kinked_args = NULL,
    V_butterfly_args = NULL,
    V_mitscherlich_args = NULL,
    V_gompertz_args = NULL,
    V_logistic_args = NULL,
    V_transcendental_args = NULL,
    V_yieldindex_args = NULL,
    undo_args = NULL,
    plot_info = NULL,
    flags = NULL,
    # private intermediate fields ----
    g = NULL,
    h2 = NULL,
    V_linear = NULL,
    V_degenerate = NULL,
    V_stepped = NULL,
    V_kinked = NULL,
    V_butterfly = NULL,
    V_mitscherlich = NULL,
    V_gompertz = NULL,
    V_logistic = NULL,
    V_transcendental = NULL,
    V_yieldindex = NULL,
    V_info = NULL,
    shat = NULL,
    # private output fields ----
    OO = NULL,
    OOhat = NULL,
    kOOhat = NULL,
    # private globals ----
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
    # private clipboard methods ----
    CopyToClipboard = function(clip)
    {
      if(!is.null(private$OUP)) { OUP$CopyToClipboard(clip) }
      else if(interactive() && clipr_available()) { write_clip(clip,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    }
  )
)

