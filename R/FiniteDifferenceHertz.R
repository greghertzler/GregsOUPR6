library(R6)
library(plotly)
library(stringr)

# roxygen ----
#' R6 Class implementing a Finite Difference method.
#'
#' @description
#' A Finite Difference method for the Ornstein-Uhlenbeck Process--without boundary
#'  conditions on the state but with arbitrary terminal values--for the pricing of
#'  index insurance, perpetual options and sequential options. From the option
#'  prices, an option envelope and a decision threshold can be calculated.  Drift,
#'  diffusion and several possible terminal values are pre-programmed.  User-defined
#'  terminal values can also be entered.
#'
#' @details # Formulas and methods:
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
#' @details # Using the formulas and methods:
#' Demonstration scripts are in files in the 'demo' directory. Identify a
#'  formula and in the console type:
#'
#'       demo(FD_FormulaName), or
#'       demo(FD_PlotFormulaName).
#'
#' @details # Discussion:
#' The Finite Difference Method can solve problems with no analytical solutions.
#'  As a numerical procedure, it will have errors. Large dif fusions will cause
#'  larger errors. To manage the errors, the domain for x should be as large as
#'  practical. The increment between times is ds and skip subdivides time into
#'  small intervals. Only calculations for ds are reported.  The increment
#'  between states is dx and ds/skip should be at most dx/100. Increasing skip
#'  or dx may reduce errors up to a point. If possible, the Finite Difference
#'  Method should be calibrated against an analytical solution for a simple
#'  problem before solving more complicated problems.
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
#'  list must be stripped of its names, the names Ohat for the option envelope
#'  and V for the terminal value, in the above example. The functions try to
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
#'  to eliminate names if you know the position in the list. Both 'Ohat' and 'V'
#'  are first in their lists:
#'
#'       FD$TerminalValue_Kinked(xo=0,vs=-1)
#'       exit <- FD$OptionEnvelope()[[1]]
#'       V <- exit+FD$TerminalValue_Kinked(xo=5,vs=1)[[1]]
#'       FD$DecisionThreshold(V=V,phi=1)
#'
#' Sorry.
#'
#' Examples in R6 don't work the same way as other R modules.  There is only
#'  one example for an R6 object, not one for each function in the object.
#'  To run examples, the devtools::run_examples() works, but the R command
#'  example("FiniteDifference") doesn't.  You can copy commands to the clipboard,
#'  paste into the console and press Enter. Examples in this help and a simple
#'  example at the bottom can be run in this way.
#'
#' A better alternative is demo(), but this works sometimes, sometimes not.
#'  If you don't see the demos for this package, go to Files and demo.  You will
#'  see the demo names and then type something like demo(FD_DecisionThreshold).

# class ----
FiniteDifference <- R6::R6Class("FiniteDifference",
  portable = FALSE,
  cloneable = FALSE,
  # portable = TRUE,
  # cloneable = TRUE,
#' @import plotly
#' @import stringr
#' @export
  #public members ----
  public = list(
    # constructor ----
    #' @description
    #' Create a new FiniteDifference object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new FiniteDifference object
    #' @examples
    #'   FD <- FiniteDifference$new()
    #'   FD$TerminalValue_Mitscherlich(plotit=FALSE)
    #'   FD$DecisionThreshold()
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
      self$default_save()
      # plot info ----
      plottype <- list(type=3)
      plotfont <- list(family="Cambria",size=14)
      plotfile <- list(format="png",width=640,height=480)
      plottheme <- list(name="light",opaque=1.0)
      if(rstudioapi::isAvailable())
      {
        if(rstudioapi::getThemeInfo()$dark) { plottheme$name <- "dark"}
      }
      plot3D <- list(walls=TRUE,floor=TRUE)
      private$plot_info <- list(plottype=plottype,plotfont=plotfont,plotfile=plotfile,plottheme=plottheme,plot3D=plot3D,plotlabels=TRUE)
      private$plot_colors <- private$rainbow(plottheme$name,plottheme$opaque)
    },
    # public set methods ----
    #' @description
    #' Set OUP parameters
    #' @param rho   rate parameter 0<=rho<inf
    #' @param mu    location parameter -inf<mu<inf
    #' @param sigma scale parameter -inf<sigma<inf
    #' @param who   identifier for object sending the parameters
    #' @return list(rho,mu,sigma)
    set_oup_params = function(rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"FD")}
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
    #' @param skip  divide the time interval and report every skip result
    #' @param who   identifier for object sending the parameters
    #' @return list(s,x,V,r,phi,theta,skip,ds,dx)
    set_x_stoch_args = function(s=NULL,x=NULL,V=NULL,r=NULL,phi=NULL,theta=NULL,skip=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,r,phi,"FD") }
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
              while(i < m & ok == TRUE)
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
                  private$O <- NULL
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
              while(i < nx & cancel == FALSE)
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
            private$O <- NULL
            private$Ohat <- NULL
            private$kOhat <- NULL
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_degenerate_args$Vmax | mn != private$V_degenerate_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_stepped_args$Vmax | mn != private$V_stepped_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_kinked_args$Vmax | mn != private$V_kinked_args$Vmin)
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
      if(!is.null(xo) & !is.null(xm))
      {
        scao <- private$extract_scalar(xo)
        scam <- private$extract_scalar(xm)
        if(!is.null(scao) & !is.null(scam))
        {
          if(scao <= scam)
          {
            if(scao != private$V_butterfly_args$xo | scam != private$V_butterfly_args$xm)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_butterfly_args$Vmax | mn != private$V_butterfly_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_mitscherlich_args$Vmax | mn != private$V_mitscherlich_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_gompertz_args$Vmax | mn != private$V_gompertz_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_logistic_args$Vmax | mn != private$V_logistic_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_transcendental_args$Vmax | mn != private$V_transcendental_args$Vmin)
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
      if(!is.null(Vmax) & !is.null(Vmin))
      {
        mx <- private$extract_scalar(Vmax)
        mn <- private$extract_scalar(Vmin)
        if(!is.null(mx) & !is.null(mn))
        {
          if(mx > mn)
          {
            if(mx != private$V_yieldindex_args$Vmax | mn != private$V_yieldindex_args$Vmin)
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
        if(is.numeric(Ix) & is.finite(Ix) & !is.na(Ix))
        {
          Ix <- as.integer(private$extract_scalar(Ix))
          if(Ix > 0 & Ix < 11) { OK <- TRUE }
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
            if(Ix > 0 & Ix < 11)
            {
              if(Ix != private$V_info$Ix)
              {
                private$V_info$Ix <- Ix
                private$V_info$name <- private$V_info$names[[Ix]]
                private$x_stoch_args[3] <- list(V=NULL)
                private$O <- NULL
                private$Ohat <- NULL
                private$kOhat <- NULL
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
              private$O <- NULL
              private$Ohat <- NULL
              private$kOhat <- NULL
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
          if(!is.null(chr))
          {
            Ix <- match(chr,private$V_info$names)
            if(!is.na(Ix))
            {
              if(Ix != private$V_info$Ix)
              {
                private$V_info$Ix <- Ix
                private$V_info$name <- chr
                private$x_stoch_args[3] <- list(V=NULL)
                private$O <- NULL
                private$Ohat <- NULL
                private$kOhat <- NULL
              }
            }
            else { private$V_info$name <- chr  }
          }
          else { message("'name' is not recognized. Names are: ",private$V_info$text) }
        }
      }
      return(private$V_info)
    },
    #' @description
    #' Set information for plotting
    #' @param type       = 1 and 2 for 2D, 3 and 4 for 3D
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
    #' @param who        identifier for object sending the parameters
    #' @return list(type,font,file,theme,3D)
    set_plot_info = function(type=NULL,fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_plot_info(type,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,"FD") }
      if(!is.null(type))
      {
        sca <- private$extract_scalar(type)
        if(!is.null(sca)) { private$plot_info$plottype$type <- sca }
        else { message("type not set.") }
      }
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
          if(chr == "png" | chr == "svg") { private$plot_info$plotfile$format <- chr }
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
      if(!is.null(theme) | !is.null(opaque))
      {
        if(!is.null(theme))
        {
          chr <- private$extract_character(theme)
          if(!is.null(chr))
          {
            if(chr == "light" | chr == "dark") { private$plot_info$plottheme$name <- chr }
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
    # public get methods ----
    #' @description
    #' Get all arguments and parameters
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
        if(is.numeric(Ix) & is.finite(Ix) & !is.na(Ix))
        {
          Ix <- as.integer(private$extract_scalar(Ix))
          if(Ix > 0 & Ix < 11) { OK <- TRUE }
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
    #' @return list(Ix,names)
    get_V_info = function() { return(private$V_info) },
    #' @description
    #' Get information for plotting options
    #' @return list(type,font,file,theme,3D,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
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
      kOhat <- private$kOhat
      if(!is.null(kOhat)) { k <- kOhat[[1]] }
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
      if(!private$vecs_equal(sseq,private$x_stoch_args$s))
      {
        private$x_stoch_args$s <- sseq
        if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(sseq,NULL,NULL,NULL,"FD") }
        private$x_stoch_args$ds <- -sby
        private$O <- NULL
        private$Ohat <- NULL
        private$kOhat <- NULL
      }
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
      if(!private$vecs_equal(xseq,private$x_stoch_args$x))
      {
        private$x_stoch_args$x <- xseq
        if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(NULL,xseq,NULL,NULL,"FD") }
        private$x_stoch_args[3] <- list(V=NULL)
        private$x_stoch_args$dx <- xby
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
        private$O <- NULL
        private$Ohat <- NULL
        private$kOhat <- NULL
      }
      return(NULL)
    },
    # public default methods ----
    #' @description
    #' Get default
    #' @return list(oup_params,x_stoch_args,V_linear_args,V_degenerate_args,V_stepped_args,V_kinked_args,V_butterfly_args,V_mitscherlich_args,V_gompertz_args,V_logistic_args,V_transcendental_args,V_yieldindex_args,V_info)
    get_default = function() { return(private$default) },
    #' @description
    #' Save to_default->current arguments are stored as defaults
    #' @return NULL
    default_save = function()
    {
      private$default <- list(oup_params=private$oup_params,
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
        V_info=private$V_info)

      return(NULL)
    },
    #' @description
    #' Read from default<-defaults replace current arguments
    #' @return NULL
    default_read = function()
    {
      private$oup_params <- private$default[[1]]
      rho <- private$oup_params$rho
      mu <- private$oup_params$mu
      sigma <- private$oup_params$sigma
      if(!is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"FD") }
      private$x_stoch_args <- private$default[[2]]
      s <- private$x_stoch_args$s
      x <- private$x_stoch_args$x
      r <- private$x_stoch_args$r
      phi <- private$x_stoch_args$phi
      if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,r,phi,"FD") }
      private$V_linear_args <- private$default[[3]]
      private$V_degenerate_args <- private$default[[4]]
      private$V_stepped_args <- private$default[[5]]
      private$V_kinked_args <- private$default[[6]]
      private$V_butterfly_args <- private$default[[7]]
      private$V_mitscherlich_args <- private$default[[8]]
      private$V_gompertz_args <- private$default[[9]]
      private$V_logistic_args <- private$default[[10]]
      private$V_transcendental_args <- private$default[[11]]
      private$V_yieldindex_args <- private$default[[12]]
      private$V_info <- private$default[[13]]
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
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL

      return(NULL)
    },
    # public calculate methods ----
    #' @description
    #' Calculate, plot and return drifts
    #' @param x       vector of n states
    #' @param rho     rate parameter 0<=<rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(g(1xn))
    Drift = function(x=NULL,rho=NULL,mu=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      x <- private$x_stoch_args[[2]]
      # calculate ----
      drift <- private$g
      if(is.null(drift))
      {
        drift <- -rho*(x-mu)
        private$g <- drift
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDrift(NULL)) }

      return(list(g=drift))
    },
    #' @description
    #' Calculate, plot and return diffusions
    #' @param x       vector of n states
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(h2(1xn))
    Diffusion = function(x=NULL,sigma=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(NULL,NULL,sigma)
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      # calculate ----
      diffusion <- private$h2
      if(is.null(diffusion))
      {
        n <- length(x)
        diffusion <- rep(sigma^2,n)
        private$h2 <- diffusion
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDiffusion(NULL,NULL)) }

      return(list(h2=diffusion))
    },
    #' @description
    #' Create and plot linear terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param vs      slope
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Linear = function(x=NULL,xo=NULL,vs=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_linear_args(xo,vs)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_linear_args[[1]]
      vs <- private$V_linear_args[[2]]
      # calculate ----
      terminalvalues <- private$V_linear
      if(is.null(terminalvalues))
      {
        terminalvalues <- vs*(x-xo)
        private$V_linear <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 1
      private$V_info$name <- "Linear"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot degenerate terminal values
    #' @param x       vector of n states
    #' @param xo      state at the spike
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Degenerate = function(x=NULL,xo=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_degenerate_args(xo,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      dx <- private$x_stoch_args[[9]]
      xo <- private$V_degenerate_args[[1]]
      Vmax <- private$V_degenerate_args[[2]]
      Vmin <- private$V_degenerate_args[[3]]
      # calculate ----
      terminalvalues <- private$V_degenerate
      if(is.null(terminalvalues))
      {
        n <- length(x)
        terminalvalues <- vector("double",n)
        j <- 0
        while(j < n)
        {
          j <- j+1
          if(x[j] >= xo-0.5*dx & x[j] < xo+0.5*dx) { terminalvalues[j] <- Vmax }
          else { terminalvalues[j] <- Vmin }
        }
        private$V_degenerate <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 2
      private$V_info$name <- "Degenerate"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot stepped terminal values
    #' @param x       vector of n states
    #' @param xo      state at the step
    #' @param vs      direction of the step
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Stepped = function(x=NULL,xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_stepped_args(xo,vs,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_stepped_args[[1]]
      vs <- private$V_stepped_args[[2]]
      Vmax <- private$V_stepped_args[[3]]
      Vmin <- private$V_stepped_args[[4]]
      # calculate ----
      terminalvalues <- private$V_stepped
      if(is.null(terminalvalues))
      {
        n <- length(x)
        terminalvalues <- vector("double",n)
        if(vs > 0)
        {
          j <- 0
          while(j < n)
          {
            j <- j+1
            if(x[j] < xo) { terminalvalues[j] <- Vmin }
            else if(x[j] == xo) { terminalvalues[j] <- 0.5*(Vmax+Vmin) }
            else { terminalvalues[j] <- Vmax }
          }
        }
        else
        {
          j <- 0
          while(j < n)
          {
            j <- j+1
            if(x[j] < xo) { terminalvalues[j] <- Vmax }
            else if(x[j] == xo) { terminalvalues[j] <- 0.5*(Vmax+Vmin) }
            else { terminalvalues[j] <- Vmin }
          }
        }
        private$V_stepped <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 3
      private$V_info$name <- "Stepped"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot kinked terminal values
    #' @param x       vector of n states
    #' @param xo      state at the kink
    #' @param vs      slope
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or
    #' @return list(V(1xn))
    TerminalValue_Kinked = function(x=NULL,xo=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_kinked_args(xo,vs,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_kinked_args[[1]]
      vs <- private$V_kinked_args[[2]]
      Vmax <- private$V_kinked_args[[3]]
      Vmin <- private$V_kinked_args[[4]]
      # calculate ----
      terminalvalues <- private$V_kinked
      if(is.null(terminalvalues))
      {
        terminalvalues <- vs*(x-xo)
        n <- length(x)
        j <- 0
        while(j < n)
        {
          j <- j+1
          if(terminalvalues[j] > Vmax) { terminalvalues[j] <- Vmax }
          else if(terminalvalues[j] < Vmin) { terminalvalues[j] <- Vmin }
        }
        private$V_kinked <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 4
      private$V_info$name <- "Kinked"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

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
    #' @param plotit  TRUE or
    #' @return list(V(1xn))
    TerminalValue_Butterfly = function(x=NULL,xo=NULL,xm=NULL,vs=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
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
      # calculate ----
      terminalvalues <- private$V_butterfly
      if(is.null(terminalvalues))
      {
        slope <- abs(vs)
        n <- length(x)
        j <- 0
        while(j < n)
        {
          j <- j+1
          left <- -slope*(x[j]-xo)
          right <- slope*(x[j]-xm)
          if(left > right) { terminalvalues[j] <- left }
          else { terminalvalues[j] <- right }
          if(terminalvalues[j] > Vmax) { terminalvalues[j] <- Vmax }
          else if(terminalvalues[j] < Vmin) { terminalvalues[j] <- Vmin }
        }
        private$V_butterfly <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 5
      private$V_info$name <- "Butterfly"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Mitscherlich terminal values
    #' @param x       vector of n states
    #' @param xo      state at the intercept
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Mitscherlich = function(x=NULL,xo=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_mitscherlich_args(xo,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xo <- private$V_mitscherlich_args[[1]]
      vr <- private$V_mitscherlich_args[[2]]
      Vmax <- private$V_mitscherlich_args[[3]]
      Vmin <- private$V_mitscherlich_args[[4]]
      # calculate ----
      terminalvalues <- private$V_mitscherlich
      if(is.null(terminalvalues))
      {
        n <- length(x)
        terminalvalues <- vector("double",n)
        j <- 0
        while(j < n)
        {
          j <- j+1
          if(vr*(x[j]-xo) > 0) { terminalvalues[j] <- Vmin+(Vmax-Vmin)*(1-exp(-vr*(x[j]-xo))) }
          else { terminalvalues[j] <- Vmin }
        }
        private$V_mitscherlich <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 6
      private$V_info$name <- "Mitscherlich"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Gompertz terminal values
    #' @param x       vector of n states
    #' @param xi      state at the inflection point
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Gompertz = function(x=NULL,xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_gompertz_args(xi,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xi <- private$V_gompertz_args[[1]]
      vr <- private$V_gompertz_args[[2]]
      Vmax <- private$V_gompertz_args[[3]]
      Vmin <- private$V_gompertz_args[[4]]
      # calculate ----
      terminalvalues <- private$V_gompertz
      if(is.null(terminalvalues))
      {
        if(vr > 0) { terminalvalues <- Vmin+(Vmax-Vmin)*exp(-exp(-vr*(x-xi))) }
        else { terminalvalues <- Vmax-(Vmax-Vmin)*exp(-exp(vr*(x-xi))) }
        private$V_gompertz <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 7
      private$V_info$name <- "Gompertz"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Create and plot Logistic terminal values
    #' @param x       vector of n states
    #' @param xi      state at the inflection point
    #' @param vr      rate of change
    #' @param Vmax    maximum terminal value
    #' @param Vmin    minimum terminal value
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Logistic = function(x=NULL,xi=NULL,vr=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_x_stoch_args(NULL,x,NULL,NULL,NULL,NULL,NULL)
      self$set_V_logistic_args(xi,vr,Vmax,Vmin)
      x <- private$x_stoch_args[[2]]
      xi <- private$V_logistic_args[[1]]
      vr <- private$V_logistic_args[[2]]
      Vmax <- private$V_logistic_args[[3]]
      Vmin <- private$V_logistic_args[[4]]
      # calculate ----
      terminalvalues <- private$V_logistic
      if(is.null(terminalvalues))
      {
        if(vr > 0) { terminalvalues <- Vmin+(Vmax-Vmin)/(1+exp(-vr*(x-xi))) }
        else { terminalvalues <- Vmax-(Vmax-Vmin)/(1+exp(vr*(x-xi))) }
        private$V_logistic <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 8
      private$V_info$name <- "Logistic"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

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
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_Transcendental = function(x=NULL,xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
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
      # calculate ----
      terminalvalues <- private$V_transcendental
      if(is.null(terminalvalues))
      {
        n <- length(x)
        terminalvalues <- vector("double",n)
        if(abs(xm-xi) < 0.0000000001 | Vmax-Vmin < 0.0000000001 | (xm-xi)*(xi-xo) < 0)
        {
          for(j in 1:n) { terminalvalues[j] <- Vmin }
        }
        else
        {
          b <- ((xm-xo)/(xm-xi))^2
          c <- (xm-xo)/(xm-xi)^2
          j <- 0
          while(j < n)
          {
            j <- j+1
            if((x[j]-xo)*(xm-xo) > 0)
            {
              lnv <- log(Vmax-Vmin)+b*log((x[j]-xo)/(xm-xo))-c*(x[j]-xm)
              terminalvalues[j] <- Vmin+exp(lnv)
            }
            else { terminalvalues[j] <- Vmin }
          }
        }
        private$V_transcendental <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 9
      private$V_info$name <- "Transcendental"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

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
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue_YieldIndex = function(x=NULL,xo=NULL,xi=NULL,xm=NULL,Vmax=NULL,Vmin=NULL,plotit=TRUE)
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
      # calculate ----
      terminalvalues <- private$V_yieldindex
      if(is.null(terminalvalues))
      {
        n <- length(x)
        terminalvalues <- vector("double",n)
        if(abs(xm-xi) < 0.0000000001 | Vmax-Vmin < 0.0000000001 | (xm-xi)*(xi-xo) < 0)
        {
          for(j in 1:n) { terminalvalues[j] <- Vmax }
        }
        else
        {
          b <- ((xm-xo)/(xm-xi))^2
          c <- (xm-xo)/(xm-xi)^2
          j <- 0
          while(j < n)
          {
            j <- j+1
            if((x[j]-xo)*(xm-xo) > 0)
            {
              lnv <- log(Vmax-Vmin)+b*log((x[j]-xo)/(xm-xo))-c*(x[j]-xm)
              terminalvalues[j] <- Vmax-exp(lnv)
              if(terminalvalues[j] < 0) { terminalvalues[j] <- 0 }
            }
            else { terminalvalues[j] <- Vmax }
          }
        }
        private$V_yieldindex <- terminalvalues
      }
      private$x_stoch_args$V <- terminalvalues
      private$V_info$Ix <- 10
      private$V_info$name <- "Yield Index"
      private$O <- NULL
      private$Ohat <- NULL
      private$kOhat <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }

      return(list(V=terminalvalues))
    },
    #' @description
    #' Retrieves and plots terminal values
    #' @param Ix      index number or name of terminal values
    #' @param name    name of terminal values
    #' @param plotit  TRUE or FALSE
    #' @return list(V(1xn))
    TerminalValue = function(Ix=NULL,name=NULL,plotit=TRUE)
    {
      # set/get ----
      self$set_V_info(Ix,name)
      V <- private$x_stoch_args[[3]]
      # retrieve and plot ----
      if(is.null(V))
      {
        Ix <- private$V_info[[1]]
        if(Ix == 1) { V <- self$TerminalValue_Linear(plotit=plotit)[[1]] }
        else if(Ix == 2) { V <- self$TerminalValue_Degenerate(plotit=plotit)[[1]] }
        else if(Ix == 3) { V <- self$TerminalValue_Stepped(plotit=plotit)[[1]] }
        else if(Ix == 4) { V <- self$TerminalValue_Kinked(plotit=plotit)[[1]] }
        else if(Ix == 5) { V <- self$TerminalValue_Butterfly(plotit=plotit)[[1]] }
        else if(Ix == 6) { V <- self$TerminalValue_Mitscherlich(plotit=plotit)[[1]] }
        else if(Ix == 7) { V <- self$TerminalValue_Gompertz(plotit=plotit)[[1]] }
        else if(Ix == 8) { V <- self$TerminalValue_Logistic(plotit=plotit)[[1]] }
        else if(Ix == 9) { V <- self$TerminalValue_Transcendental(plotit=plotit)[[1]] }
        else if(Ix == 10) { V <- self$TerminalValue_YieldIndex(plotit=plotit)[[1]] }
      }
      else
      {
        if(plotit == TRUE) { print(self$PlotTerminalValue(NULL)) }
      }
      return(list(V=V))
    },
    #' @description
    #' Calculate and plot option prices
    #' @param s       vector of m times
    #' @param x       vector of n states
    #' @param V       vector of n terminal values
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param r       discount rate 0<=r<inf
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    divide the time interval and report every skip result
    #' @param plotit  TRUE or FALSE
    #' @return list(O(mxn))
    Option = function(s=NULL,x=NULL,V=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,theta=NULL,skip=NULL,plotit=TRUE)
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
      # calculate ----
      c <- private$O
      if(is.null(c))
      {
        m <- length(s)
        n <- length(x)
        g <- -rho*(x-mu)
        h2 <- rep(sigma^2,n)
        if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
        mskip <- (m-1)*skip+1
        dsskip <- ds/skip
        c <- matrix(0.0,m,n)
        cskip <- matrix(0.0,mskip,n)
        # create and factor A matrix
        A <- private$OUPOptionAtoLU(n,g,h2,r,theta,dsskip,dx)
        # initialize coefficients for terminal value
        message("Option, 1:1")
        cskip[1,] <- V
        c[1,] <- cskip[1,]
        # solve backwards in time
        i <- 1
        im <- 1
        while(i < mskip)
        {
          iskip <- 0
          while(iskip < skip)
          {
            ncnt <- str_length(as.integer(i))
            back <- strrep("\b",ncnt+1)
            iskip <- iskip+1
            i <- i+1
            message(paste(sep="",back,i))
            # create b for this time and solve
            b <- private$OUPOptionb(i,n,cskip,g,h2,r,theta,dsskip,dx)
            cskip[i,] <- private$OUPLUSolve(i,n,A,b)
          }
          ncnt <- str_length(as.integer(i))
          ntime <- str_length(as.integer(im))
          back <- strrep("\b",ncnt+ntime+2)
          im <- im+1
          message(paste(sep="",back,im,":",i))
          c[im,] <- cskip[i,]
        }
        private$O <- c
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotOption(NULL,NULL)) }

      return(list(O=c))
    },
    #' @description
    #' Calculate and plot the envelope of option prices
    #' @param x       vector of n states
    #' @param V       vector of n terminal values
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param r       discount rate 0<=r<inf
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    divide the time interval and report every skip result
    #' @param plotit  TRUE or FALSE
    #' @return list(Ohat(1xn))
    OptionEnvelope = function(x=NULL,V=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,theta=NULL,skip=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(NULL,x,V,r,NULL,theta,skip)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      # calculate ----
      Oenv <- private$Ohat
      if(is.null(Oenv))
      {
        n <- length(x)
        g <- -rho*(x-mu)
        h2 <- rep(sigma^2,n)
        if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
        mskip <- 500*skip+1
        dsskip <- ds/skip
        Oenv <- vector("double",n)
        tsenv <- vector("double",n)
        cskip <- matrix(0.0,mskip,n)
        # create and factor A matrix
        A <- private$OUPOptionAtoLU(n,g,h2,r,theta,dsskip,dx)
        # initialize coefficients for terminal value
        message("Option Envelope, 1")
        cskip[1,] <- V
        Oenv <- cskip[1,]
        tsenv <- rep(0,n)
        # solve backwards in time
        i <- 1
        up <- n
        while(i < mskip & up > 0)
        {
          ncnt <- str_length(as.integer(i))
          back <- strrep("\b",ncnt+1)
          i <- i+1
          message(paste(sep="",back,i))
          # create b for this time and solve
          b <- private$OUPOptionb(i,n,cskip,g,h2,r,theta,dsskip,dx)
          cskip[i,] <- private$OUPLUSolve(i,n,A,b)
          # check for maximum
          up <- 0
          for(j in 1:n)
          {
            if(cskip[i,j] > cskip[i-1,j])
            {
              up <- up+1
              if(cskip[i,j] > Oenv[j])
              {
                Oenv[j] <- cskip[i,j]
                tsenv[j] <- i*dsskip
              }
            }
          }
        }
        private$Ohat <- Oenv
        private$tshat <- tsenv
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotOptionEnvelope(NULL,NULL)) }

      return(list(Ohat=Oenv))
    },
    #' @description
    #' Calculate and plot the decision threshold
    #' @param x       vector of n states
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param V       vector of n terminal values
    #' @param r       discount rate 0<=r<inf
    #' @param phi     search direction for exit or entry options
    #' @param theta   weight of current time in time stepping 0.5<=theta<=1
    #' @param skip    divide the time interval and report every skip result
    #' @param plotit  TRUE or FALSE
    #' @return list(k,Ohat)
    DecisionThreshold = function(x=NULL,V=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,phi=NULL,theta=NULL,skip=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_x_stoch_args(NULL,x,V,r,phi,theta,skip)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- as.double(private$x_stoch_args[[2]])
      V <- private$x_stoch_args[[3]]
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      # calculate ----
      decision <- private$kOhat
      if(is.null(decision))
      {
        n <- length(x)
        g <- -rho*(x-mu)
        h2 <- rep(sigma^2,n)
        if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
        decision <- vector("double",2)
        O <- vector("double",13)
        P <- vector("double",13)
        Oenv <- self$OptionEnvelope(plotit=FALSE)[[1]]
        # search for enter to right
        if(phi > 0 | (phi == 0 & V[n] > V[1]))
        {
          j <- 6
          hit <- FALSE
          decision[1] <- x[n]
          decision[2] <- Oenv[n]
          # moving right, find last point where envelope exceeds terminal value
          while(j < n-2 & hit == FALSE)
          {
            j <- j+1
            if(Oenv[j] > V[j] & Oenv[j+1] <= V[j+1])
            {
              hit <- TRUE
              thisdiff <- Oenv[j]-V[j]
              # extrapolate O[xx] using 7 point Lagrange polynomial
              O[0+7] <- Oenv[j]
              O[-1+7] <- Oenv[j-1]
              O[-2+7] <- Oenv[j-2]
              O[-3+7] <- Oenv[j-3]
              O[-4+7] <- Oenv[j-4]
              O[-5+7] <- Oenv[j-5]
              O[-6+7] <- Oenv[j-6]
              doublehit <- FALSE
              dxx <- (x[j+1]-x[j+0])/100
              xx <- x[j+0]
              k <- 0
              while(k < 200 & doublehit == FALSE)
              {
                k <- k+1
                xx <- xx+dxx
                P[0+7] <- ((xx-x[j-1])*(xx-x[j-2]) * (xx-x[j-3])*(xx-x[j-4])*(xx-x[j-5])*(xx-x[j-6])) /
                  ((x[j-0]-x[j-1])*(x[j-0]-x[j-2])*(x[j-0]-x[j-3])*(x[j-0]-x[j-4])*(x[j-0]-x[j-5])*(x[j-0]-x[j-6]))
                P[-1+7] <- ((xx-x[j-0])*(xx-x[j-2])*(xx-x[j-3])*(xx-x[j-4])*(xx-x[j-5])*(xx-x[j-6])) /
                  ((x[j-1]-x[j-0])*(x[j-1]-x[j-2])*(x[j-1]-x[j-3])*(x[j-1]-x[j-4])*(x[j-1]-x[j-5])*(x[j-1]-x[j-6]))
                P[-2+7] <- ((xx-x[j-0])*(xx-x[j-1])*(xx-x[j-3])*(xx-x[j-4])*(xx-x[j-5])*(xx-x[j-6])) /
                  ((x[j-2]-x[j-0])*(x[j-2]-x[j-1])*(x[j-2]-x[j-3])*(x[j-2]-x[j-4])*(x[j-2]-x[j-5])*(x[j-2]-x[j-6]))
                P[-3+7] <- ((xx-x[j-0])*(xx-x[j-1])*(xx-x[j-2])*(xx-x[j-4])*(xx-x[j-5])*(xx-x[j-6])) /
                  ((x[j-3]-x[j-0])*(x[j-3]-x[j-1])*(x[j-3]-x[j-2])*(x[j-3]-x[j-4])*(x[j-3]-x[j-5])*(x[j-3]-x[j-6]))
                P[-4+7] <- ((xx-x[j-0])*(xx-x[j-1])*(xx-x[j-2])*(xx-x[j-3])*(xx-x[j-5])*(xx-x[j-6])) /
                  ((x[j-4]-x[j-0])*(x[j-4]-x[j-1])*(x[j-4]-x[j-2])*(x[j-4]-x[j-3])*(x[j-4]-x[j-5])*(x[j-4]-x[j-6]))
                P[-5+7] <- ((xx-x[j-0])*(xx-x[j-1])*(xx-x[j-2])*(xx-x[j-3])*(xx-x[j-4])*(xx-x[j-6])) /
                  ((x[j-5]-x[j-0])*(x[j-5]-x[j-1])*(x[j-5]-x[j-2])*(x[j-5]-x[j-3])*(x[j-5]-x[j-4])*(x[j-5]-x[j-6]))
                P[-6+7] <- ((xx-x[j-0])*(xx-x[j-1])*(xx-x[j-2])*(xx-x[j-3])*(xx-x[j-4])*(xx-x[j-5])) /
                  ((x[j-6]-x[j-0])*(x[j-6]-x[j-1])*(x[j-6]-x[j-2])*(x[j-6]-x[j-3])*(x[j-6]-x[j-4])*(x[j-6]-x[j-5]))
                O[1+7] <- O[0+7]*P[0+7]+O[-1+7]*P[-1+7]+O[-2+7]*P[-2+7]+O[-3+7]*P[-3+7]+O[-4+7]*P[-4+7]+O[-5+7]*P[-5+7]+O[-6+7]*P[-6+7]
                #compare with interpolated V
                VV <- V[j-0]+(V[j+1]-V[j-0]) / (x[j+1]-x[j-0])*(xx-x[j-0])
                prevdiff <- thisdiff
                thisdiff <- O[1+7]-VV
                if(thisdiff < 0 | thisdiff > prevdiff)
                {
                  doublehit <- TRUE
                  decision[1] <- xx
                  decision[2] <- O[1+7]
                }
              }
            }
          }
        }
        # search for exit to left
        else
        {
          j <- n-5
          hit <- FALSE
          decision[1] <- x[1]
          decision[2] <- Oenv[1]
          # moving left, find last point where envelope exceeds terminal value
          while(j > 3 & hit == FALSE)
          {
            j <- j-1
            if(Oenv[j] > V[j] & Oenv[j-1] <= V[j-1])
            {
              hit <- TRUE
              thisdiff <- Oenv[j]-V[j]
              # extrapolate O(xx) using 7 point Lagrange polynomial
              O[0+7] <- Oenv[j]
              O[1+7] <- Oenv[j+1]
              O[2+7] <- Oenv[j+2]
              O[3+7] <- Oenv[j+3]
              O[4+7] <- Oenv[j+4]
              O[5+7] <- Oenv[j+5]
              O[6+7] <- Oenv[j+6]
              doublehit <- FALSE
              dxx <- (x[j+0]-x[j-1])/100
              xx <- x[j-0]
              k <- 200
              while(k > 0 & doublehit == FALSE)
              {
                k <- k-1
                xx <- xx-dxx
                P[0+7] <- ((xx-x[j+1])*(xx-x[j+2])*(xx-x[j+3])*(xx-x[j+4])*(xx-x[j+5])*(xx-x[j+6]))/
                  ((x[j+0]-x[j+1])*(x[j+0]-x[j+2])*(x[j+0]-x[j+3])*(x[j+0]-x[j+4])*(x[j+0]-x[j+5])*(x[j+0]-x[j+6]))
                P[1+7] <- ((xx-x[j+0])*(xx-x[j+2])*(xx-x[j+3])*(xx-x[j+4])*(xx-x[j+5])*(xx-x[j+6]))/
                  ((x[j+1]-x[j+0])*(x[j+1]-x[j+2])*(x[j+1]-x[j+3])*(x[j+1]-x[j+4])*(x[j+1]-x[j+5])*(x[j+1]-x[j+6]))
                P[2+7] <- ((xx-x[j+0])*(xx-x[j+1])*(xx-x[j+3])*(xx-x[j+4])*(xx-x[j+5])*(xx-x[j+6]))/
                  ((x[j+2]-x[j+0])*(x[j+2]-x[j+1])*(x[j+2]-x[j+3])*(x[j+2]-x[j+4])*(x[j+2]-x[j+5])*(x[j+2]-x[j+6]))
                P[3+7] <- ((xx-x[j+0])*(xx-x[j+1])*(xx-x[j+2])*(xx-x[j+4])*(xx-x[j+5])*(xx-x[j+6]))/
                  ((x[j+3]-x[j+0])*(x[j+3]-x[j+1])*(x[j+3]-x[j+2])*(x[j+3]-x[j+4])*(x[j+3]-x[j+5])*(x[j+3]-x[j+6]))
                P[4+7] <- ((xx-x[j+0])*(xx-x[j+1])*(xx-x[j+2])*(xx-x[j+3])*(xx-x[j+5])*(xx-x[j+6]))/
                  ((x[j+4]-x[j+0])*(x[j+4]-x[j+1])*(x[j+4]-x[j+2])*(x[j+4]-x[j+3])*(x[j+4]-x[j+5])*(x[j+4]-x[j+6]))
                P[5+7] <- ((xx-x[j+0])*(xx-x[j+1])*(xx-x[j+2])*(xx-x[j+3])*(xx-x[j+4])*(xx-x[j+6]))/
                  ((x[j+5]-x[j+0])*(x[j+5]-x[j+1])*(x[j+5]-x[j+2])*(x[j+5]-x[j+3])*(x[j+5]-x[j+4])*(x[j+5]-x[j+6]))
                P[6+7] <- ((xx-x[j+0])*(xx-x[j+1])*(xx-x[j+2])*(xx-x[j+3])*(xx-x[j+4])*(xx-x[j+5]))/
                  ((x[j+6]-x[j+0])*(x[j+6]-x[j+1])*(x[j+6]-x[j+2])*(x[j+6]-x[j+3])*(x[j+6]-x[j+4])*(x[j+6]-x[j+5]))
                O[-1+7] <- O[0+7]*P[0+7]+O[1+7]*P[1+7]+O[2+7]*P[2+7]+O[3+7]*P[3+7]+O[4+7]*P[4+7]+O[5+7]*P[5+7]+O[6+7]*P[6+7]
                #compare with interpolated V
                VV <- V[j+0]+(V[j-1]-V[j+0])/(x[j+0]-x[j-1])*(x[j+0]-xx)
                prevdiff <- thisdiff
                thisdiff <- O[-1+7]-VV
                if(thisdiff < 0 | thisdiff > prevdiff)
                {
                  doublehit <- TRUE
                  decision[1] <- xx
                  decision[2] <- O[-1+7]
                }
              }
            }
          }
        }
        private$kOhat <- decision
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDecisionThreshold(NULL)) }
      return(list(k=decision[1],Ohat=decision[2]))
    },
    # public plot methods ----
    #' @description
    #' Plot drifts
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDrift = function(title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      x <- private$x_stoch_args[[2]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      drift <- private$g
      if(is.null(drift)) { drift <- self$Drift(plotit=FALSE)[[1]] }
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        n <- length(x)
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
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot diffusions
    #' @param type  = 2 and 3 for 2D
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDiffusion = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # get ----
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      drift <- private$g
      if(is.null(drift)) { drift <- self$Drift(plotit=FALSE)[[1]] }
      diffusion <- private$h2
      if(is.null(diffusion)) { diffusion <- self$Diffusion(plotit=FALSE)[[1]]  }
      n <- length(x)
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        diffusion <- diffusion[Ixbeg:Ixend]
        n <- length(x)
      }
      sqrt <- diffusion^0.5
      driftplus <- drift+sqrt
      driftminus <- drift-sqrt
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
      if(type < 2.5)
      {
        if(is.null(yaxis)) { yaxis <- "<i>g</i>(<i>x</i>)&plusmn;<i>h</i>" }
        lookleft <- list(text=yaxis)
        horz=list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert=list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
        diffusionline <- list(color=mgn$d,width=4)
        driftline <- list(color=cyn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Diffusion2Dg")
        legendpos <- list(orientation="h",x=1.05,y=1.0,xanchor="right")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=drift,name="<i>g</i>(<i>x</i>)",mode="lines",line=driftline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=driftplus,name="<i>g</i>(<i>x</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=driftminus,name="<i>g</i>(<i>x</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",showlegend=FALSE,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions) %>%
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
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot terminal values
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotTerminalValue = function(title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # get ----
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
      V_info <- private$V_info
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      n <- length(V)
      Ix <- private$V_info$Ix
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        n <- length(x)
      }
      # plot ----
      # OUP_FD_Terminal2D
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- ""
        if(is.null(title)) { title <- private$V_info$name }
        if(Ix == 1)
        {
          xo <- private$V_linear_args[[1]]
          vs <- private$V_linear_args[[2]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>v</i><sub>1</sub>=",format(vs,digits=4),")",esml)
        }
        else if(Ix == 2)
        {
          xo <- private$V_degenerate_args[[1]]
          Vmax <- private$V_degenerate_args[[2]]
          Vmin <- private$V_degenerate_args[[3]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 3)
        {
          xo <- private$V_stepped_args[[1]]
          vs <- private$V_stepped_args[[2]]
          Vmax <- private$V_stepped_args[[3]]
          Vmin <- private$V_stepped_args[[4]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>v</i><sub>1</sub>=",format(vs,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 4)
        {
          xo <- private$V_kinked_args[[1]]
          vs <- private$V_kinked_args[[2]]
          Vmax <- private$V_kinked_args[[3]]
          Vmin <- private$V_kinked_args[[4]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>v</i><sub>1</sub>=",format(vs,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 5)
        {
          xo <- private$V_butterfly_args[[1]]
          xm <- private$V_butterfly_args[[2]]
          vs <- private$V_butterfly_args[[3]]
          Vmax <- private$V_butterfly_args[[4]]
          Vmin <- private$V_butterfly_args[[5]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>x</i><sub>m</sub>=",format(xm,digits=4),",<i>v</i><sub>s</sub>=",format(vs,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 6)
        {
          xo <- private$V_mitscherlich_args[[1]]
          vr <- private$V_mitscherlich_args[[2]]
          Vmax <- private$V_mitscherlich_args[[3]]
          Vmin <- private$V_mitscherlich_args[[4]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>v</i><sub>2</sub>=",format(vr,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 7)
        {
          xi <- private$V_gompertz_args[[1]]
          vr <- private$V_gompertz_args[[2]]
          Vmax <- private$V_gompertz_args[[3]]
          Vmin <- private$V_gompertz_args[[4]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>1</sub>=",format(xi,digits=4),",<i>v</i><sub>2</sub>=",format(vr,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 8)
        {
          xi <- private$V_logistic_args[[1]]
          vr <- private$V_logistic_args[[2]]
          Vmax <- private$V_logistic_args[[3]]
          Vmin <- private$V_logistic_args[[4]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>1</sub>=",format(xi,digits=4),",<i>v</i><sub>2</sub>=",format(vr,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 9)
        {
          xo <- private$V_transcendental_args[[1]]
          xi <- private$V_transcendental_args[[2]]
          xm <- private$V_transcendental_args[[3]]
          Vmax <- private$V_transcendental_args[[4]]
          Vmin <- private$V_transcendental_args[[5]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>x</i><sub>1</sub>=",format(xi,digits=4),",<i>x</i><sub>2</sub>=",format(xm,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
        else if(Ix == 10)
        {
          xo <- private$V_yieldindex_args[[1]]
          xi <- private$V_yieldindex_args[[2]]
          xm <- private$V_yieldindex_args[[3]]
          Vmax <- private$V_yieldindex_args[[4]]
          Vmin <- private$V_yieldindex_args[[5]]
          syms <- paste(sep="",bsml,"(<i>x</i><sub>0</sub>=",format(xo,digits=4),",<i>x</i><sub>1</sub>=",format(xi,digits=4),",<i>x</i><sub>2</sub>=",format(xm,digits=4),",<i>V</i><sub>max</sub>=",format(Vmax,digits=4),",<i>V</i><sub>min</sub>=",format(Vmin,digits=4),")",esml)
        }
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
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot options
    #' @param type  = 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(s)
      n <- length(x)
      options <- private$O
      if(is.null(options)) { options <- self$Option(plotit=FALSE)[[1]] }
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m | Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        options <- options[Ixend:Ixbeg,]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend]
        n <- length(x)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family:Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),bsym,",<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>q</i>=",esym,format(theta,digits=4),",<i>ds</i>",bsym,"=",esym,format(ds,digits=4),",<i>skip</i>",bsym,"=",esym,skip,",<i>dx</i>",bsym,"=",esym,format(dx,digits=2),")",esml)
        if(is.null(title)) { title <- "Option"}
      }
      else if(is.null(title)) { title <- ""}
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_FD_Option2D
      if(type < 2.5)
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
        ds <- as.integer((m-1)/10)
        if(ds < 1) { ds <- 1 }
        optionline <- list(color=red$d,width=4)
        lineopacity <- 1
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Option2D")
        fig <- plot_ly()
        i <- 1
        while(i <= m)
        {
          fig <- add_trace(fig,type="scatter",x=x,y=options[i,],name=paste(sep="","&#x1D546;(",s[i],"<i>,x</i>)"),mode="lines",line=optionline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+ds
          lineopacity <- lineopacity-0.05
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 3D
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
        hover <- list(bgcolor=red$e,font=list(color=red$b))
        if(V[n] > V[1]) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(V[n] == V[1]) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_FD_Option3DSurface and OUP_FD_Option3DSurfaceScatter
        if(type < 4.5)
        {
          gradient <- list(c(0,red$c),c(1,red$c))
          rise <- list(x=0,y=-300,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=x,y=s,z=options,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
          # OUP_FD_Option3DSurface
          if(type < 3.5) { imageoptions=list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_Option3DSurface") }
          # OUP_FD_Option3DSurfaceScatter
          else
          {
            optionline <- list(color=red$e,width=8)
            imageoptions=list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_Option3DSurfaceScatter")
            ss <- vector("double",n)
            ds <- as.integer((m-1)/10)
            if(ds < 1) { ds <- 1 }
            i <- 1
            while(i <= m)
            {
              for(j in 1:n) { ss[j] <- s[i] }
              fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],mode="lines",line=optionline,hoverinfo="text",text=coordinates[i,])
              i <- i+ds
            }
          }
        }
        # OUP_FD_Option3DScatter
        else
        {
          optionline <- list(color=red$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_Option3DScatter")
          fig <- plot_ly()
          ss <- vector("double",n)
          ds <- as.integer((m-1)/10)
          if(ds < 1) { ds <- 1 }
          i <- 1
          while(i <= m)
          {
            for(j in 1:n) { ss[j] <- s[i] }
            fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=options[i,],mode="lines",line=optionline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,])
            i <- i+ds
            lineopacity <- lineopacity-0.05
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the option envelope
    #' @param type  = 3 for 2D, 4 for 3D
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
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      s <- private$x_stoch_args[[1]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
      r <- private$x_stoch_args[[4]]
      phi <- private$x_stoch_args[[5]]
      theta <- private$x_stoch_args[[6]]
      skip <- private$x_stoch_args[[7]]
      ds <- private$x_stoch_args[[8]]
      dx <- private$x_stoch_args[[9]]
      type <- private$plot_info$plottype$type
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
      m <- length(s)
      n <- length(x)
      Oenv <- private$Ohat
      if(is.null(Oenv)) { Oenv <- self$OptionEnvelope(plotit=FALSE)[[1]] }
      tsenv <- s[1]-private$tshat
      options <- private$O
      if(is.null(options)) { options <- self$Option(plotit=FALSE)[[1]] }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        Oenv <- Oenv[Ixbeg:Ixend]
        tsenv <- tsenv[Ixbeg:Ixend]
        options <- options[,Ixbeg:Ixend]
        n <- length(x)
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
      if(type < 3.5)
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
        Oenvline <- list(color=red$d,width=4)
        terminalline <- list(color=gry$c,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Envelope2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=V,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=Oenv,name="\u00D4(<i>x</i>)",mode="lines",line=Oenvline,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions) %>%
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
        Ohold <- vector("double",n)
        Oexercise <- vector("double",n)
        coordinatesenv <- vector("double",n)
        coordinates <- matrix("",m,n)
        for(j in 1:n)
        {
          if(tsenv[j] == s[1])
          {
            Ohold[j] <- NA
            Oexercise[j] <- Oenv[j]
          }
          else
          {
            Ohold[j] <- Oenv[j]
            Oexercise[j] <- NA
          }
          coordinatesenv[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(Oenv[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          for(i in 1:m) { coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
        }
        hover <- list(bgcolor=red$e,font=list(color=red$b))
        if(V[n] > V[1]) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(V[n] == V[1]) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(x[1],x[n]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(s[m],s[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        Oholdline <- list(color=red$e,width=10)
        Oexerciseline <- list(color=red$e,width=8)
        Oenvline <- list(dash="dash",color=red$e,width=6)
        gradient <- list(c(0,red$c),c(1,red$c))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_FD_OptionEnvelope3D")
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=x,y=s,z=options,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
          add_trace(.,type="scatter3d",x=x,y=tsenv,z=Oenv,mode="lines",line=Oenvline,hoverinfo="text",text=coordinatesenv) %>%
          add_trace(.,type="scatter3d",x=x,y=tsenv,z=Ohold,mode="lines",line=Oholdline,hoverinfo="text",text=coordinatesenv) %>%
          add_trace(.,type="scatter3d",x=x,y=tsenv,z=Oexercise,mode="lines",line=Oexerciseline,hoverinfo="text",text=coordinatesenv) %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,hoverlabel=hover,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the decision threshold
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param xbeg    begin value for state axis
    #' @param xend    end value for state axis
    #' @return plot
    PlotDecisionThreshold = function(title=NULL,xaxis=NULL,yaxis=NULL,xbeg=NULL,xend=NULL)
    {
      # get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      x <- private$x_stoch_args[[2]]
      V <- private$x_stoch_args[[3]]
      if(is.null(V)) { V <- self$TerminalValue(plotit=FALSE)[[1]] }
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
      n <- length(x)
      Oenv <- private$Ohat
      if(is.null(Oenv)) { Oenv <- self$OptionEnvelope(plotit=FALSE)[[1]] }
      decision <- private$kOhat
      if(is.null(decision)) { decision <- self$DecisionThreshold(plotit=FALSE) }
      k <- decision[[1]]
      Ohat <- decision[[2]]
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        V <- V[Ixbeg:Ixend]
        Oenv <- Oenv[Ixbeg:Ixend]
        n <- length(x)
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
      if(phi > 0 | (phi == 0 & V[n] > V[1])) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero",side="right") }
      else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
      Oenvline <- list(color=red$d,width=4)
      terminalline <- list(color=gry$c,width=4)
      oline <- list(dash="dot",color=red$d,width=4)
      kline <- list(dash="dot",color=red$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_FD_Decision2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=V,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=x,y=Oenv,name="\u00D4(<i>x</i>)",mode="lines",line=Oenvline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(k,k),y=c(0,Ohat),name="<i>k</i>",mode="lines",line=kline,hoverinfo="x+y")
      if(phi > 0 | (phi == 0 & V[n] > V[1]))
      {
        fig <- add_trace(fig,type="scatter",x=c(x[n],k),y=c(Ohat,Ohat),name="\u00D4",mode="lines",line=oline,hoverinfo="x+y")
        kOhat <- list(x=k,y=Ohat,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(Ohat,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",showarrow=FALSE)
      }
      else
      {
        fig <- add_trace(fig,type="scatter",x=c(x[1],k),y=c(Ohat,Ohat),name="\u00D4",mode="lines",line=oline,hoverinfo="x+y")
        kOhat <- list(x=k,y=Ohat,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(Ohat,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",showarrow=FALSE)
      }
      fig <- config(fig,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,annotations=kOhat,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

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
    plot_info = NULL,
    default = NULL,
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
    tshat = NULL,
    # private output fields ----
    O = NULL,
    Ohat = NULL,
    kOhat = NULL,
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
        if(is.numeric(input[1]) & is.finite(input[1]) & !is.na(input[1])) { sca <- input[1] }
        else { message(paste(sep="",input," is not a real number:")) }
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
        while(i < n & OK == TRUE)
        {
          i <- i+1
          if(!is.numeric(input[i]) | is.infinite(input[i]) | is.na(input[i])) { OK <- FALSE }
        }
        if(OK == TRUE)
        {
          if(updown > 0) { vec <- sort(input) }
          else if(updown == 0) { vec <- input }
          if(updown == 0) { vec <- input }
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
        if(is.character(input)) { chr <- input[1] }
        else { message(paste("non-character input:",input)) }
      }
      return(chr)
    },
    extract_boolean = function(input)
    {
      bool <- NULL
      if(!is.null(input))
      {
        if(is.list(input)) { input <- input[[1]] }
        if(input[1] == TRUE | input[1] == FALSE) { bool <- input[1] }
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
              if(!is.null(red) & !is.null(grn) & !is.null(blu)) { clr <- paste(sep="","rgb(",red,",",grn,",",blu,")") }
              else { message("unrecognized hexadecimal color:",rgb)}
            }
            else if(n > 9)
            {
              if(tolower(str_sub(rgb,1,4)) == "rgb(" & str_sub(rgb,n,n) == ")")
              {
                OK <- TRUE
                hit <- FALSE
                k <- 4
                j <- k
                while(OK == TRUE & hit == FALSE & k < n)
                {
                  k <- k+1
                  token <- str_sub(rgb,k,k)
                  if(token == "," | token == ")")
                  {
                    hit <- TRUE
                    if(k-j == 1) { OK <- FALSE }
                  }
                  else if(token != "0" & token != "1" & token != "2" & token != "3" & token != "4" & token != "5" & token != "6" & token != "7" & token != "8" & token != "9") { OK <- FALSE }
                }
                if(OK == TRUE)
                {
                  red <- as.integer(str_sub(rgb,j+1,k-1))
                  if(red >= 0 & red <= 255)
                  {
                    hit <- FALSE
                    j <- k
                    while(OK == TRUE & hit == FALSE & k < n)
                    {
                      k <- k+1
                      token <- str_sub(rgb,k,k)
                      if(token == "," | token == ")")
                      {
                        hit <- TRUE
                        if(k-j == 1) { OK <- FALSE }
                      }
                      else if(token != "0" & token != "1" & token != "2" & token != "3" & token != "4" & token != "5" & token != "6" & token != "7" & token != "8" & token != "9") { OK <- FALSE }
                    }
                    if(OK == TRUE)
                    {
                      grn <- as.integer(str_sub(rgb,j+1,k-1))
                      if(grn >= 0 & grn <= 255)
                      {
                        hit <- FALSE
                        j <- k
                        while(OK == TRUE & hit == FALSE & k < n)
                        {
                          k <- k+1
                          token <- str_sub(rgb,k,k)
                          if(token == "," | token == ")")
                          {
                            hit <- TRUE
                            if(k-j == 1) { OK <- FALSE }
                          }
                          else if(token != "0" & token != "1" & token != "2" & token != "3" & token != "4" & token != "5" & token != "6" & token != "7" & token != "8" & token != "9") { OK <- FALSE }
                        }
                        if(OK == TRUE)
                        {
                          blu <- as.integer(str_sub(rgb,j+1,k-1))
                          if(blu >= 0 & blu <= 255) { clr <- paste(sep="","rgb(",red,",",grn,",",blu,")") }
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
    vecs_equal = function(vec1,vec2)
    {
      tolerance <- 1e-12
      allequal <- TRUE
      n1 <- length(vec1)
      n2 <- length(vec2)
      if(n1 == n2)
      {
        j <- 0
        while(j < n1 & allequal == TRUE )
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
    OUPOptionAtoLU = function(n,g,h2,r,theta,ds,dx)
    {
      # create A matrix
      # first row has 3 entries
      A <- matrix(0.0,n,n)
      A[1,1] <- 1/ds+theta*r+0.5*theta*(3*g[1]/dx-h2[1]/dx^2)
      A[1,2] <- -theta*(2*g[1]/dx-h2[1]/dx^2)
      A[1,3] <- 0.5*theta*(g[1]/dx-h2[1]/dx^2)
      # middle rows are tridiagonal
      j <-1
      while(j < n-1)
      {
        j <- j+1
        A[j,j-1] <- 0.5*theta*(g[j]/dx-h2[j]/dx^2)
        A[j,j] <- 1/ds+theta*r+theta*h2[j]/dx^2
        A[j,j+1] <- -0.5*theta*(g[j]/dx+h2[j]/dx^2)
      }
      # last row has 3 entries
      A[n,n-2] <- -0.5*theta*(g[n]/dx+h2[n]/dx^2)
      A[n,n-1] <- theta*(2*g[n]/dx+h2[n]/dx^2)
      A[n,n] <- 1/ds+theta*r-0.5*theta*(3*g[n]/dx+h2[n]/dx^2)

      # LU factorisation
      # begin in row 2
      j <- 0
      while(j < n-3)
      {
        j <- j+1
        i <- j+1
        A[i,j] <- A[i,j]/A[j,j]
        k <- j
        while(k < j+2)
        {
          k <- k+1
          A[i,k] <- A[i,k]-A[i,j]*A[j,k]
        }
      }
      # last rows done this way because last row has 3 elements
      j <- n-3
      while(j < n-1)
      {
        j <- j+1
        i <- j
        while(i < n)
        {
          i <- i+1
          A[i,j] <- A[i,j]/A[j,j]
          k <- j
          while(k < n)
          {
            k <- k+1
            A[i,k] <- A[i,k]-A[i,j]*A[j,k]
          }
        }
      }
      # LU is in A with implicit 1 on diagonal of L
      return(A)
    },
    OUPOptionb = function(i,n,c,g,h2,r,theta,ds,dx)
    {
      b <- vector("double",n)
      b[1] <- (1/ds-(1-theta)*r)*c[i-1,1]+
              0.5*(1-theta)*(g[1]*(-3*c[i-1,1]+4*c[i-1,2]-c[i-1,3])/dx+
              h2[1]*(c[i-1,1]-2*c[i-1,2]+c[i-1,3])/dx^2)
      j <- 1
      while(j < n-1)
      {
        j <- j+1
        b[j] <- (1/ds-(1-theta)*r)*c[i-1,j]+
                0.5*(1-theta)*(g[j]*(-c[i-1,j-1]+c[i-1,j+1])/dx+
                h2[j]*(c[i-1,j-1]-2*c[i-1,j]+c[i-1,j+1])/dx^2)
      }
      b[n] <- (1/ds-(1-theta)*r)*c[i-1,n]+
              0.5*(1-theta)*(g[n]*(c[i-1,n-2]-4*c[i-1,n-1]+3*c[i-1,n])/dx+
              h2[n]*(c[i-1,n-2]-2*c[i-1,n-1]+c[i-1,n])/dx^2)
      return(b)
    },
    OUPLUSolve = function(i,n,A,b)
    {
      # forward substitution to find q as solution of Lq=b
      j <- 1
      while(j < n-1)
      {
        j <- j+1
        b[j] <- b[j]-A[j,j-1]*b[j-1]
      }
      b[n] <- b[n]-A[n,n-2]*b[n-2]-A[n,n-1] * b[n-1]
      # backward substitution to find c as solution of Uc=q
      c <- vector("double",n)
      c[n] <- b[n]/A[n,n]
      j <- n
      while(j > 2)
      {
        j <- j-1
        c[j] <- (b[j]-A[j,j+1]*c[j+1])/A[j,j]
      }
      c[1] <- (b[1]-A[1,2]*c[2]-A[1,3]*c[3])/A[1,1]
      return(c)
    },
    # private plot methods ----
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
              while(j < n & hit == FALSE)
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
              while(j > 1 & hit == FALSE)
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
              while(j > 1 & hit == FALSE)
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
              while(j < n & hit == FALSE)
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
    }
  )
)
