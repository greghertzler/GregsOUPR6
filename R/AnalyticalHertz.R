library(R6)
library(plotly)
library(stringr)

# roxygen ----
#' R6 Class implementing Analytical formulas for an Ornstein-Uhlenbeck Process
#'
#' @description
#' A set of functions to calculate probabilities, option prices, visiting
#'  times, first passage times and decision thresholds--everything for
#'  simple exit and entry options in a Real Options Analysis.
#'
#' @details # Formulas:
#'     z stochastic
#'       Drift
#'       Diffusion
#'     y stochastic
#'       Mean
#'       MeanToConverge
#'       Variance
#'       VarianceToConverge
#'       Density
#'       Probability
#'       DoubleIntegral
#'     x stochastic
#'       Option
#'       OptionEnvelope
#'       DecisionThreshold
#'       Obligation
#'     t stochastic
#'       PassageTimeMode
#'       PassageTimeMedian
#'       PassateTimeMean
#'       PassageTimeVariance
#'       PassageTimePercentiles
#'       PassageTimeDensity
#'       PassageTimeProbability
#'
#' @details # Plots
#'       PlotDrift
#'       PlotDiffusion
#'       PlotMean
#'       PlotMeanToConverge
#'       PlotVariance
#'       PlotVarianceToConverge
#'       PlotDensity
#'       PlotProbability
#'       PlotDoubleIntegral
#'       PlotOption
#'       PlotOptionEnvelope
#'       PlotDecisionThreshold
#'       PlotObligation
#'       PlotPassageTimeModeMedianMean
#'       PlotPassageTimeVariance
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
#'
#' @details # Using the formulas:
#' Demonstration scripts are files in the 'demo' directory. Identify a formula
#'  and in the console type:
#'
#'       demo(A_FormulaName), or
#'       demo(A_PlotFormulaName).
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
#'  Difference method, an analytical formula can calculate an option price
#'  up to 10,000 times quicker. Compared with Monte Carlo simulation, an
#'  analytical formula can calculate passage time probabilities up to
#'  100,000 times quicker.
#'
#' The bottom line of a real options analysis is the decision threshold
#'  and the mode, median and mean passage times. The many functions available
#'  are used in the calculations, but to go straight to the bottom line:
#'
#'       A <- Analytical$new()
#'       A$DecisionThreshold()
#'       A$PassageTimeMean()
#'       A$PassageTimeVariance()
#'       A$PlotPassageTimePercentiles(type=1)
#'       A$PlotPassageTimePercentiles(type=2)
#'
#' Perhaps a better measure of passage times are the percentiles. The mean
#'  does not exist if rho is zero and the variance does not exist if rho is
#'  small. The mode and median always exist.  Percentiles, which include the
#'  median, always exist.
#'
#' The passage time calculations are challenging. For example, the passage
#'  time density is too complicated to explain in this short discussion.
#'  But don't be surprised if it is bi-modal or even negative. Care must be
#'  taken in interpreting passage times.
#'
#' Examples in R6 don't work the same way as other R modules.  There is only
#'  one example for an R6 object, not one for each function in the object.
#'  To run examples, the devtools::run_examples() works, but the R command
#'  example("Analytical") doesn't.  You can copy commands to the clipboard,
#'  paste into the console and press Enter. Examples in this help and a simple
#'  example at the bottom can be run in this way.
#'
#' A better alternative is demo(), but this works sometimes, sometimes not.
#'  If you don't see the demos for this package, go to Files and demo.  You will
#'  see the demo names and then type something like demo(A_DecisionThreshold).

# class ----
#' @import plotly
#' @import stringr
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
    #' @examples
    #'   A <- Analytical$new()
    #'   A$DecisionThreshold()
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
      private$t_stoch_args <- list(t=seq(from=0,to=10,by=0.1),k=0,s=0,x=15,z=seq(from=-30,to=30,by=6),omega=1,Ppct=0.75)
      private$psiphi <- -1
      private$undo_args <- list(list(oup_params=private$oup_params,z_stoch_args=private$z_stoch_args,y_stoch_args=private$y_stoch_args,x_stoch_args=private$x_stoch_args,t_stoch_args=private$t_stoch_args))
      private$undoIx <- 1
      # plot info ----
      plottype <- list(type=3,pmax=0.066,ptmax=1.1)
      plotfont <- list(family="Cambria",size=14)
      plotfile <- list(format="png",width=640,height=480)
      plottheme <- list(name="dark",opaque=1.0)
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
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"A") }
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
    #' @return list(t,y,s,x,psi,eps)
    set_y_stoch_args = function(t=NULL,y=NULL,s=NULL,x=NULL,psi=NULL,eps=NULL)
    {
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
    #' @param who identifier for object sending the parameters
    #' @return list(s,x,t,y,r,phi)
    set_x_stoch_args = function(s=NULL,x=NULL,t=NULL,y=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,r,phi,"A") }
      hit <- FALSE
      if(!is.null(s))
      {
        vec <- private$extract_vector(s,-1)
        if(!is.null(vec))
        {
          if(!private$vecs_equal(vec,private$x_stoch_args$s))
          {
            private$x_stoch_args$s <- vec
            private$Oneg <- NULL
            private$Opos <- NULL
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
          if(n > 100 | is.null(private$OUP))
          {
            if(!private$vecs_equal(vec,private$x_stoch_args$x))
            {
              private$x_stoch_args$x <- vec
              private$Oneg <- NULL
              private$Opos <- NULL
              private$Oshatneg <- NULL
              private$Oshatpos <- NULL
              private$Oscarfneg <- NULL
              private$Oscarfpos <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
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
            private$Oneg <- NULL
            private$Opos <- NULL
            private$Oshatneg <- NULL
            private$Oshatpos <- NULL
            private$Oscarfneg <- NULL
            private$Oscarfpos <- NULL
            private$KOneg <- NULL
            private$KOpos <- NULL
            private$BCneg <- NULL
            private$BCpos <- NULL
            hit <- TRUE
          }
        }
        else { message("c not set.") }
      }
      sn <- private$x_stoch_args$s[1]
      if(private$x_stoch_args$t < sn)
      {
        private$x_stoch_args$t <- sn
        message(paste(sep="","t has been set to ",sn,"."))
        private$Oneg <- NULL
        private$Opos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
      }
      if(hit == TRUE) { private$psiphi <- private$x_stoch_args$phi}
      return(private$x_stoch_args)
    },
    #' @description
    #' Set t stochastic arguments
    #' @param t     vector of m times s<=t<inf
    #' @param k     decision threshold -inf<k<inf
    #' @param s     initial time -inf<s<inf
    #' @param x     initial state -inf<x<inf
    #' @param z     vector of n alternate initial states -inf<z<inf
    #' @param omega degree of irreversibility 0<=omega<=1
    #' @param Ppct  passage time probability for a percentile  0.01<=Ppct<=0.99
    #' @return list(t,k,s,x,z,omega)
    set_t_stoch_args = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,Ppct=NULL)
    {
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
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
            private$tmode <- NULL
            private$tmedian <- NULL
            private$tmean <- NULL
            private$tvariance <- NULL
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
            private$tmodes <- NULL
            private$tmedians <- NULL
            private$tmeans <- NULL
            private$tvariances <- NULL
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
            private$tmode <- NULL
            private$tmodes <- NULL
            private$tmedian <- NULL
            private$tmedians <- NULL
            private$tmean <- NULL
            private$tmeans <- NULL
            private$tvariance <- NULL
            private$tvariances <- NULL
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
        private$tmode <- NULL
        private$tmodes <- NULL
        private$tmedian <- NULL
        private$tmedians <- NULL
        private$tmean <- NULL
        private$tmeans <- NULL
        private$tvariance <- NULL
        private$tvariances <- NULL
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
    #' Set information for plotting
    #' @param type       = 1 and 2 for 2D, 3 and 4 for 3D
    #' @param pmax       maximum vertical axis for transition density plot
    #' @param ptmax      maximum vertical axis for passage time density plot
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
    set_plot_info = function(type=NULL,pmax=NULL,ptmax=NULL,fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_plot_info(type,pmax,ptmax,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,"A")}
      if(!is.null(type))
      {
        sca <- private$extract_scalar(type)
        if(!is.null(sca)) { private$plot_info$plottype$type <- sca }
        else { message("type not set.") }
      }
      if(!is.null(pmax))
      {
        if(is.numeric(pmax))
        {
          if(is.infinite(pmax) | is.nan(pmax)) { private$plot_info$plottype$pmax <- NaN }
          else
          {
            sca <- private$extract_scalar(pmax)
            if(sca <= 0) { sca <- NaN }
            private$plot_info$plottype$pmax <- sca
          }
        }
        else { message("pmax not set.") }
      }
      if(!is.null(ptmax))
      {
        if(is.numeric(ptmax))
        {
          if(is.infinite(ptmax) | is.nan(ptmax)) { private$plot_info$plottype$ptmax <- NaN }
          else
          {
            sca <- private$extract_scalar(ptmax)
            if(sca <= 0) { sca <- NaN }
            private$plot_info$plottype$ptmax <- sca
          }
        }
        else { message("ptmax not set.") }
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
    #' Get all arguments, time series and information
    #' @return list(oup_params,z_stoch_args,y_stoch_args,x_stoch_args,t_stoch_args,plot_info)
    get_all = function()
    {
      all <- list(oup_params = private$oup_params,
        z_stoch_args = private$z_stoch_args,
        y_stoch_args = private$y_stoch_args,
        x_stoch_args = private$x_stoch_args,
        t_stoch_args = private$t_stoch_args,
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
    #' Get information for plotting
    #' @return list(type,font,file,theme,3D,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
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
      if(!private$vecs_equal(zseq,private$z_stoch_args$z))
      {
        private$z_stoch_args$z <- zseq
        private$g <- NULL
        private$h2 <- NULL
      }
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
      eps <- private$y_stoch_args[[6]]
      # state
      y <- 4*abs(mu)
      if(y < 1) { y <- 1}
      yscale <- 1
      while(y > yscale) { yscale <- 10*yscale }
      y <- round(y/yscale,1)*yscale
      yby <- y/100
      yscale <- 1
      while(abs(mu) > yscale) { yscale <- 10*yscale }
      x <- round(mu/yscale,2)*yscale
      if(mu < 0) { yfrom <- x-30*yby }
      else { yfrom <- x-70*yby }
      yto <- yfrom+100*yby
      yseq <- seq(from=yfrom,to=yto,by=yby)
      if(!private$vecs_equal(yseq,private$y_stoch_args$y))
      {
        private$y_stoch_args$y <- yseq
        private$p <- NULL
        private$Pneg <- NULL
        private$Ppos <- NULL
        private$PPneg <- NULL
        private$PPpos <- NULL
      }
      if(-x != private$y_stoch_args$x)
      {
        private$y_stoch_args$x <- -x
        private$G <- NULL
        private$Gteps <- NULL
        private$p <- NULL
        private$Pneg <- NULL
        private$Ppos <- NULL
        private$PPneg <- NULL
        private$PPpos <- NULL
      }
      # time
      t <- 100
      if(eps > 0 & rho > 0) { t <- -1.6*log(eps)/rho }
      if(t > 100) { t <- 100 }
      else if(t < 1) { t <- 1}
      tscale <- 1
      while(t > tscale) { tscale <- 10*tscale }
      t <- round(t/tscale,1)*tscale
      tfrom <- s
      tto <- s+t
      tby <- t/100
      tseq <- seq(from=tfrom,to=tto,by=tby)
      if(!private$vecs_equal(tseq,private$y_stoch_args$t))
      {
        private$y_stoch_args$t <- tseq
        private$G <- NULL
        private$Gteps <- NULL
        private$H2 <- NULL
        private$H2teps <- NULL
        private$p <- NULL
        private$Pneg <- NULL
        private$Ppos <- NULL
        private$PPneg <- NULL
        private$PPpos <- NULL
      }
      # density
      mean <- mu+(-x-mu)*exp(-rho*t)
      pmax <- 2.5*private$OUPDensity(0,-x,t,mean,rho,mu,sigma,yby)
      pscale <- 0.01
      while(pmax > pscale) { pscale <- 10*pscale }
      pmax <- round(pmax/pscale,2)*pscale
      private$plot_info$plottype$pmax <- pmax
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
      k <- self$DecisionThreshold(plotit=FALSE)[[1]]
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
      if(!private$vecs_equal(xseq,private$x_stoch_args$x))
      {
        private$x_stoch_args$x <- xseq
        if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(NULL,xseq,NULL,NULL,"A") }
        private$Oneg <- NULL
        private$Opos <- NULL
        private$Oshatneg <- NULL
        private$Oshatpos <- NULL
        private$Oscarfneg <- NULL
        private$Oscarfpos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
      }
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
        if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(sseq,NULL,NULL,NULL,"A") }
        private$Oneg <- NULL
        private$Opos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
      }
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
        zby <- z/5
        if(k < x) { zfrom <- k-3*zby }
        else { zfrom <- k-7*zby }
        zto <- zfrom+10*zby
      }
      else
      {
        z <- 2*abs(sigma)
        if(z < 1) { z <- 1}
        zscale <- 1
        while(z > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,1)*zscale
        zby <- z/5
        zscale <- 1
        while(abs(z) > zscale) { zscale <- 10*zscale }
        z <- round(z/zscale,2)*zscale
        zfrom <- z-5*zby
        zto <- zfrom+10*zby
      }
      zseq <- seq(from=zfrom,to=zto,by=zby)
      if(!private$vecs_equal(zseq,private$t_stoch_args$z))
      {
        private$t_stoch_args$z <- zseq
        private$tmodes <- NULL
        private$tmedians <- NULL
        private$tmeans <- NULL
        private$tvariances <- NULL
        private$tpercentiles <- NULL
        private$pt <- NULL
        private$Pt <- NULL
      }
      # time
      tpercentiles <- self$PassageTimePercentiles(plotit=FALSE)[[2]]
      tuppers <- tpercentiles[[2]]
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
      if(!private$vecs_equal(tseq,private$t_stoch_args$t))
      {
        private$t_stoch_args$t <- tseq
        private$ptx <- NULL
        private$pt <- NULL
        private$Ptx <- NULL
        private$Pt <- NULL
      }
      # density
      tmode <- self$PassageTimeMode(plotit=FALSE)[[1]]
      if(is.finite(tmode))
      {
        ptmax <- 1.2*private$OUPPassageTimeDensity(s,x,tmode,k,omega,rho,mu,sigma,tby)
        pscale <- 0.01
        while(ptmax > pscale) { pscale <- 10*pscale }
        ptmax <- round(ptmax/pscale,2)*pscale
        private$plot_info$plottype$ptmax <- ptmax
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
        decision <- private$KOneg
        if(is.null(decision))
        {
          k <- mu
          x <- k+sigma*(2/rho)^0.5
        }
        else
        {
          x <- self$DecisionThreshold(phi=1,plotit=FALSE)[[1]]
          k <- self$DecisionThreshold(phi=-1,plotit=FALSE)[[1]]
        }
      }
      else
      {
        decision <- private$KOpos
        if(is.null(decision))
        {
          k <- mu
          x <- k-sigma*(2/rho)^0.5
        }
        else
        {
          x <- self$DecisionThreshold(phi=-1,plotit=FALSE)[[1]]
          k <- self$DecisionThreshold(phi=1,plotit=FALSE)[[1]]
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
      if(k != private$t_stoch_args$k | !private$vecs_equal(zseq,private$t_stoch_args$z))
      {
        private$t_stoch_args$k <- k
        private$t_stoch_args$z <- zseq
        private$tmode <- NULL
        private$tmodes <- NULL
        private$tmedian <- NULL
        private$tmedians <- NULL
        private$tmean <- NULL
        private$tmeans <- NULL
        private$tvariance <- NULL
        private$tvariances <- NULL
        private$tpercentile <- NULL
        private$tpercentiles <- NULL
        private$ptx <- NULL
        private$pt <- NULL
        private$Ptx <- NULL
        private$Pt <- NULL
      }
      if(x != private$t_stoch_args$x)
      {
        private$t_stoch_args$x <- x
        private$tmode <- NULL
        private$tmedian <- NULL
        private$tmean <- NULL
        private$tvariance <- NULL
        private$tpercentile <- NULL
        private$ptx <- NULL
        private$Ptx <- NULL
      }
      # z state
      if(!private$vecs_equal(xyzseq,private$z_stoch_args$z))
      {
        private$z_stoch_args$z <- xyzseq
        private$g <- NULL
        private$h2 <- NULL
      }
      # y state
      if(!private$vecs_equal(xyzseq,private$y_stoch_args$y))
      {
        private$y_stoch_args$y <- xyzseq
        private$p <- NULL
        private$Pneg <- NULL
        private$Ppos <- NULL
        private$PPneg <- NULL
        private$PPpos <- NULL
      }
      if(psi != private$psiphi) {  private$y_stoch_args$psi <- private$psiphi }
      # x state
      if(!private$vecs_equal(xyzseq,private$x_stoch_args$x))
      {
        private$x_stoch_args$x <- xyzseq
        if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(NULL,xyzseq,NULL,NULL,"A") }
        private$Oneg <- NULL
        private$Opos <- NULL
        private$Oshatneg <- NULL
        private$Oshatpos <- NULL
        private$Oscarfneg <- NULL
        private$Oscarfpos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
      }
      if(phi != private$psiphi)
      {
        private$x_stoch_args$phi <- private$psiphi
        private$Oneg <- NULL
        private$Opos <- NULL
        private$Oshatneg <- NULL
        private$Oshatpos <- NULL
        private$Oscarfneg <- NULL
        private$Oscarfpos <- NULL
        private$KOneg <- NULL
        private$KOpos <- NULL
        private$BCneg <- NULL
        private$BCpos <- NULL
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
        z_stoch_args=private$z_stoch_args,
        y_stoch_args=private$y_stoch_args,
        x_stoch_args=private$x_stoch_args,
        t_stoch_args=private$t_stoch_args))
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
                not_equal <- FALSE
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
            t_stoch_args=private$t_stoch_args)))
        n <- n+1
      }

      return(n)
    },
    #' @description
    #' Replace current arguments from undo list
    #' @return c(index of this argument set, number of argument sets)
    undo_undo = function()
    {
      n <- length(private$undo_args)
      this_Ix <- private$undoIx
      these_undo_args <- private$undo_args[[this_Ix]]
      private$oup_params <- these_undo_args[[1]]
      rho <- private$oup_params$rho
      mu <- private$oup_params$mu
      sigma <- private$oup_params$sigma
      if(!is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"A") }
      private$z_stoch_args <- these_undo_args[[2]]
      private$y_stoch_args <- these_undo_args[[3]]
      private$x_stoch_args <- these_undo_args[[4]]
      s <- private$x_stoch_args$s
      x <- private$x_stoch_args$x
      r <- private$x_stoch_args$r
      phi <- private$x_stoch_args$phi
      if(!is.null(private$OUP)) { private$OUP$send_x_stoch_args(s,x,r,phi,"A") }
      private$t_stoch_args <- these_undo_args[[5]]
      next_Ix <- this_Ix+1
      if(next_Ix > n) { next_Ix <- 1 }
      private$undoIx <- next_Ix
      private$g <- NULL
      private$h2 <- NULL
      private$G <- NULL
      private$Gteps <- NULL
      private$H2 <- NULL
      private$H2teps <- NULL
      private$p <- NULL
      private$Pneg <- NULL
      private$Ppos <- NULL
      private$PPneg <- NULL
      private$PPpos <- NULL
      private$Oneg <- NULL
      private$Opos <- NULL
      private$Oshatneg <- NULL
      private$Oshatpos <- NULL
      private$Oscarfneg <- NULL
      private$Oscarfpos <- NULL
      private$KOneg <- NULL
      private$KOpos <- NULL
      private$BCneg <- NULL
      private$BCpos <- NULL
      private$tmode <- NULL
      private$tmodes <- NULL
      private$tmedian <- NULL
      private$tmedians <- NULL
      private$tmean <- NULL
      private$tmeans <- NULL
      private$tvariance <- NULL
      private$tvariances <- NULL
      private$tpercentile <- NULL
      private$tpercentiles <- NULL
      private$ptx <- NULL
      private$pt <- NULL
      private$Ptx <- NULL
      private$Pt <- NULL

      return(c(this_Ix,n))
    },
    # public calculate methods ----
    #' @description
    #' Calculate, plot and return drifts
    #' @param z       vector of n states
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(g(1xn))
    Drift = function(z=NULL,rho=NULL,mu=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_z_stoch_args(z)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      z <- private$z_stoch_args[[1]]
      # calculate ----
      drifts <- private$g
      if(is.null(drifts))
      {
        drifts <- -rho*(z-mu)
        private$g <- drifts
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDrift(NULL)) }

      return(list(g=drifts))
    },
    #' @description
    #' Calculate, plot and return diffusions
    #' @param z       vector of n states
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(h2(1xn))
    Diffusion = function(z=NULL,sigma=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(NULL,NULL,sigma)
      self$set_z_stoch_args(z)
      sigma <- private$oup_params[[3]]
      z <- private$z_stoch_args[[1]]
      # calculate ----
      diffusions <- private$h2
      if(is.null(diffusions))
      {
        n <- length(z)
        diffusions <- rep(sigma^2,n)
        private$h2 <- diffusions
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDiffusion(NULL,NULL)) }

      return(list(h2=diffusions))
    },
    #' @description
    #' Calculate and plot means
    #' @param t       vector of m times s<=t<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(G(mx1))
    Mean = function(t=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,mu,NULL)
      self$set_y_stoch_args(t,NULL,s,x,NULL,NULL)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      # calculate ----
      means <- private$G
      if(is.null(means))
      {
        means <- mu+(x-mu)*exp(-rho*(t-s))
        private$G <- means
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotMean(NULL,NULL,NULL)) }

      return(list(G=means))
    },
    #' @description
    #' Calculate and plot time for mean to converge
    #' @param s       initial time -inf<s<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param eps     proportion remaining after convergence 0<=eps<=1
    #' @param plotit  TRUE or FALSE
    #' @return list(Gteps)
    MeanToConverge = function(s=NULL,rho=NULL,eps=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,NULL,NULL)
      self$set_y_stoch_args(NULL,NULL,s,NULL,NULL,eps)
      rho <- private$oup_params[[1]]
      s <- private$y_stoch_args[[3]]
      eps <- private$y_stoch_args[[6]]
      # calculate ----
      timeeps <- private$Gteps
      if(is.null(timeeps))
      {
        timeeps <- s
        if(eps <= 0 | rho < 0.0000000001) { timeeps <- Inf }
        else if(eps < 1) { timeeps <- s-log(eps)/rho }
        private$Gteps <- timeeps
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotMeanToConverge(NULL)) }

      return(list(Gteps=timeeps))
    },
    #' @description
    #' Calculate and plot variances
    #' @param t       vector of m times s<=t<inf
    #' @param s       initial time -inf<s<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(H2(mx1))
    Variance = function(t=NULL,s=NULL,rho=NULL,sigma=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,NULL,sigma)
      self$set_y_stoch_args(t,NULL,s,NULL,NULL,NULL)
      rho <- private$oup_params[[1]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      # calculate ----
      variances <- private$H2
      if(is.null(variances))
      {
        if(rho < 0.0000000001) { variances <- sigma^2*(t-s) }
        else { variances <- sigma^2/(2*rho)*(1-exp(-2*rho*(t-s))) }
        private$H2 <- variances
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotVariance(NULL,NULL,NULL)) }

      return(list(H2=variances))
    },
    #' @description
    #' Calculate and plot time for variance to converge
    #' @param s       initial time -inf<s<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param eps     proportion remaining after convergence 0<=eps<=1
    #' @param plotit  TRUE or FALSE
    #' @return list(H2teps)
    VarianceToConverge = function(s=NULL,rho=NULL,eps=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_oup_params(rho,NULL,NULL)
      self$set_y_stoch_args(NULL,NULL,s,NULL,NULL,eps)
      rho <- private$oup_params[[1]]
      s <- private$y_stoch_args[[3]]
      eps <- private$y_stoch_args[[6]]
      # calculate ----
      timeeps <- private$H2teps
      if(is.null(timeeps))
      {
        timeeps <- s
        if(eps <= 0 | rho < 0.0000000001) { timeeps <- Inf }
        else if(eps < 1) { timeeps <- s-log(eps)/(2*rho) }
        private$H2teps <- timeeps
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotVarianceToConverge(NULL)) }

      return(list(H2teps=timeeps))
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
    #' @param plotit  TRUE or FALSE
    #' @return list(p(mxn))
    Density = function(t=NULL,y=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      densities <- private$p
      if(is.null(densities))
      {
        m <- length(t)
        n <- length(y)
        if(n > 1) { dy <- (y[n]-y[1])/(n-1) }
        else { dy <- 0.1 }
        densities <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            densities[i,j] <- private$OUPDensity(s,x,t[i],y[j],rho,mu,sigma,dy)
          }
        }
        private$p <- densities
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDensity(NULL,NULL,NULL)) }

      return(list(p=densities))
    },
    #' @description
    #' Calculate and plot transition probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param y       vector of nstates -inf<y<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param plotit  TRUE or FALSE
    #' @return list(P(mxn))
    Probability = function(t=NULL,y=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,psi=NULL,plotit=TRUE)
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
      # calculate ----
      if(psi <= 0) { probabilities <- private$Pneg }
      else { probabilities <- private$Ppos}
      if(is.null(probabilities))
      {
        m <- length(t)
        n <- length(y)
        probabilities <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            probabilities[i,j] <- private$OUPProbability(s,x,t[i],y[j],rho,mu,sigma,psi)
          }
        }
        if(psi <= 0) { private$Pneg <- probabilities }
        else { private$Ppos <- probabilities }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotProbability(NULL,NULL)) }

      return(list(P=probabilities))
    },
    #' @description
    #' Calculate and plot double integrals of transition densities
    #' @param t       vector of m times s<=t<inf
    #' @param y       vector of n states -inf<y<inf
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param psi     <=0 for integral -inf to y, >0 for integral y to inf
    #' @param plotit  TRUE or FALSE
    #' @return list(PP(mxn))
    DoubleIntegral = function(t=NULL,y=NULL,s=NULL,x=NULL,rho=NULL,mu=NULL,sigma=NULL,psi=NULL,plotit=TRUE)
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
      # calculate ----
      if(psi <= 0 ) { doubleintegrals <- private$PPneg }
      else { doubleintegrals <- private$PPpos }
      if(is.null(doubleintegrals))
      {
        m <- length(t)
        n <- length(y)
        doubleintegrals <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            doubleintegrals[i,j] <- private$OUPDoubleIntegral(s,x,t[i],y[j],rho,mu,sigma,psi)
          }
        }
        if(psi <= 0) { private$PPneg <- doubleintegrals }
        else { private$PPpos <- doubleintegrals }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDoubleIntegral(NULL,NULL)) }

      return(list(PP=doubleintegrals))
    },
    #' @description
    #' Calculate and plot option prices
    #' @param s       vector of m times -inf<s<t
    #' @param x       vector of n states -inf<x<inf
    #' @param t       terminal time -inf<t<inf
    #' @param y       terminal state -inf<y<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param plotit  TRUE or FALSE
    #' @return list(O(mxn))
    Option = function(s=NULL,x=NULL,t=NULL,y=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,plotit=TRUE)
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
      # calculate ----
      if(phi <= 0) { options <- private$Oneg }
      else { options <- private$Opos }
      if(is.null(options))
      {
        m <- length(s)
        n <- length(x)
        options <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            options[i,j] <- private$OUPOption(s[i],x[j],t,y,rho,mu,sigma,r,phi,b,c)
          }
        }
        if(phi <= 0) { private$Oneg <- options }
        else { private$Opos <- options }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotOption(NULL,NULL)) }

      return(list(O=options))
    },
    #' @description
    #' Calculate and plot the envelope of option prices
    #' @param x       vector of n states -inf<x<inf
    #' @param y       terminal state -inf<y<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param plotit  TRUE or FALSE
    #' @return list(Ohat(1xn),shat(1xn))
    OptionEnvelope = function(x=NULL,y=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,plotit=TRUE)
    {
      # set / gets ----
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
      # calculate ----
      if(phi <= 0) { Oshat <- private$Oshatneg }
      else { Oshat <- private$Oshatpos }
      if(is.null(Oshat))
      {
        n <- length(x)
        Oshat <- matrix(0.0,2,n)
        env <- private$OUPOptionMax(y,y,rho,mu,sigma,r,phi,b,c,0,1100,0,1)
        tsguess <- env[2]
        m <- 1
        while(x[m] < y & m < n) { m <- m+1 }
        j <- m+1
        while(j > 1)
        {
          j <- j-1
          env <- private$OUPOptionMax(x[j],y,rho,mu,sigma,r,phi,b,c,env[2],1100,0,1)
          Oshat[1,j] <- env[1]
          Oshat[2,j] <- env[2]
        }
        env[2] <- tsguess
        j <- m
        while(j < n)
        {
          j <- j+1
          env <- private$OUPOptionMax(x[j],y,rho,mu,sigma,r,phi,b,c,env[2],1100,0,1)
          Oshat[1,j] <- env[1]
          Oshat[2,j] <- env[2]
        }
        if(phi <= 0) { private$Oshatneg <- Oshat }
        else { private$Oshatpos <- Oshat }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotOptionEnvelope(NULL,NULL)) }

      return(list(Ohat=Oshat[1,],shat=Oshat[2,]))
    },
    #' @description
    #' Calculate and plot the decision threshold
    #' @param x       vector of n states -inf<x<inf
    #' @param y       terminal state -inf<y<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for exit option, >0 for entry option
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param plotit  TRUE or FALSE
    #' @return list(k,O)
    DecisionThreshold = function(x=NULL,y=NULL,rho=NULL,mu=NULL,sigma=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,plotit=TRUE)
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
      # calculate ----
      if(phi <= 0)
      {
        decision <- private$KOneg
        if(is.null(decision))
        {
          decision <- private$OUPThresholdSearch(y,rho,mu,sigma,r,phi,b,c)
          private$KOneg <- decision
        }
      }
      else
      {
        decision <- private$KOpos
        if(is.null(decision))
        {
          decision <- private$OUPThresholdSearch(y,rho,mu,sigma,r,phi,b,c)
          private$KOpos <- decision
        }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotDecisionThreshold(NULL)) }

      return(list(k=decision[1],O=decision[2]))
    },
    #' @description
    #' Calculate and plot obligations and prohibitions, ie. benefit/cost analysis
    #' @param s       vector of m times -inf<s<t
    #' @param x       vector of n states -inf<x<inf
    #' @param t       terminal time -inf<t<inf
    #' @param y       terminal state -inf<y<inf
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param r       discount rate -inf<r<inf
    #' @param phi     <=0 for obligation, >0 for prohibition
    #' @param b       lump-sum benefit for entry option
    #' @param c       lump-sum cost for exit option
    #' @param plotit  TRUE or FALSE
    #' @return list(BC(mxn))
    Obligation = function(s=NULL,x=NULL,t=NULL,y=NULL,rho=NULL,mu=NULL,r=NULL,phi=NULL,b=NULL,c=NULL,plotit=TRUE)
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
      # calculate ----
      if(phi <= 0) { obligations <- private$BCneg }
      else { obligations <- private$BCpos }
      if(is.null(obligations))
      {
        m <- length(s)
        n <- length(x)
        obligations <- matrix(0.0,m,n)
        if(phi > 0) { ltgt <- -1 }
        else( ltgt <- 1)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            G <- mu+(x[j]-mu)*exp(-rho*(t-s[i]))
            obligations[i,j] <- ltgt*exp(-r*(t-s[i]))*(G-y+b+c)
          }
        }
        if(phi <= 0) { private$BCneg <- obligations }
        else { private$BCpos <- obligations }
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotObligation(NULL,NULL)) }

      return(list(BC=obligations))
    },
    #' @description
    #' Calculate and plot passage time modes
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(tmode(1xn),tmodes(mxn))
    PassageTimeMode = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      medians <- self$PassageTimeMedian(plotit=FALSE)
      tmedian <- medians[[1]]
      tmedians <- medians[[2]]
      tmode <- private$tmode
      if(is.null(tmode))
      {
        tmode <- private$OUPPassageTimeModeSearch(s,x,k,omega,rho,mu,sigma,tmedian)
        private$tmode <- tmode
      }
      tmodes <- private$tmodes
      if(is.null(tmodes))
      {
        n <- length(z)
        tmodes <- vector("double",n)
        for(j in 1:n)
        {
          tmodes[j] <- private$OUPPassageTimeModeSearch(s,z[j],k,omega,rho,mu,sigma,tmedians[j])
        }
        private$tmodes <- tmodes
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeModeMedianMean(NULL,NULL,NULL)) }

      return(list(tmode=tmode,tmodes=tmodes))
    },
    #' @description
    #' Calculate and plot passage time medians
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(tmedian(1xn),tmedians(mxn))
    PassageTimeMedian = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      tmedian <- private$tmedian
      if(is.null(tmedian))
      {
        tmedian <- private$OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,0.5)
        private$tmedian <- tmedian
      }
      tmedians <- private$tmedians
      if(is.null(tmedians))
      {
        n <- length(z)
        tmedians <- vector("double",n)
        for(j in 1:n)
        {
          tmedians[j] <- private$OUPPassageTimePctSearch(s,z[j],k,omega,rho,mu,sigma,0.5)
        }
        private$tmedians <- tmedians
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeModeMedianMean(NULL,NULL,NULL)) }

      return(list(tmedian=tmedian,tmedians=tmedians))
    },
    #' @description
    #' Calculate and plot passage time means
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(tmean(1xn),tmeans(mxn))
    PassageTimeMean = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      medians <- self$PassageTimeMedian(plotit=FALSE)
      tmedian <- medians[[1]]
      tmedians <- medians[[2]]
      tmean <- private$tmean
      if(is.null(tmean))
      {
        tmean <- private$OUPPassageTimeMeanIntegrate(s,x,k,omega,rho,mu,sigma,tmedian)
        private$tmean <- tmean
      }
      tmeans <- private$tmeans
      if(is.null(tmeans))
      {
        n <- length(z)
        tmeans <- vector("double",n)
        for(j in 1:n)
        {
          tmeans[j] <- private$OUPPassageTimeMeanIntegrate(s,z[j],k,omega,rho,mu,sigma,tmedians[j])
        }
        private$tmeans <- tmeans
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeModeMedianMean(NULL,NULL,NULL)) }

      return(list(tmean=tmean,tmeans=tmeans))
    },
    #' @description
    #' Calculate and plot passage time variances
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return (tvariance(1xn),tvariances(mxn))
    PassageTimeVariance = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      medians <- self$PassageTimeMedian(plotit=FALSE)
      tmedian <- medians[[1]]
      tmedians <- medians[[2]]
      means <- self$PassageTimeMean(plotit=FALSE)
      tmean <- means[[1]]
      tmeans <- means[[2]]
      tvariance <- private$tvariance
      if(is.null(tvariance))
      {
        tvariance <- private$OUPPassageTimeVarianceIntegrate(s,x,k,omega,rho,mu,sigma,tmedian,tmean)
        private$tvariance <- tvariance
      }
      tvariances <- private$tvariances
      if(is.null(tvariances))
      {
        n <- length(z)
        tvariances <- vector("double",n)
        for(j in 1:n)
        {
          tvariances[j] <- private$OUPPassageTimeVarianceIntegrate(s,z[j],k,omega,rho,mu,sigma,tmedians[j],tmeans[j])
        }
        private$tvariances <- tvariances
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeVariance(NULL,NULL)) }

      return(list(tvariance=tvariance,tvariances=tvariances))
    },
    #' @description
    #' Calculate and plot passage time percentiles
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param Ppct    passage time probability for a percentile 0.01<=Ppct<=0.99
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(tpercentile(1xn),tpercentiles(mxn))
    PassageTimePercentiles = function(k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,Ppct=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      tpercentile <- private$tpercentile
      if(is.null(tpercentile))
      {
        if(Ppct > 0.5)
        {
          tupper <- private$OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,Ppct)
          tlower <- private$OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,1-Ppct)
        }
        else
        {
          tupper <- private$OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,1-Ppct)
          tlower <- private$OUPPassageTimePctSearch(s,x,k,omega,rho,mu,sigma,Ppct)
        }
        tpercentile <- list(tlower=tlower,tupper=tupper)
        private$tpercentile <- tpercentile
      }
      tpercentiles <- private$tpercentiles
      if(is.null(tpercentiles))
      {
        n <- length(z)
        tlowers <- vector("double",n)
        tuppers <- vector("double",n)
        if(Ppct > 0.5)
        {
          for(j in 1:n)
          {
            tuppers[j] <- private$OUPPassageTimePctSearch(s,z[j],k,omega,rho,mu,sigma,Ppct)
            tlowers[j] <- private$OUPPassageTimePctSearch(s,z[j],k,omega,rho,mu,sigma,1-Ppct)
          }
        }
        else
        {
          for(j in 1:n)
          {
            tuppers[j] <- private$OUPPassageTimePctSearch(s,z[j],k,omega,rho,mu,sigma,1-Ppct)
            tlowers[j] <- private$OUPPassageTimePctSearch(s,z[j],k,omega,rho,mu,sigma,Ppct)
          }
        }
        tpercentiles <- list(tlowers=tlowers,tuppers=tuppers)
        private$tpercentiles <- tpercentiles
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimePercentiles(NULL,NULL,NULL)) }

      return(list(tpercentile=tpercentile,tpercentiles=tpercentiles))
    },
    #' @description
    #' Calculate and plot passage time densities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(ptx(1xn),pt(mxn))
    PassageTimeDensity = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) & is.null(pt))
      {
        m <- length(t)
        n <- length(z)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        ptx <- vector("double",m)
        pt <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          ptx[i] <- private$OUPPassageTimeDensity(s,x,t[i],k,omega,rho,mu,sigma,dt)
          for(j in 1:n)
          {
            pt[i,j] <- private$OUPPassageTimeDensity(s,z[j],t[i],k,omega,rho,mu,sigma,dt)
          }
        }
        private$ptx <- ptx
        private$pt <- pt
      }
      if(is.null(ptx))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        ptx <- vector("double",m)
        for(i in 1:m)
        {
          ptx[i] <- private$OUPPassageTimeDensity(s,x,t[i],k,omega,rho,mu,sigma,dt)
        }
        private$ptx <- ptx
      }
      if(is.null(pt))
      {
        m <- length(t)
        n <- length(z)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        pt <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            pt[i,j] <- private$OUPPassageTimeDensity(s,z[j],t[i],k,omega,rho,mu,sigma,dt)
          }
        }
        private$pt <- pt
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeDensity(NULL,NULL,NULL)) }

      return(list(ptx=ptx,pt=pt))
    },
    #' @description
    #' Calculate and plot passage time probabilities
    #' @param t       vector of m times s<=t<inf
    #' @param k       decision threshold -inf<k<int
    #' @param s       initial time -inf<s<inf
    #' @param x       initial state -inf<x<inf
    #' @param z       vector of alternate initial states -inf<x<inf
    #' @param omega   degree of irreversibility 0<=omega<=1
    #' @param rho     rate parameter 0<=rho<inf
    #' @param mu      location parameter -inf<mu<inf
    #' @param sigma   scale parameter -inf<sigma<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(Ptx(1xn),Pt(mxn))
    PassageTimeProbability = function(t=NULL,k=NULL,s=NULL,x=NULL,z=NULL,omega=NULL,rho=NULL,mu=NULL,sigma=NULL,plotit=TRUE)
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
      # calculate ----
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) & is.null(Pt))
      {
        m <- length(t)
        n <- length(z)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        Ptx <- vector("double",m)
        Pt <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          Ptx[i] <- private$OUPPassageTimeProbability(s,x,t[i],k,omega,rho,mu,sigma)
          for(j in 1:n)
          {
            Pt[i,j] <- private$OUPPassageTimeProbability(s,z[j],t[i],k,omega,rho,mu,sigma)
          }
        }
        private$Ptx <- Ptx
        private$Pt <- Pt
      }
      if(is.null(Ptx))
      {
        m <- length(t)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        Ptx <- vector("double",m)
        for(i in 1:m)
        {
          Ptx[i] <- private$OUPPassageTimeProbability(s,x,t[i],k,omega,rho,mu,sigma)
        }
        private$Ptx <- Ptx
      }
      if(is.null(Pt))
      {
        m <- length(t)
        n <- length(z)
        if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
        else { dt <- 0.05 }
        Pt <- matrix(0.0,m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            Pt[i,j] <- private$OUPPassageTimeProbability(s,z[j],t[i],k,omega,rho,mu,sigma)
          }
        }
        private$Pt <- Pt
      }
      # plot ----
      if(plotit == TRUE) { print(self$PlotPassageTimeProbability(NULL,NULL)) }

      return(list(Ptx=Ptx,Pt=Pt))
    },
    # public plot methods ----
    #' @description
    #' Plot drifts
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotDrift = function(title=NULL,xaxis=NULL,yaxis=NULL,zbeg=NULL,zend=NULL)
    {
      # get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      z <- private$z_stoch_args[[1]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      drift <- self$Drift(plotit=FALSE)[[1]]
      n <- length(z)
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        n <- length(z)
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
        add_trace(.,type="scatter",x=z,y=drift,name="<i>g</i>(<i>z</i>)",mode="lines",line=driftline) %>%
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
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotDiffusion = function(type=NULL,title=NULL,xaxis=NULL,yaxis=NULL,zbeg=NULL,zend=NULL)
    {
      # get ----
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      z <- private$z_stoch_args[[1]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      drift <- self$Drift(plotit=FALSE)[[1]]
      diffusion <- self$Diffusion(plotit=FALSE)[[1]]
      sqrt <- diffusion^0.5
      driftplus <- drift+sqrt
      driftminus <- drift-sqrt
      n <- length(z)
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        drift <- drift[Ixbeg:Ixend]
        diffusion <- diffusion[Ixbeg:Ixend]
        driftplus <- driftplus[Ixbeg:Ixend]
        driftminus <- driftminus[Ixbeg:Ixend]
        n <- length(z)
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
      if(type < 2.5)
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
          add_trace(.,type="scatter",x=z,y=drift,name="<i>g</i>(<i>z</i>)",mode="lines",line=driftline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=z,y=driftplus,name="<i>g</i>(<i>z</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=z,y=driftminus,name="<i>g</i>(<i>z</i>)&plusmn;<i>h</i>",mode="lines",line=diffusionline,legendgroup="g+h",showlegend=FALSE,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions) %>%
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
          add_trace(.,type="scatter",x=z,y=diffusion,name="<i>h</i><sup>2</sup>",mode="lines",line=diffusionline) %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot means
    #' @param type  = 3 for 2D, 4 and 5 for 3D
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
      self$set_plot_info(type,pmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      type <- private$plot_info$plottype$type
      pmax <- private$plot_info$plottype$pmax
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
      m <- length(t)
      n <- length(y)
      if(n > 1) { dy <- (y[n]-y[1])/(n-1) }
      else { dy <- 1 }
      means <- self$Mean(plotit=FALSE)[[1]]
      densities <- self$Density(plotit=FALSE)[[1]]
      probabilities <- self$Probability(plotit=FALSE)[[1]]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,]
        probabilities <- probabilities[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend]
        probabilities <- probabilities[,Ixbeg:Ixend]
        n <- length(y)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),")",esml)
        if(is.null(title)) { title <- "Mean" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      #2D
      if(type < 3.5)
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
        # OUP_A_Mean2D
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        meanline <- list(color=cyn$d,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Mean2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline) %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
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
        # OUP_A_Mean3Ddensity
        if(type < 4.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
          pmeans <- vector("double",m)
          coordinatemeans <- vector("double",m)
          coordinatepmeans <- vector("double",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            pmeans[i] <- private$OUPDensity(s,x,t[i],means[i],rho,mu,sigma,dy)
            coordinatepmeans[i] <- paste(sep="","<i>p</i>(<i>t,G</i>)=",format(pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>p</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>y</i>=",format(y[j],digits=4)) }
          }
          pmesh <- MeshWall(means,t,pmeans)
          ygap <- 0.03*(y[n]-y[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(y[1]-ygap,y[n]+ygap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,range=c(0-0.03*pmax,pmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=blu$e,font=list(color=blu$b))
          pmeanline <- list(color=cyn$d,width=8)
          densityline <- list(color=blu$d,width=6)
          gradient <- list(c(0,cyn$b),c(1,cyn$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Mean3Ddensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=pmeans,name="<i>p</i>(<i>t,G</i>)",mode="lines",line=pmeanline,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=pmesh$xvertex,y=pmesh$yvertex,z=pmesh$zvertex,i=pmesh$ivertex,j=pmesh$jvertex,k=pmesh$kvertex,intensity=pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans",visible="legendonly",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          q <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly")
          i <- i+dt
          q <- q+1
          while(i <= m)
          {
            lineopacity <- lineopacity-0.05
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly",showlegend=FALSE)
            i <- i+dt
            q <- q+1
          }
        }
        # OUP_A_Mean3Dprobability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
          Pmeans <- vector("double",m)
          coordinatemeans <- vector("double",m)
          coordinatePmeans <- vector("double",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            Pmeans[i] <- private$OUPProbability(s,x,t[i],means[i],rho,mu,sigma,psi)
            coordinatePmeans[i] <- paste(sep="","<i>P</i>(<i>t,G</i>)=",format(Pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>y</i>=",y[j]) }
          }
          Pmesh <- MeshWall(means,t,Pmeans)
          ygap <- 0.03*(y[n]-y[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(y[1]-ygap,y[n]+ygap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=cyn$c,linewidth=3,gridcolor=cyn$c,gridwidth=2,backgroundcolor=cyn$b,showbackground=floor,range=c(0-0.03,1+0.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=grn$e,font=list(color=grn$b))
          Pmeanline <- list(color=cyn$d,width=8)
          probabilityline <- list(color=grn$d,width=6)
          gradient <- list(c(0,cyn$b),c(1,cyn$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Mean3Dprobability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=Pmeans,name="<i>P</i>(<i>t,G</i>)",mode="lines",line=Pmeanline,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Pmesh$xvertex,y=Pmesh$yvertex,z=Pmesh$zvertex,i=Pmesh$ivertex,j=Pmesh$jvertex,k=Pmesh$kvertex,intensity=Pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans",visible="legendonly",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          q <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly")
          i <- i+dt
          q <- q+1
          while(i <= m)
          {
            lineopacity <- lineopacity-0.05
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly",showlegend=FALSE)
            i <- i+dt
            q <- q+1
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the mean to converge
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotMeanToConverge = function(title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # get ----
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      eps <- private$y_stoch_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      cyn <- private$plot_colors$cyn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      means <- self$Mean(plotit=FALSE)[[1]]
      timeeps <- self$MeanToConverge(plotit=FALSE)[[1]]
      meanteps <- x
      if(rho > 0) { meanteps <- mu+(x-mu)*exp(-rho*(timeeps-s)) }
      m <- length(t)
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        m <- length(t)
      }
      # plot ----
      # OUP_A_MeanToConverge2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",<i>x</i>",bsym,"=",esym,format(x,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>e</i>=",esym,format(eps,digits=4),")",esml)
        if(is.null(title)) { title <- "Mean To Converge" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "<i>G</i>(<i>t</i>|<i>s,x</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      meanline <- list(color=cyn$d,width=4)
      meantepsline <- list(dash="dot",color=cyn$d,width=4)
      tepsline <- list(dash="dot",color=cyn$d,width=4)
      horzline <- list(color=gry$d,width=1)
      if(x < mu) { tepsG <- list(x=timeeps,y=meanteps,text=paste(sep="","&nbsp;<i>t</i>",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>G</i>",bsym,"=",esym,format(meanteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="top",align="left",showarrow=FALSE) }
      else { tepsG <- list(x=timeeps,y=meanteps,text=paste(sep="","&nbsp;<i>t</i>",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>G</i>",bsym,"=",esym,format(meanteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE) }
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_MeanToConverge2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),name=" ",mode="lines",line=horzline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(x,x),name=" ",mode="lines",line=horzline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(t[1],timeeps),y=c(meanteps,meanteps),name=paste(sep="","<i>G</i>(<i>t</i><sub>",eps,"</sub>)"),mode="lines",line=meantepsline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(timeeps,timeeps),y=c(mu,x),name=paste(sep="","<i>t</i><sub>",eps,"</sub>"),mode="lines",line=tepsline,hoverinfo="x+y") %>%
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,annotations=tepsG,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot variances
    #' @param type  = 2 and 3 for 2D, 4 and 5 for 3D
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
      self$set_plot_info(type,pmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      type <- private$plot_info$plottype$type
      pmax <- private$plot_info$plottype$pmax
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
      m <- length(t)
      n <- length(y)
      if(n > 1) { dy <- (y[n]-y[1])/(n-1) }
      else { dy <- 1 }
      means <- self$Mean(plotit=FALSE)[[1]]
      variances <- self$Variance(plotit=FALSE)[[1]]
      sqrts <- variances^0.5
      meansplus <- means+sqrts
      meansminus <- means-sqrts
      densities <- self$Density(plotit=FALSE)[[1]]
      probabilities <- self$Probability(plotit=FALSE)[[1]]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        means <- means[Ixbeg:Ixend]
        variances <- variances[Ixbeg:Ixend]
        meansplus <- meansplus[Ixbeg:Ixend]
        meansminus <- meansminus[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,]
        probabilities <- probabilities[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend]
        probabilities <- probabilities[,Ixbeg:Ixend]
        n <- length(y)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- "Variance" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # 2D
      if(type < 3.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
        lookdown <- list(text=xaxis)
        #OUP_A_Variance2DG
        if(type < 2.5)
        {
          if(is.null(yaxis)) { yaxis <- "<i>G</i>(<i>t</i>|<i>s,x</i>)&plusmn;<i>H</i>(<i>t</i>|<i>s</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          varianceline <- list(color=mgn$d,width=4)
          meanline <- list(color=cyn$d,width=4)
          horzline <- list(color=gry$d,width=1)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Diffusion2DG")
          ymax <- max(meansplus)
          ymin <- min(meansminus)
          legendy <- (mu-ymin) / (ymax-ymin)
          if(legendy > 0.6) { legendy <- 0.15 }
          else if(legendy < 0.4) { legendy <- 0.95 }
          else if(x > mu) { legendy <- legendy-0.05 }
          else { legendy <- legendy+0.15 }
          legendpos <- list(orientation="h",x=1.0,y=legendy,xanchor="right")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter",x=t,y=means,name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=meansplus,name="<i>G</i>(<i>t</i>)&plusmn;<i>H</i>(<i>t</i>)",mode="lines",line=varianceline,legendgroup="G+H",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=meansminus,name="<i>G</i>(<i>t</i>)&plusmn;<i>H</i>(<i>t</i>)",mode="lines",line=varianceline,legendgroup="G+H",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(mu,mu),mode="lines",line=horzline,showlegend=FALSE,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_Variance2D
        else
        {
          if(is.null(yaxis)) { yaxis <- "<i>H</i><sup>2</sup>(<i>t</i>|<i>s</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
          varianceline <- list(color=mgn$d,width=4)
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Variance2D")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter",x=t,y=variances,name="<i>H</i><sup>2</sup>(<i>t</i>)",mode="lines",line=varianceline) %>%
            config(.,toImageButtonOptions=imageoptions) %>%
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
        # OUP_A_Variance3Ddensity
        if(type < 4.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
          pmeansplus <- vector("double",m)
          pmeans <- vector("double",m)
          pmeansminus <- vector("double",m)
          coordinatemeansplus <- vector("double",m)
          coordinatemeans <- vector("double",m)
          coordinatemeansminus <- vector("double",m)
          coordinatepmeansplus <- vector("double",m)
          coordinatepmeans <- vector("double",m)
          coordinatepmeansminus <- vector("double",m)
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
          pmeshplus <- MeshWall(meansplus,t,pmeansplus)
          pmesh <- MeshWall(means,t,pmeans)
          pmeshminus <- MeshWall(meansminus,t,pmeansminus)
          ygap <- 0.03*(y[n]-y[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(y[1]-ygap,y[n]+ygap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,range=c(0-0.03*pmax,pmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=blu$e,font=list(color=blu$b))
          pplusline <- list(color=mgn$d,width=8)
          pmeanline <- list(color=cyn$d,width=8)
          pminusline <- list(color=mgn$d,width=8)
          densityline <- list(color=blu$d,width=6)
          gradientplus <- list(c(0,mgn$b),c(1,mgn$b))
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientminus <- list(c(0,mgn$b),c(1,mgn$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Variance3Ddensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)",mode="lines",line=plusline,hoverinfo="text",text=coordinatemeansplus) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)",mode="lines",line=minusline,hoverinfo="text",text=coordinatemeansminus) %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=pmeansplus,name="<i>p</i>(<i>G</i>+<i>H</i>)",mode="lines",line=pplusline,hoverinfo="text",text=coordinatepmeansplus,legendgroup="pmeansplus",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=pmeshplus$xvertex,y=pmeshplus$yvertex,z=pmeshplus$zvertex,i=pmeshplus$ivertex,j=pmeshplus$jvertex,k=pmeshplus$kvertex,intensity=pmeshplus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientplus,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatepmeansplus,legendgroup="pmeansplus",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=pmeans,name="<i>p</i>(<i>t,G</i>)",mode="lines",line=pmeanline,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=pmesh$xvertex,y=pmesh$yvertex,z=pmesh$zvertex,i=pmesh$ivertex,j=pmesh$jvertex,k=pmesh$kvertex,intensity=pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatepmeans,legendgroup="pmeans",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=pmeansminus,name="<i>p</i>(<i>G</i>-<i>H</i>)",mode="lines",line=pminusline,hoverinfo="text",text=coordinatepmeansminus,legendgroup="pmeansminus",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=pmeshminus$xvertex,y=pmeshminus$yvertex,z=pmeshminus$zvertex,i=pmeshminus$ivertex,j=pmeshminus$jvertex,k=pmeshminus$kvertex,intensity=pmeshminus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientminus,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatepmeansminus,legendgroup="pmeansminus",visible="legendonly",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          q <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],name="<i>p</i>(<i>t,y</i>)",mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly")
          i <- i+dt
          q <- q+1
          while(i <= m)
          {
            lineopacity <- lineopacity-0.05
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="p",visible="legendonly",showlegend=FALSE)
            i <- i+dt
            q <- q+1
          }
        }
        # OUP_A_Variance3Dprobability
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
          coordinatemeansplus <- vector("double",m)
          coordinatemeans <- vector("double",m)
          coordinatemeansminus <- vector("double",m)
          coordinatePmeansplus <- vector("double",m)
          coordinatePmeans <- vector("double",m)
          coordinatePmeansminus <- vector("double",m)
          coordinates <- matrix("",m,n)
          for(i in 1:m)
          {
            coordinatemeansplus[i] <- paste(sep="","<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)=",format(meansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeans[i] <- paste(sep="","<i>G</i>(<i>t</i>)=",format(means[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatemeansminus[i] <- paste(sep="","<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)=",format(meansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4))
            coordinatePmeansplus[i] <- paste(sep="","<i>P</i>(<i>t,G+<i>H</i>)=",format(Pmeansplus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G+H</i>=",format(meansplus[i],digits=4))
            coordinatePmeans[i] <- paste(sep="","<i>P</i>(<i>t,G</i>)=",format(Pmeans[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G</i>=",format(means[i],digits=4))
            coordinatePmeansminus[i] <- paste(sep="","<i>P</i>(<i>t,G-<i>H</i>)=",format(Pmeansminus[i],digits=4),"<br><i>t</i>=",format(t[i],digits=4),"<br><i>G-H</i>=",format(meansminus[i],digits=4))
            for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
          }
          Pmeshplus <- MeshWall(meansplus,t,Pmeansplus)
          Pmesh <- MeshWall(means,t,Pmeans)
          Pmeshminus <- MeshWall(meansminus,t,Pmeansminus)
          ygap <- 0.03*(y[n]-y[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(y[1]-ygap,y[n]+ygap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=mgn$c,linewidth=3,gridcolor=mgn$c,gridwidth=2,backgroundcolor=mgn$b,showbackground=floor,range=c(0-0.03,1+0.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=grn$e,font=list(color=grn$b))
          Pplusline <- list(color=mgn$d,width=8)
          Pmeanline <- list(color=cyn$d,width=8)
          Pminusline <- list(color=mgn$d,width=8)
          probabilityline <- list(color=grn$d,width=6)
          gradientplus <- list(c(0,mgn$b),c(1,mgn$b))
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientminus <- list(c(0,mgn$b),c(1,mgn$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Variance3Dprobability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)+<i>H</i>(<i>t</i>)",mode="lines",line=plusline,hoverinfo="text",text=coordinatemeansplus) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=rep(0,m),name="<i>G</i>(<i>t</i>)-<i>H</i>(<i>t</i>)",mode="lines",line=minusline,hoverinfo="text",text=coordinatemeansminus) %>%
            add_trace(.,type="scatter3d",x=meansplus,y=t,z=Pmeansplus,name="<i>P</i>(<i>t,G</i>+<i>H</i>)",mode="lines",line=Pplusline,hoverinfo="text",text=coordinatePmeansplus,legendgroup="Pmeansplus",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Pmeshplus$xvertex,y=Pmeshplus$yvertex,z=Pmeshplus$zvertex,i=Pmeshplus$ivertex,j=Pmeshplus$jvertex,k=Pmeshplus$kvertex,intensity=Pmeshplus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientplus,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePmeansplus,legendgroup="Pmeansplus",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=means,y=t,z=Pmeans,name="<i>P</i>(<i>t,G</i>)",mode="lines",line=Pmeanline,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Pmesh$xvertex,y=Pmesh$yvertex,z=Pmesh$zvertex,i=Pmesh$ivertex,j=Pmesh$jvertex,k=Pmesh$kvertex,intensity=Pmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePmeans,legendgroup="Pmeans",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=meansminus,y=t,z=Pmeansminus,name="<i>P</i>(<i>t,G</i>-<i>H</i>)",mode="lines",line=Pminusline,hoverinfo="text",text=coordinatePmeansminus,legendgroup="Pmeansminus",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Pmeshminus$xvertex,y=Pmeshminus$yvertex,z=Pmeshminus$zvertex,i=Pmeshminus$ivertex,j=Pmeshminus$jvertex,k=Pmeshminus$kvertex,intensity=Pmeshminus$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientminus,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePmeansminus,legendgroup="Pmeansminus",visible="legendonly",showlegend=FALSE)
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          lineopacity <- 1
          i <- 1
          q <- 1
          for(j in 1:n) { tt[j] <- t[i] }
          fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],name="<i>P</i>(<i>t,y</i>)",mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly")
          i <- i+dt
          q <- q+1
          while(i <= m)
          {
            lineopacity <- lineopacity-0.05
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,],legendgroup="P",visible="legendonly",showlegend=FALSE)
            i <- i+dt
            q <- q+1
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the variance to converge
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @return plot
    PlotVarianceToConverge = function(title=NULL,xaxis=NULL,yaxis=NULL,tbeg=NULL,tend=NULL)
    {
      # get ----
      rho <- private$oup_params[[1]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      s <- private$y_stoch_args[[3]]
      eps <- private$y_stoch_args[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      variances <- self$Variance(plotit=FALSE)[[1]]
      asymvar <- sigma^2/(2*rho)
      timeeps <- self$VarianceToConverge(plotit=FALSE)[[1]]
      varianceteps <- Inf
      if(rho > 0) { varianceteps <- sigma^2/(2*rho)*(1-exp(-2*rho*(timeeps-s))) }
      m <- length(t)
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        variances <- variances[Ixbeg:Ixend]
        m <- length(t)
      }
      # plot ----
      # OUP_A_VarianceToConverge2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>s</i>",bsym,"=",esym,format(s,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",",bsym,"<i>e</i>=",esym,format(eps,digits=4),")",esml)
        if(is.null(title)) { title <- "Variance To Converge" }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>t</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "<i>H</i><sup>2</sup>(<i>t</i>|<i>s</i>)" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero")
      varianceline <- list(color=mgn$d,width=4)
      variancetepsline <- list(dash="dot",color=mgn$d,width=4)
      tepsline <- list(dash="dot",color=mgn$d,width=4)
      horzline <- list(color=gry$d,width=1)
      tepsH2 <- list(x=timeeps,y=varianceteps,text=paste(sep="","&nbsp;<i>t</i>",bsym,"=",esym,format(timeeps,digits=4),"<br>&nbsp;<i>H</i><sup>2</sup>",bsym,"=",esym,format(varianceteps,digits=4)),xref="x",yref="y",xanchor="left",yanchor="top",align="left",showarrow=FALSE)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_VarianceToConverge2D")
      fig <- plot_ly()  %>%
        add_trace(.,type="scatter",x=t,y=variances,name="<i>H</i><sup>2</sup>(<i>t</i>)",mode="lines",line=varianceline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(t[1],t[m]),y=c(asymvar,asymvar),name=" ",mode="lines",line=horzline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(t[1],timeeps),y=c(varianceteps,varianceteps),name=paste(sep="","<i>H</i><sup>2</sup>(t<sub>",eps,"</sub>)"),mode="lines",line=variancetepsline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(timeeps,timeeps),y=c(0,asymvar),name=paste(sep="","<i>t</i><sub>",eps,"</sub>"),mode="lines",line=tepsline,hoverinfo="x+y") %>%
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,annotations=tepsH2,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot transition densities
    #' @param type  = 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type,pmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      type <- private$plot_info$plottype$type
      pmax <- private$plot_info$plottype$pmax
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(t)
      n <- length(y)
      densities <- self$Density(plotit=FALSE)[[1]]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        densities <- densities[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        densities <- densities[,Ixbeg:Ixend]
        n <- length(y)
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
      if(type < 2.5)
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
          fig <- add_trace(fig,type="scatter",x=y,y=densities[i,],name=paste(sep="","<i>p</i>(",t[i],"<i>,y</i>)"),mode="lines",line=densityline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
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
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p</i>(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>p</i>(<i>t,y</i>)=",format(densities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        hover <- list(bgcolor=blu$e,font=list(color=blu$b))
        if(x < mu) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(x == mu) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(pmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(0-0.03*pmax,pmax),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_Density3DSurface and OUP_A_Density3DSurfaceScatter
        if(type < 4.5)
        {
          gradient <- list(c(0,blu$c),c(1,blu$c))
          rise <- list(x=0,y=-300,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=y,y=t,z=densities,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
          # OUP_A_Density3DSurface
          if(type < 3.5) { imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Density3DSurface") }
          # OUP_A_Density3DSurfaceScatter
          else
          {
            densityline <- list(color=blu$e,width=8)
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Density3DSurfaceScatter")
            tt <- vector("double",n)
            dt <- as.integer((m-1)/10)
            if(dt < 1) { dt <- 1 }
            i <- 1
            while(i <= m)
            {
              for(j in 1:n) { tt[j] <- t[i] }
              fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,hoverinfo="text",text=coordinates[i,])
              i <- i+dt
            }
          }
        }
        # OUP_A_Density3DScatter
        else
        {
          densityline <- list(color=blu$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Density3DScatter")
          fig <- plot_ly()
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          i <- 1
          while(i <= m)
          {
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=densities[i,],mode="lines",line=densityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,])
            i <- i+dt
            lineopacity <- lineopacity-0.05
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot transition probabilities
    #' @param type  = 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(t)
      n <- length(y)
      probabilities <- self$Probability(plotit=FALSE)[[1]]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        probabilities <- probabilities[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        probabilities <- probabilities[,Ixbeg:Ixend]
        n <- length(y)
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
      if(type < 2.5)
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
          fig <- add_trace(fig,type="scatter",x=y,y=probabilities[i,],name=paste(sep="","<i>P</i>(",t[i],"<i>,y</i>)"),mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
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
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P</i>(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P</i>(<i>t,y</i>)=",format(probabilities[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        hover <- list(bgcolor=grn$e,font=list(color=grn$b))
        if(psi > 0) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(psi == 0) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview, zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_Probability3DSurface and OUP_A_Probability3DSurfaceScatter
        if(type < 4.5)
        {
          gradient <- list(c(0,grn$c),c(1,grn$c))
          rise <- list(x=0,y=-300,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=y,y=t,z=probabilities,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
          # OUP_A_Probability3DSurface
          if(type < 3.5) { imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Probability3DSurface") }
          # OUP_A_Probability3DSurfaceScatter
          else
          {
            probabilityline <- list(color=grn$e,width=8)
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Probability3DSurfaceScatter")
            tt <- vector("double",n)
            dt <- as.integer((m-1)/10)
            if(dt < 1) { dt <- 1 }
            i <- 1
            while(i <= m)
            {
              for(j in 1:n) { tt[j] <- t[i] }
              fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,hoverinfo="text",text=coordinates[i,])
              i <- i+dt
            }
          }
        }
        # OUP_A_Probability3DScatter
        else
        {
          probabilityline <- list(color=grn$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Probability3DScatter")
          fig <- plot_ly()
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          i <- 1
          while(i <= m)
          {
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=probabilities[i,],mode="lines",line=probabilityline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,])
            i <- i+dt
            lineopacity <- lineopacity-0.05
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot double integrals of transition densities
    #' @param type  = 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$y_stoch_args[[1]]
      y <- private$y_stoch_args[[2]]
      s <- private$y_stoch_args[[3]]
      x <- private$y_stoch_args[[4]]
      psi <- private$y_stoch_args[[5]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(t)
      n <- length(y)
      doubleintegrals <- self$DoubleIntegral(plotit=FALSE)[[1]]
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(y,ybeg,yend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        y <- y[Ixbeg:Ixend]
        doubleintegrals <- doubleintegrals[,Ixbeg:Ixend]
        n <- length(y)
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
      if(type < 2.5)
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
          fig <- add_trace(fig,type="scatter",x=y,y=doubleintegrals[i,],name=paste(sep="","\u2119(",t[i],"<i>,y</i>)"),mode="lines",line=doubleintegralline,opacity=lineopacity,hoverinfo="x+y")
          i <- i+dt
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
        if(is.null(xaxis)) { xaxis <- "<i>y</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "\u2119(<i>t,y</i>|<i>s,x</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","\u2119(<i>t,y</i>)=",format(doubleintegrals[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>y</i>=",y[j]) }
        }
        hover <- list(bgcolor=red$e,font=list(color=red$b))
        if(psi > 0) { spy <- list(x=0.8,y=-2.3,z=0.5) }
        else if(psi == 0) { spy <- list(x=0,y=-2.4,z=0.5) }
        else { spy <- list(x=-0.8,y=-2.3,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_DoubleIntegral3DSurface and OUP_A_DoubleIntegral3DSurfaceScatter
        if(type < 4.5)
        {
          gradient <- list(c(0,red$c),c(1,red$c))
          rise <- list(x=0,y=-300,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=y,y=t,z=doubleintegrals,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
          # OUP_A_DoubleIntegral3DSurface
          if(type < 3.5) { imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_DoubleIntegral3DSurface") }
          # OUP_A_DoubleIntegral3DSurfaceScatter
          else
          {
            doubleintegralline <- list(color=red$e,width=8)
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_DoubleIntegral3DSurfaceScatter")
            tt <- vector("double",n)
            dt <- as.integer((m-1)/10)
            if(dt < 1) { dt <- 1 }
            i <- 1
            while(i <= m)
            {
              for(j in 1:n) { tt[j] <- t[i] }
              fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=doubleintegrals[i,],mode="lines",line=doubleintegralline,hoverinfo="text",text=coordinates[i,])
              i <- i+dt
            }
          }
        }
        # OUP_A_DoubleIntegral3DScatter
        else
        {
          doubleintegralline <- list(color=red$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_DoubleIntegral3DScatter")
          fig <- plot_ly()
          tt <- vector("double",n)
          dt <- as.integer((m-1)/10)
          if(dt < 1) { dt <- 1 }
          i <- 1
          while(i <= m)
          {
            for(j in 1:n) { tt[j] <- t[i] }
            fig <- add_trace(fig,type="scatter3d",x=y,y=tt,z=doubleintegrals[i,],mode="lines",line=doubleintegralline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,])
            i <- i+dt
            lineopacity <- lineopacity-0.05
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
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
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
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
      options <- self$Option(plotit=FALSE)[[1]]
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
        syms <- paste(sep="","(<i>t</i>",bsml,bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),")",esml)
        if(is.null(title)) { title <- "Option" }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Option2D
      if(type < 2.5)
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
          fig <- add_trace(fig,type="scatter",x=x,y=options[i,],name=paste(sep="","\uD835\uDD46(",s[i],"<i>,x</i>)"),mode="lines",line=optionline,opacity=lineopacity,hoverinfo="x+y")
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
        if(is.null(zaxis)) { zaxis <- "\uD835\uDD46(<i>s,x</i>|<i>t,y</i>)" }
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",s[i],"<br><i>x</i>=",x[j])
          }
        }
        hover <- list(bgcolor=red$e,font=list(color=red$b))
        if(phi > 0) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_Option3DSurface and OUP_A_Option3DSurfaceScatter
        if(type < 4.5)
        {
          gradient <- list(c(0,red$c),c(1,red$c))
          rise <- list(x=0,y=-300,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=x,y=s,z=options,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
          # OUP_A_Option3DSurface
          if(type < 3.5) { imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Option3DSurface") }
          # OUP_A_Option3DSurfaceScatter
          else
          {
            optionline <- list(color=red$e,width=8)
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Option3DSurfaceScatter")
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
        # OUP_A_Option3DScatter
        else
        {
          optionline <- list(color=red$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Option3DScatter")
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
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot the option envelope
    #' @param type  = 3 for 2D, 4, 5 and 6 for 3D
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
      t <- private$x_stoch_args[[3]]
      y <- private$x_stoch_args[[4]]
      r <- private$x_stoch_args[[5]]
      phi <- private$x_stoch_args[[6]]
      b <- private$x_stoch_args[[7]]
      c <- private$x_stoch_args[[8]]
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
      Oshat <- self$OptionEnvelope(plotit=FALSE)
      Ohat <- Oshat[[1]]
      shat <- t-Oshat[[2]]
      options <- self$Option(plotit=FALSE)[[1]]
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
        Ohat <- Ohat[Ixbeg:Ixend]
        shat <- shat[Ixbeg:Ixend]
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
        syms <- paste(sep="",bsml,"(<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),")",esml)
        if(is.null(title)) { title="Option Envelope"  }
      }
      else if(is.null(title)) { title="" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Envelope2D
      if(type < 3.5)
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
        Ohatline <- list(color=red$d,width=4)
        terminalline <- list(color=gry$c,width=4)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Envelope2D")
        fig <- plot_ly() %>%
          add_trace(.,type="scatter",x=x,y=terminal,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
          add_trace(.,type="scatter",x=x,y=Ohat,name="\u00D4(<i>x</i>)",mode="lines",line=Ohatline,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>x</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>s</i>" }
        if(is.null(zaxis)) { zaxis <- "\u00D4(<i>x</i>|<i>y</i>)" }
        Ohold <- vector("double",n)
        Oexercise <- vector("double",n)
        coordinateshat <- vector("double",n)
        coordinates <- matrix("",m,n)
        for(j in 1:n)
        {
          if(shat[j] == t)
          {
            Ohold[j] <- NA
            Oexercise[j] <- Ohat[j]
          }
          else
          {
            Ohold[j] <- Ohat[j]
            Oexercise[j] <- NA
          }
          coordinateshat[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(Ohat[j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          for(i in 1:m) { coordinates[i,j] <- paste(sep="","\uD835\uDD46(<i>s,x</i>)=",format(options[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
        }
        hover <- list(bgcolor=red$e,font=list(color=red$b))
        if(phi > 0) { spy <- list(x=-0.4,y=-2.3,z=0.1) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.1) }
        else { spy <- list(x=0.4,y=-2.3,z=0.1) }
        xview <- list(title=xaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(x[1],x[n]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$a,showbackground=walls,range=c(s[m],s[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=red$c,linewidth=3,gridcolor=red$c,gridwidth=2,backgroundcolor=red$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview, zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        Oholdline <- list(color=red$e,width=10)
        Oexerciseline <- list(color=red$e,width=8)
        Ohatline <- list(dash="dash",color=ylw$e,width=8)
        Oscarfline <- list(color=red$a,width=10)
        gradient <- list(c(0,red$c),c(1,red$c))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        fig <- plot_ly() %>%
          add_trace(.,type="surface",x=x,y=s,z=options,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates)
      # OUP_A_Envelope3D
        if(type < 4.5)
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_OptionEnvelope3D")
          fig <- add_trace(fig,type="scatter3d",x=x,y=shat,z=Ohat,mode="lines",line=Ohatline,hoverinfo="text",text=coordinateshat) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=Ohold,mode="lines",line=Oholdline,hoverinfo="text",text=coordinateshat) %>%
          add_trace(.,type="scatter3d",x=x,y=shat,z=Oexercise,mode="lines",line=Oexerciseline,hoverinfo="text",text=coordinateshat)
        }
      # OUP_A_Envelope3DScarf and OUP_A_Envelope3DMittens
        else if(type < 6.5)
        {
          if(phi <= 0) {
            Oscarf <- private$Oscarfneg
            if(is.null(Oscarf)) { Oscarf <- private$OUPOptionDtZero() }
            sscarf <- t-private$sscarfneg
          }
          else
          {
            Oscarf <- private$Oscarfpos
            if(is.null(Oscarf)) { Oscarf <- private$OUPOptionDtZero() }
            sscarf <- t-private$sscarfpos
          }
          Oscarf <- Oscarf[,Ixbeg:Ixend]
          sscarf <- sscarf[,Ixbeg:Ixend]
          coordinatesscarf1 <- vector("double",n)
          coordinatesscarf2 <- vector("double",n)
          for(j in 1:n)
          {
            coordinatesscarf1[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(Oscarf[1,j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
            coordinatesscarf2[j] <- paste(sep="","\u00D4(<i>x</i>)=",format(Oscarf[2,j],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          }
      # OUP_A_Envelope3DScarf
          if(type < 5.5)
          {
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_OptionEnvelope3DScarf")
            fig <- add_trace(fig,type="scatter3d",x=x,y=sscarf[1,],z=Oscarf[1,],mode="lines",line=Oscarfline,hoverinfo="text",text=coordinatesscarf1) %>%
              add_trace(.,type="scatter3d",x=x,y=sscarf[2,],z=Oscarf[2,],mode="lines",line=Oscarfline,hoverinfo="text",text=coordinatesscarf2)
          }
      # OUP_A_Envelope3DMittens
          else
          {
            for(j in 1:n)
            {
              if(shat[j] < t) { Oscarf[2,j] <- NA }
            }
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_OptionEnvelope3DMittens")
            fig <- add_trace(fig,type="scatter3d",x=x,y=shat,z=Ohat,mode="lines",line=Ohatline,hoverinfo="text",text=coordinateshat) %>%
              add_trace(.,type="scatter3d",x=x,y=shat,z=Ohold,mode="lines",line=Oholdline,hoverinfo="text",text=coordinateshat) %>%
              add_trace(.,type="scatter3d",x=x,y=shat,z=Oexercise,mode="lines",line=Oexerciseline,hoverinfo="text",text=coordinateshat) %>%
              add_trace(.,type="scatter3d",x=x,y=sscarf[1,],z=Oscarf[1,],mode="lines",line=Oscarfline,hoverinfo="text",text=coordinatesscarf1) %>%
              add_trace(.,type="scatter3d",x=x,y=sscarf[2,],z=Oscarf[2,],mode="lines",line=Oscarfline,hoverinfo="text",text=coordinatesscarf2)
          }
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
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
      n <- length(x)
      decision <- self$DecisionThreshold(plotit=FALSE)
      k <- decision[[1]]
      O <- decision[[2]]
      Ohat <- self$OptionEnvelope(plotit=FALSE)[[1]]
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        Ohat <- Ohat[Ixbeg:Ixend]
        n <- length(x)
      }
      # plot ----
      # OUP_A_Decision2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),")",esml)
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
      Ohatline <- list(color=red$d,width=4)
      terminalline <- list(color=gry$c,width=4)
      oline <- list(dash="dot",color=red$d,width=4)
      kline <- list(dash="dot",color=red$d,width=4)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_Decision2D")
      fig <- plot_ly() %>%
        add_trace(.,type="scatter",x=x,y=terminal,name="<i>V</i>(<i>x</i>)",mode="lines",line=terminalline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=x,y=Ohat,name="\u00D4(<i>x</i>)",mode="lines",line=Ohatline,hoverinfo="x+y") %>%
        add_trace(.,type="scatter",x=c(k,k),y=c(0,O),name="<i>k</i>",mode="lines",line=kline,hoverinfo="x+y")
      if(phi > 0)
      {
        fig <- add_trace(fig,type="scatter",x=c(x[n],k),y=c(O,O),name="\u00D4",mode="lines",line=oline,hoverinfo="x+y")
        KO <- list(x=k,y=O,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(O,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
      }
      else
      {
        fig <- add_trace(fig,type="scatter",x=c(x[1],k),y=c(O,O),name="\u00D4",mode="lines",line=oline,hoverinfo="x+y")
        KO <- list(x=k,y=O,text=paste(sep="","<i>k</i>",bsym,"=",esym,format(k,digits=4),"<br>\u00D4",bsym,"=",esym,format(O,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
      }
      fig <- config(fig,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,annotations=KO,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot obligations
    #' @param type  = 2 for 2D, 3, 4, 5 and 6 for 3D
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
      self$set_plot_info(type)
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
      obligations <- self$Obligation(plotit=FALSE)[[1]]
      Oshat <- self$OptionEnvelope(plotit=FALSE)
      Ohat <- Oshat[[1]]
      shat <- t-Oshat[[2]]
      Inx <- xedni(s,sbeg,send)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg < m | Ixend > 1)
      {
        s <- s[Ixend:Ixbeg]
        obligations <- obligations[Ixend:Ixbeg,]
        m <- length(s)
      }
      Inx <- index(x,xbeg,xend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        x <- x[Ixbeg:Ixend]
        Ohat <- Ohat[Ixbeg:Ixend]
        shat <- shat[Ixbeg:Ixend]
        Obenv <- Obenv[Ixbeg:Ixend]
        obligations <- obligations[,Ixbeg:Ixend]
        n <- length(x)
      }
      # plot ----
      if(labels == TRUE)
      {
        bsml <- "<span style='font-size: 10pt;'>"
        esml <- "</span>"
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        syms <- paste(sep="","(<i>t</i>",bsml,bsym,"=",esym,format(t,digits=4),",<i>y</i>",bsym,"=",esym,format(y,digits=4),",",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",<i>r</i>",bsym,"=",esym,format(r,digits=4),",",bsym,"<i>f</i>=",esym,format(phi,digits=4),")",esml)
        if(is.null(title))
        {
          if(phi > 0) { title <- "Prohibition" }
          else { title <- "Obligation" }
        }
      }
      else if(is.null(title)) { title <- "" }
      lookup <- list(text=title,yref="container",y=0.95)
      # OUP_A_Obligation2D
      if(type < 2.5)
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
            fig <- add_trace(fig,type="scatter",x=x,y=obligations[i,],name=paste(sep="","\u2102(",s[i],"<i>,x</i>)"),mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="x+y")
            i <- i+ds
            lineopacity <- lineopacity-0.05
          }
        }
        else
        {
          while(i <= m)
          {
            fig <- add_trace(fig,type="scatter",x=x,y=obligations[i,],name=paste(sep="","\uD835\uDD39(",s[i],"<i>,x</i>)"),mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="x+y")
            i <- i+ds
            lineopacity <- lineopacity-0.05
          }
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
        if(is.null(zaxis))
        {
          if(phi>0) { zaxis <- "\u2102(<i>s,x</i>|<i>t,y</i>)" }
          else { zaxis <- "\uD835\uDD39(<i>s,x</i>|<i>t,y</i>)" }
        }
        zeroes <- matrix(0.0,m,n)
        coordinates <- matrix("",m,n)
        coordinateszero <- matrix("",m,n)
        for(i in 1:m)
        {
          for(j in 1:n)
          {
            if(phi > 0) { coordinates[i,j] <- paste(sep="","\u2102(<i>s,x</i>)=",format(obligations[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
            else { coordinates[i,j] <- paste(sep="","\uD835\uDD39(<i>s,x</i>)=",format(obligations[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4)) }
            coordinateszero[i,j] <- paste(sep="","0=",format(zeroes[i,j],digits=4),"<br><i>s</i>=",format(s[i],digits=4),"<br><i>x</i>=",format(x[j],digits=4))
          }
        }
        hover <- list(bgcolor=ylw$e,font=list(color=ylw$b))
        if(phi > 0) { spy <- list(x=0.6,y=-2.3,z=0.3) }
        else if(phi == 0) { spy <- list(x=0,y=-2.2,z=0.3) }
        else { spy <- list(x=-0.6,y=-2.3,z=0.3) }
        obgradient <- list(c(0,ylw$c),c(1,ylw$b))
        zgradient <- list(c(0,ylw$b),c(1,ylw$b))
        rise <- list(x=0,y=-300,z=0)
        shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
        xview <- list(title=xaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$a,showbackground=walls,range=c(x[1],x[n]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$a,showbackground=walls,range=c(s[m],s[1]),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=ylw$c,linewidth=3,gridcolor=ylw$c,gridwidth=2,backgroundcolor=ylw$b,showbackground=floor,tickmode="auto",nticks=5,zeroline=TRUE,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_Obligation3DSurface and OUP_A_Obligation3DSurfaceScatter
        if(type < 4.5)
        {
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=x,y=s,z=obligations,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
            add_trace(.,type="surface",x=x,y=s,z=zeroes,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=zgradient,reversescale=reverse,opacity=0.2,hoverinfo="text",text=coordinateszero)
          # OUP_A_Obligation3DSurface
          if(type < 3.5) { imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Obligation3DSurface") }
          # OUP_A_Obligation3DSurfaceScatter
          else
          {
            obligationline <- list(color=ylw$e,width=8)
            imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Obligation3DSurfaceScatter")
            ss <- vector("double",n)
            ds <- as.integer((m-1)/10)
            if(ds < 1) { ds <- 1 }
            i <- 1
            while(i <= m)
            {
              for(j in 1:n) { ss[j] <- s[i] }
              fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=obligations[i,],mode="lines",line=obligationline,hoverinfo="text",text=coordinates[i,])
              i <- i+ds
            }
          }
        }
        # OUP_A_Obligation3DScatter
        else if(type < 5.5)
        {
          obligationline <- list(color=ylw$d,width=6)
          lineopacity <- 1
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Obligation3DScatter")
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=x,y=s,z=zeroes,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateszero)
          ss <- vector("double",n)
          ds <- as.integer((m-1)/10)
          if(ds < 1) { ds <- 1 }
          i <- 1
          while(i <= m)
          {
            for(j in 1:n) { ss[j] <- s[i] }
            fig <- add_trace(fig,type="scatter3d",x=x,y=ss,z=obligations[i,],mode="lines",line=obligationline,opacity=lineopacity,hoverinfo="text",text=coordinates[i,])
            i <- i+ds
            lineopacity <- lineopacity-0.05
          }
        }
        # OUP_A_Obligation3DEnvelope
        else
        {
          obhat <- vector("double",n)
          obhold <- vector("double",n)
          obexercise <- vector("double",n)
          ophat <- vector("double",n)
          ophold <- vector("double",n)
          opexercise <- vector("double",n)
          zero <- rep(0,n)
          coordinatesobhat <- vector("double",n)
          coordinatesophat <- vector("double",n)
          coordinateszero <- vector("double",n)
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
              ophat[j] <- Ohat[j]+obhat[j]
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
                if(j > 1)
                {
                  obhat[j] <- obhat[j-1]
                  ophat[j] <- ophat[j-1]
                }
              }
              else
              {
                obhold[j] <- obhat[j]
                obexercise[j] <- NA
                ophold[j] <- ophat[j]
                opexercise[j] <- NA
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
                if(j < n)
                {
                  obhat[j] <- obhat[j+1]
                  ophat[j] <- ophat[j+1]
                }
              }
              else
              {
                obhold[j] <- obhat[j]
                obexercise[j] <- NA
                ophold[j] <- ophat[j]
                opexercise[j] <- NA
              }
              j <- j-1
            }
          }
          holdmesh <- MeshVertical(x,shat,ophold,obhold)
          exercisemesh <- MeshVertical(x,shat,opexercise,obexercise)
          obholdline <- list(color=red$c,width=10)
          obexerciseline <- list(color=red$c,width=8)
          obhatline <- list(dash="dash",color=ylw$c,width=8)
          opholdline <- list(color=red$c,width=10)
          opexerciseline <- list(color=red$c,width=8)
          ophatline <- list(dash="dash",color=ylw$c,width=8)
          zeroline <- list(color=red$c,width=8)
          opgradient <- list(c(0,red$c),c(1,red$c))
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_Obligation3DEnvelope")
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=x,y=s,z=obligations,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=obgradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=obhat,mode="lines",line=obhatline,hoverinfo="text",text=coordinatesobhat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=obhold,mode="lines",line=obholdline,hoverinfo="text",text=coordinatesobhat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=obexercise,mode="lines",line=obexerciseline,hoverinfo="text",text=coordinatesobhat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=ophat,mode="lines",line=ophatline,hoverinfo="text",text=coordinatesophat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=ophold,mode="lines",line=opholdline,hoverinfo="text",text=coordinatesophat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=opexercise,mode="lines",line=opexerciseline,hoverinfo="text",text=coordinatesophat) %>%
            add_trace(.,type="scatter3d",x=x,y=shat,z=zero,mode="lines",line=zeroline,hoverinfo="text",text=coordinateszero) %>%
            add_trace(.,type="mesh3d",x=holdmesh$xvertex,y=holdmesh$yvertex,z=holdmesh$zvertex,i=holdmesh$ivertex,j=holdmesh$jvertex,k=holdmesh$kvertex,intensity=holdmesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=opgradient,reversescale=reverse,opacity=0.7,hoverinfo="text",text=coordinatesophat) %>%
            add_trace(.,type="mesh3d",x=exercisemesh$xvertex,y=exercisemesh$yvertex,z=exercisemesh$zvertex,i=exercisemesh$ivertex,j=exercisemesh$jvertex,k=exercisemesh$kvertex,intensity=exercisemesh$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=opgradient,reversescale=reverse,opacity=0.7,hoverinfo="text",text=coordinatesophat)
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time modes, medians and means
    #' @param type  = 1, 2, 3 and 4 for 2D, 5 and 6 for 3D
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
    PlotPassageTimeModeMedianMean = function(type=NULL,xaxis=NULL,yaxis=NULL,zaxis=NULL,ptmax=NULL,title=NULL,tbeg=NULL,tend=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      self$set_plot_info(type,NULL,ptmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      type <- private$plot_info$plottype$type
      ptmax <- private$plot_info$plottype$ptmax
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
      m <- length(t)
      n <- length(z)
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      tmedian <- private$tmedian
      tmedians <- private$tmedians
      if(is.null(tmedian) | is.null(tmedians))
      {
        medians <- self$PassageTimeMedian(plotit=FALSE)
        tmedian <- medians[[1]]
        tmedians <- medians[[2]]
      }
      tmode <- private$tmode
      tmodes <- private$tmodes
      if(is.null(tmode) | is.null(tmodes))
      {
        modes <- self$PassageTimeMode(plotit=FALSE)
        tmode <- modes[[1]]
        tmodes <- modes[[2]]
      }
      tmean <- private$tmean
      tmeans <- private$tmeans
      if(is.null(tmean) | is.null(tmeans))
      {
        means <- self$PassageTimeMean(plotit=FALSE)
        tmean <- means[[1]]
        tmeans <- means[[2]]
      }
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) | is.null(pt))
      {
        densities <- self$PassageTimeDensity(plotit=FALSE)
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) | is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(plotit=FALSE)
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        tmedians <- tmedians[Ixbeg:Ixend]
        tmodes <- tmodes[Ixbeg:Ixend]
        tmeans <- tmeans[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend]
        Pt <- Pt[,Ixbeg:Ixend]
        n <- length(z)
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
      if(type < 2.5)
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
        meandotline <- list(color=cyn$d,dash="dot",width=3)
        mediandashline <- list(color=grn$d,dash="dash",width=3)
        mediandotline <- list(color=grn$d,dash="dot",width=3)
        modedashline <- list(color=blu$d,dash="dot",width=3)
        modedotline <- list(color=blu$d,dash="dot",width=3)
        fig <- plot_ly()
        # OUP_A_PassageTimeModeMedianMean2Dpt
        if(type < 1.5)
        {
          if(is.null(yaxis)) { yaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          lookleft <- list(text=yaxis)
          horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
          if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",rangemode="tozero") }
          else
          {
            mindensity <- 9
            for(i in 1:m)
            {
              if(ptx[i] < mindensity) { mindensity <- ptx[i] }
            }
            if( mindensity > 0) { mindensity <- 0 }
            vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(mindensity,ptmax))
          }
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2Dpt")
          legendpos <- list(x=1.0,y=1.0,xanchor="right",tracegroupgap=0)
          if(m > 1) { dt = t[2]-t[1] }
          else { dt = 0.1 }
          if(rho > 0)
          {
            ptmean <- private$OUPPassageTimeDensity(s,x,tmean,k,omega,rho,mu,sigma,dt)
            fig <- add_trace(fig,type="scatter",x=c(tmean,tmean),y=c(0,ptmean),name="mean",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(s,tmean),y=c(ptmean,ptmean),name="<i>p<sub>t</sub></i>(mean)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y")
          }
          ptmedian <- private$OUPPassageTimeDensity(s,x,tmedian,k,omega,rho,mu,sigma,dt)
          ptmode <- private$OUPPassageTimeDensity(s,x,tmode,k,omega,rho,mu,sigma,dt)
          fig <-  add_trace(fig,type="scatter",x=c(tmedian,tmedian),y=c(0,ptmedian),name="median",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tmedian),y=c(ptmedian,ptmedian),name="<i>p<sub>t</sub></i>(median)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tmode,tmode),y=c(0,ptmode),name="mode",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tmode),y=c(ptmode,ptmode),name="<i>p<sub>t</sub></i>(mode)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
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
          if(rho > 0)
          {
            Ptmean <- private$OUPPassageTimeProbability(s,x,tmean,k,omega,rho,mu,sigma)
            fig <- add_trace(fig,type="scatter",x=c(tmean,tmean),y=c(0,Ptmean),name="mean",mode="lines",line=meandashline,legendgroup="mean",hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(s,tmean),y=c(Ptmean,Ptmean),name="<i>P<sub>t</sub></i>(mean)",mode="lines",line=meandotline,legendgroup="mean",showlegend=FALSE,hoverinfo="x+y")
          }
          Ptmedian <- private$OUPPassageTimeProbability(s,x,tmedian,k,omega,rho,mu,sigma)
          Ptmode <- private$OUPPassageTimeProbability(s,x,tmode,k,omega,rho,mu,sigma)
          fig <-  add_trace(fig,type="scatter",x=c(tmedian,tmedian),y=c(0,Ptmedian),name="median",mode="lines",line=mediandashline,legendgroup="median",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tmedian),y=c(Ptmedian,Ptmedian),name="<i>P<sub>t</sub></i>(median)",mode="lines",line=mediandotline,legendgroup="median",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tmode,tmode),y=c(0,Ptmode),name="mode",mode="lines",line=modedashline,legendgroup="mode",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tmode),y=c(Ptmode,Ptmode),name="<i>P<sub>t</sub></i>(mode)",mode="lines",line=modedotline,legendgroup="mode",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=probabilityline,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 2D continued
      else if(type < 4.5)
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
        # OUP_A_PassageTimeModeMedianMean2Dx
        if(type < 3.5)
        {
          meandashline <- list(color=cyn$d,dash="longdash",width=3)
          meandotline <- list(color=cyn$d,dash="dot",width=3)
          mediandashline <- list(color=grn$d,dash="dash",width=3)
          mediandotline <- list(color=grn$d,dash="dot",width=3)
          modedashline <- list(color=blu$d,dash="dot",width=3)
          modedotline <- list(color=blu$d,dash="dot",width=3)
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
            xmmm <- list(x=x,y=legendy,text=paste(sep="","mean",bsym,"=",esym,format(tmean,digits=4),"<br>median",bsym,"=",esym,format(tmedian,digits=4),"<br>mode",bsym,"=",esym,format(tmode,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE,hoverinfo="x+y")
          }
          else
          {
            if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tmean,tmean),name="mean",mode="lines",line=meandashline,hoverinfo="x+y") }
            fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tmedian,tmedian),name="median",mode="lines",line=mediandashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(tmode,tmode),name="mode",mode="lines",line=modedashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","mean",bsym,"=",esym,format(tmean,digits=4),"<br>median",bsym,"=",esym,format(tmedian,digits=4),"<br>mode",bsym,"=",esym,format(tmode,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
          }
          if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=z,y=tmeans,name="mean(<i>z</i>)",mode="lines",line=meanline,hoverinfo="x+y") }
          fig <- add_trace(fig,type="scatter",x=z,y=tmedians,name="median(<i>z</i>)",mode="lines",line=medianline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tmodes,name="mode(<i>z</i>)",mode="lines",line=modeline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=xmmm,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimeModeMedianMean2D
        else
        {
          legendx <- (k-z[1])/(z[n]-z[1])
          if(legendx < 0.15) { legendx <- 0.15 }
          else if(legendx > 0.85) { legendx <- 0.85 }
          legendpos <- list(x=legendx,y=1.0,xanchor="center",yanchor="top")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeModeMedianMean2D")
          if(is.finite(tmean)) { fig <- add_trace(fig,type="scatter",x=z,y=tmeans,name="mean(<i>z</i>)",mode="lines",line=meanline,hoverinfo="x+y") }
          fig <-  add_trace(fig,type="scatter",x=z,y=tmedians,name="median(<i>z</i>)",mode="lines",line=medianline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tmodes,name="mode(<i>z</i>)",mode="lines",line=modeline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(z[1],z[n]),y=c(s,s),mode="lines",line=horzline,showlegend=FALSE,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        meanline <- list(color=cyn$e,width=8)
        medianline <- list(color=grn$e,width=8)
        modeline <- list(color=blu$e,width=8)
        legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
        # OUP_A_PassageTimeModeMedianMean3Ddensity
        if(type < 5.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          if(k > mu) { spy <- list(x=2.35,y=0.85,z=0.5) }
          else if(k == mu) { spy <- list(x=0,y=2.55,z=0.5) }
          else { spy <- list(x=-2.35,y=0.85,z=0.5) }
          mindensity <- 9
          ptmeans <- vector("double",n)
          ptmedians <- vector("double",n)
          ptmodes <- vector("double",n)
          coordinatemeans <- vector("double",n)
          coordinatemedians <- vector("double",n)
          coordinatemodes <- vector("double",n)
          coordinateptmeans <- vector("double",n)
          coordinateptmedians <- vector("double",n)
          coordinateptmodes <- vector("double",n)
          coordinatex <- vector("double",m)
          coordinates <- matrix("",m,n)
          for(j in 1:n)
          {
            if(rho > 0)
            {
              coordinatemeans[j] <- paste(sep="","mean=",format(tmeans[j],digits=4),"<br><i>x</i>=",z[j])
              ptmeans[j] <- private$OUPPassageTimeDensity(s,z[j],tmeans[j],k,omega,rho,mu,sigma,dt)
              coordinateptmeans[j] <- paste(sep="","<i>p<sub>t</sub></i>(mean)=",format(ptmeans[j],digits=4),"<br><i>t</i>=",tmeans[j],"<br><i>x</i>=",z[j])
              if(ptmeans[j] < mindensity) { mindensity <- ptmeans[j] }
            }
            coordinatemedians[j] <- paste(sep="","median=",format(tmedians[j],digits=4),"<br><i>x</i>=",z[j])
            ptmedians[j] <- private$OUPPassageTimeDensity(s,z[j],tmedians[j],k,omega,rho,mu,sigma,dt)
            coordinateptmedians[j] <- paste(sep="","<i>p<sub>t</sub></i>(median)=",format(ptmedians[j],digits=4),"<br><i>t</i>=",tmedians[j],"<br><i>x</i>=",z[j])
            if(ptmedians[j] < mindensity) { mindensity <- ptmedians[j] }
            coordinatemodes[j] <- paste(sep="","mode=",format(tmodes[j],digits=4),"<br><i>x</i>=",z[j])
            ptmodes[j] <- private$OUPPassageTimeDensity(s,z[j],tmodes[j],k,omega,rho,mu,sigma,dt)
            coordinateptmodes[j] <- paste(sep="","<i>p<sub>t</sub></i>(mode)=",format(ptmodes[j],digits=4),"<br><i>t</i>=",tmodes[j],"<br><i>x</i>=",z[j])
            if(ptmodes[j] < mindensity) { mindensity <- ptmodes[j] }
          }
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
          ptmeshmeans <- MeshWall(z,tmeans,ptmeans)
          ptmeshmedians <- MeshWall(z,tmedians,ptmedians)
          ptmeshmodes <- MeshWall(z,tmodes,ptmodes)
          zgap <- 0.03*(z[n]-z[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(z[1]-zgap,z[n]+zgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(mindensity-0.03*(ptmax-mindensity),ptmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=blu$e,font=list(color=blu$b))
          ptmeanline <- list(color=cyn$d,width=8)
          ptmedianline <- list(color=grn$d,width=8)
          ptmodeline <- list(color=blu$d,width=8)
          ptxline <- list(color=gry$d,width=6)
          ptline <- list(color=blu$d,width=6)
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientmedian <- list(c(0,grn$b),c(1,grn$b))
          gradientmode <- list(c(0,blu$b),c(1,blu$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeModeMedianMean3Ddensity")
          fig <- plot_ly()
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=rep(0,n),name="mean(<i>z</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=rep(0,n),name="median(<i>z</i>)",mode="lines",line=medianline,hoverinfo="text",text=coordinatemedians) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=rep(0,n),name="mode(<i>z</i>)",mode="lines",line=modeline,hoverinfo="text",text=coordinatemodes)
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=ptmeans,name="<i>p<sub>t</sub></i>(mean)",mode="lines",line=ptmeanline,hoverinfo="text",text=coordinateptmeans,legendgroup="ptmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmeans$xvertex,y=ptmeshmeans$yvertex,z=ptmeshmeans$zvertex,i=ptmeshmeans$ivertex,j=ptmeshmeans$jvertex,k=ptmeshmeans$kvertex,intensity=ptmeshmeans$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateptmeans,legendgroup="ptmeans",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=ptmedians,name="<i>p<sub>t</sub></i>(median)",mode="lines",line=ptmedianline,hoverinfo="text",text=coordinateptmedians,legendgroup="ptmedians",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmedians$xvertex,y=ptmeshmedians$yvertex,z=ptmeshmedians$zvertex,i=ptmeshmedians$ivertex,j=ptmeshmedians$jvertex,k=ptmeshmedians$kvertex,intensity=ptmeshmedians$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmedian,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateptmedians,legendgroup="ptmedians",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=ptmodes,name="<i>p<sub>t</sub></i>(mode)",mode="lines",line=ptmodeline,hoverinfo="text",text=coordinateptmodes,legendgroup="ptmodes",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshmodes$xvertex,y=ptmeshmodes$yvertex,z=ptmeshmodes$zvertex,i=ptmeshmodes$ivertex,j=ptmeshmodes$jvertex,k=ptmeshmodes$kvertex,intensity=ptmeshmodes$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmode,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateptmodes,legendgroup="ptmodes",visible="legendonly",showlegend=FALSE)
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="pt",visible="legendonly")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="pt",visible="legendonly",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,visible="legendonly")
        }
        # OUP_A_PassageTimeModeMedianMean3Dprobability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          if(k > mu) { spy <- list(x=0.85,y=-2.35,z=0.5) }
          else if(k == mu) { spy <- list(x=0,y=2.55,z=0.5) }
          else { spy <- list(x=-0.85,y=-2.35,z=0.5) }
          Ptmeans <- vector("double",n)
          Ptmedians <- vector("double",n)
          Ptmodes <- vector("double",n)
          coordinatemeans <- vector("double",n)
          coordinatemedians <- vector("double",n)
          coordinatemodes <- vector("double",n)
          coordinatePtmeans <- vector("double",n)
          coordinatePtmedians <- vector("double",n)
          coordinatePtmodes <- vector("double",n)
          coordinatex <- vector("double",m)
          coordinates <- matrix("",m,n)
          kindex <- 0
          for(j in 1:n)
          {
            PInf <- private$OUPPassageTimeProbabilityInf(z[j],k,omega,rho,mu,sigma)
            if(rho > 0)
            {
              coordinatemeans[j] <- paste(sep="","mean=",format(tmeans[j],digits=4),"<br><i>x</i>=",z[j])
              if(sigma > 0) { Ptmeans[j] <- private$OUPPassageTimeProbability(s,z[j],tmeans[j],k,omega,rho,mu,sigma) }
              else { Ptmeans[j] <- 0.5*PInf }
              coordinatePtmeans[j] <- paste(sep="","<i>P<sub>t</sub></i>(mean)=",format(Ptmeans[j],digits=4),"<br><i>t</i>=",tmeans[j],"<br><i>x</i>=",z[j])
            }
            coordinatemedians[j] <- paste(sep="","median=",format(tmedians[j],digits=4),"<br><i>x</i>=",z[j])
            Ptmedians[j] <- 0.5*PInf
            coordinatePtmedians[j] <- paste(sep="","<i>P<sub>t</sub></i>(median)=",format(Ptmedians[j],digits=4),"<br><i>t</i>=",tmedians[j],"<br><i>x</i>=",z[j])
            coordinatemodes[j] <- paste(sep="","mode=",format(tmodes[j],digits=4),"<br><i>x</i>=",z[j])
            if( sigma > 0) { Ptmodes[j] <- private$OUPPassageTimeProbability(s,z[j],tmodes[j],k,omega,rho,mu,sigma) }
            else { Ptmodes[j] <- 0.5*PInf }
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
          }
          if(kindex > 0 & sigma > 0)
          {
            mediansmall <- private$OUPPassageTimePctSearch(s,k-0.01,k,omega,rho,mu,sigma,0.5)
            medianbig <- private$OUPPassageTimePctSearch(s,k+0.01,k,omega,rho,mu,sigma,0.5)
            if(rho > 0)
            {
              meansmall <- private$OUPPassageTimeMeanIntegrate(s,k-0.01,k,omega,rho,mu,sigma,mediansmall)
              meanbig <- private$OUPPassageTimeMeanIntegrate(s,k+0.01,k,omega,rho,mu,sigma,medianbig)
              Ptmeansmall <- private$OUPPassageTimeProbability(s,k-0.01,meansmall,k,omega,rho,mu,sigma)
              Ptmeanbig <- private$OUPPassageTimeProbability(s,k+0.01,meanbig,k,omega,rho,mu,sigma)
              Ptmeans[kindex] <- 0.5*(Ptmeansmall+Ptmeanbig)
            }
            modesmall <- private$OUPPassageTimeModeSearch(s,k-0.01,k,omega,rho,mu,sigma,mediansmall)
            modebig <- private$OUPPassageTimeModeSearch(s,k+0.01,k,omega,rho,mu,sigma,medianbig)
            Ptmodesmall <- private$OUPPassageTimeProbability(s,k-0.01,modesmall,k,omega,rho,mu,sigma)
            Ptmodebig <- private$OUPPassageTimeProbability(s,k+0.01,modebig,k,omega,rho,mu,sigma)
            Ptmodes[kindex] <- 0.5*(Ptmodesmall+Ptmodebig)
          }
          Ptmeshmeans <- MeshWall(z,tmeans,Ptmeans)
          Ptmeshmedians <- MeshWall(z,tmedians,Ptmedians)
          Ptmeshmodes <- MeshWall(z,tmodes,Ptmodes)
          zgap <- 0.03*(z[n]-z[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(z[1]-zgap,z[n]+zgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(0-0.03,1+0.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=grn$e,font=list(color=grn$b))
          Ptmeanline <- list(color=cyn$d,width=8)
          Ptmedianline <- list(color=grn$d,width=8)
          Ptmodeline <- list(color=blu$d,width=8)
          Ptxline <- list(color=gry$d,width=6)
          Ptline <- list(color=grn$d,width=6)
          gradientmean <- list(c(0,cyn$b),c(1,cyn$b))
          gradientmedian <- list(c(0,grn$b),c(1,grn$b))
          gradientmode <- list(c(0,blu$b),c(1,blu$b))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeModeMedianMean3Dprobability")
          fig <- plot_ly()
          if(rho > 0) { fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=rep(0,n),name="mean(<i>z</i>)",mode="lines",line=meanline,hoverinfo="text",text=coordinatemeans) }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=rep(0,n),name="median(<i>z</i>)",mode="lines",line=medianline,hoverinfo="text",text=coordinatemedians) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=rep(0,n),name="mode(<i>z</i>)",mode="lines",line=modeline,hoverinfo="text",text=coordinatemodes)
          if(rho > 0)
          {
            fig <- add_trace(fig,type="scatter3d",x=z,y=tmeans,z=Ptmeans,name="<i>P<sub>t</sub></i>(mean)",mode="lines",line=Ptmeanline,hoverinfo="text",text=coordinatePtmeans,legendgroup="Ptmeans",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmeans$xvertex,y=Ptmeshmeans$yvertex,z=Ptmeshmeans$zvertex,i=Ptmeshmeans$ivertex,j=Ptmeshmeans$jvertex,k=Ptmeshmeans$kvertex,intensity=Ptmeshmeans$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmean,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePtmeans,legendgroup="Ptmeans",visible="legendonly",showlegend=FALSE)
          }
          fig <- add_trace(fig,type="scatter3d",x=z,y=tmedians,z=Ptmedians,name="<i>P<sub>t</sub></i>(median)",mode="lines",line=Ptmedianline,hoverinfo="text",text=coordinatePtmedians,legendgroup="Ptmedians",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmedians$xvertex,y=Ptmeshmedians$yvertex,z=Ptmeshmedians$zvertex,i=Ptmeshmedians$ivertex,j=Ptmeshmedians$jvertex,k=Ptmeshmedians$kvertex,intensity=Ptmeshmedians$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmedian,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePtmedians,legendgroup="Ptmedians",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tmodes,z=Ptmodes,name="<i>P<sub>t</sub></i>(mode)",mode="lines",line=Ptmodeline,hoverinfo="text",text=coordinatePtmodes,legendgroup="Ptmodes",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshmodes$xvertex,y=Ptmeshmodes$yvertex,z=Ptmeshmodes$zvertex,i=Ptmeshmodes$ivertex,j=Ptmeshmodes$jvertex,k=Ptmeshmodes$kvertex,intensity=Ptmeshmodes$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientmode,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePtmodes,legendgroup="Ptmodes",visible="legendonly",showlegend=FALSE)
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",visible="legendonly")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",visible="legendonly",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Pt",visible="legendonly",showlegend=FALSE) }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,visible="legendonly")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time variances
    #' @param type  = 3 and 4 for 2D
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @param zbeg    begin value for state axis
    #' @param zend    end value for state axis
    #' @return plot
    PlotPassageTimeVariance = function(type=NULL,xaxis=NULL,yaxis=NULL,title=NULL,zbeg=NULL,zend=NULL)
    {
      # set/get ----
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      n <- length(z)
      tvariance <- private$tvariance
      tvariances <- private$tvariances
      if(is.null(tvariance) | is.null(tvariances))
      {
        variances <- self$PassageTimeVariance(plotit=FALSE)
        tvariance <- variances[[1]]
        tvariances <- variances[[2]]
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        tvariances <- tvariances[Ixbeg:Ixend]
        n <- length(z)
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
          if(rho > 0) { title <- "Passage Time Variance" }
          else { title <- "Passage Time" }
        }
        if(is.null(xaxis)) { xaxis <- paste(sep="","<i>z</i><br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- "<i>z</i><br>" }
      }
      if(is.null(yaxis)) { yaxis <- "(<i>t-s</i>)<sup>2</sup>" }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      minv <- 0
      maxv <- 1
      if(is.finite(tvariances[1]) & is.finite(tvariances[n]))
      {
        minv <- min(tvariances)
        maxv <- max(tvariances)
      }
      maxv <- 1.2*(maxv-minv)
      maxvscale <- 1
      while(maxv > maxvscale) { maxvscale <- 10*maxvscale }
      maxv <- round(maxv/maxvscale,2)*maxvscale+minv
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(z[1],z[n]),zeroline=FALSE)
      if(x > k) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(minv,maxv),side="right") }
      else { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",range=c(minv,maxv)) }
      varianceline <- list(color=mgn$d,width=4,shape="spline",smoothing=1.3)
      fig <- plot_ly()
      # OUP_A_PassageTimeVariance2Dx
      if(type < 3.5)
      {
        variancedashline <- list(color=mgn$d,dash="longdash",width=3)
        variancedotline <- list(color=mgn$d,dash="dot",width=3)
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeVariance2Dx")
        fig <- add_trace(fig,type="scatter",x=c(x,x),y=c(0,tvariance),name="variance",mode="lines",line=variancedotline,hoverinfo="x+y")
        if(x > k)
        {
          fig <- add_trace(fig,type="scatter",x=c(x,z[n]),y=c(tvariance,tvariance),name="variance",mode="lines",line=variancedashline,hoverinfo="x+y")
          xv <- list(x=x,y=tvariance,text=paste(sep="","variance",bsym,"=",esym,format(tvariance,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
        }
        else
        {
          fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tvariance,tvariance),name="variance",mode="lines",line=variancedashline,hoverinfo="x+y")
          xv <- list(x=x,y=tvariance,text=paste(sep="","variance",bsym,"=",esym,format(tvariance,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
        }
        fig <- add_trace(fig,type="scatter",x=z,y=tvariances,name="variance(<i>z</i>)",mode="lines",line=varianceline,hoverinfo="x+y") %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=xv,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # OUP_A_PassageTimeVariance2D
      else
      {
        imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeVariance2D")
        fig <- add_trace(fig,type="scatter",x=z,y=tvariances,name="variance(<i>z</i>)",mode="lines",line=varianceline) %>%
          config(.,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time lower, median and upper percentiles
    #' @param type  = 1, 2, 3 and 4 for 2D, 5 and 6 for 3D
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
      self$set_plot_info(type,NULL,ptmax)
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
      type <- private$plot_info$plottype$type
      ptmax <- private$plot_info$plottype$ptmax
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
      if(m > 1) { dt <- (t[m]-t[1])/(m-1) }
      else { dt <- 0.05 }
      thalf <- private$tmedian
      thalfs <- private$tmedians
      if(is.null(thalf) | is.null(thalfs))
      {
        halfs <- self$PassageTimeMedian(plotit=FALSE)
        thalf <- halfs[[1]]
        thalfs <- halfs[[2]]
      }
      tpercentile <- private$tpercentile
      tpercentiles <- private$tpercentiles
      if(is.null(tpercentile) | is.null(tpercentiles))
      {
        percentiles <- self$PassageTimePercentiles(plotit=FALSE)
        tpercentile <- percentiles[[1]]
        tpercentiles <- percentiles[[2]]
      }
      tlower <- tpercentile[[1]]
      tupper <- tpercentile[[2]]
      tlowers <- tpercentiles[[1]]
      tuppers <- tpercentiles[[2]]
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) | is.null(pt))
      {
        densities <- self$PassageTimeDensity(plotit=FALSE)
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) | is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(plotit=FALSE)
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        thalfs <- thalfs[Ixbeg:Ixend]
        tlowers <- tlowers[Ixbeg:Ixend]
        tuppers <- tuppers[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend]
        Pt <- Pt[,Ixbeg:Ixend]
        n <- length(z)
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
      if(type < 2.5)
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
        upperdotline <- list(color=grn$d,dash="dot",width=3)
        halfdashline <- list(color=grn$d,dash="dash",width=3)
        halfdotline <- list(color=grn$d,dash="dot",width=3)
        lowerdashline <- list(color=grn$d,dash="dot",width=3)
        lowerdotline <- list(color=grn$d,dash="dot",width=3)
        fig <- plot_ly()
        # OUP_A_PassageTimePercentiles2Dpt
        if(type < 1.5)
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
          if(m > 1) { dt = t[2]-t[1] }
          else { dt = 0.1 }
          ptupper <- private$OUPPassageTimeDensity(s,x,tupper,k,omega,rho,mu,sigma,dt)
          pthalf <- private$OUPPassageTimeDensity(s,x,thalf,k,omega,rho,mu,sigma,dt)
          ptlower <- private$OUPPassageTimeDensity(s,x,tlower,k,omega,rho,mu,sigma,dt)
          fig <- add_trace(fig,type="scatter",x=c(tupper,tupper),y=c(0,ptupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tupper),y=c(ptupper,ptupper),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(thalf,thalf),y=c(0,pthalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,thalf),y=c(pthalf,pthalf),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tlower,tlower),y=c(0,ptlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tlower),y=c(ptlower,ptlower),name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
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
          Ptupper <- private$OUPPassageTimeProbability(s,x,tupper,k,omega,rho,mu,sigma)
          Pthalf <- private$OUPPassageTimeProbability(s,x,thalf,k,omega,rho,mu,sigma)
          Ptlower <- private$OUPPassageTimeProbability(s,x,tlower,k,omega,rho,mu,sigma)
          fig <- add_trace(fig,type="scatter",x=c(tupper,tupper),y=c(0,Ptupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,legendgroup="upper",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tupper),y=c(Ptupper,Ptupper),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=upperdotline,legendgroup="upper",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(thalf,thalf),y=c(0,Pthalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,legendgroup="half",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,thalf),y=c(Pthalf,Pthalf),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=halfdotline,legendgroup="half",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(tlower,tlower),y=c(0,Ptlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,legendgroup="lower",hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=c(s,tlower),y=c(Ptlower,Ptlower),name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=lowerdotline,legendgroup="lower",showlegend=FALSE,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=t,y=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=probabilityline,hoverinfo="x+y")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
      }
      # 2D continued
      else if(type < 4.5)
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
        # OUP_A_PassageTimePercentiles2Dx
        if(type < 3.5)
        {
          upperdashline <- list(color=grn$d,dash="longdash",width=3)
          upperdotline <- list(color=grn$d,dash="dot",width=3)
          halfdashline <- list(color=grn$d,dash="dash",width=3)
          halfdotline <- list(color=grn$d,dash="dot",width=3)
          lowerdashline <- list(color=grn$d,dash="dot",width=3)
          lowerdotline <- list(color=grn$d,dash="dot",width=3)
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
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>",bsym,"=",esym,format(tupper,digits=4),"<br><i>t</i><sub>",format(Phalf,digits=4),"</sub>",bsym,"=",esym,format(thalf,digits=4),"<br><i>t</i><sub>",format(Plower,digits=4),"</sub>",bsym,"=",esym,format(tlower,digits=4)),xref="x",yref="y",xanchor="right",yanchor="bottom",align="right",showarrow=FALSE)
          }
          else
          {
            fig <- add_trace(fig,type="scatter",x=c(z[1],x),y=c(tupper,tupper),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(thalf,thalf),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfdashline,hoverinfo="x+y") %>%
              add_trace(.,type="scatter",x=c(z[1],x),y=c(tlower,tlower),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerdashline,hoverinfo="x+y")
            xmmm <- list(x=x,y=legendy,text=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>",bsym,"=",esym,format(tupper,digits=4),"<br><i>t</i><sub>",format(Phalf,digits=4),"</sub>",bsym,"=",esym,format(thalf,digits=4),"<br><i>t</i><sub>",format(Plower,digits=4),"</sub>",bsym,"=",esym,format(tlower,digits=4)),xref="x",yref="y",xanchor="left",yanchor="bottom",align="left",showarrow=FALSE)
          }
          fig <- add_trace(fig,type="scatter",x=z,y=tuppers,name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>"),mode="lines",line=upperline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=thalfs,name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>"),mode="lines",line=halfline,hoverinfo="x+y") %>%
            add_trace(.,type="scatter",x=z,y=tlowers,name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>"),mode="lines",line=lowerline,hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=xmmm,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimePercentiles2D
        else
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
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
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
        # OUP_A_PassageTimePercentiles3Ddensity
        if(type < 5.5)
        {
          if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          if(k > mu) { spy <- list(x=2.35,y=0.85,z=0.5) }
          else if(k == mu) { spy <- list(x=0,y=2.55,z=0.5) }
          else { spy <- list(x=-2.35,y=0.85,z=0.5) }
          mindensity <- 0
          ptuppers <- vector("double",n)
          pthalfs <- vector("double",n)
          ptlowers <- vector("double",n)
          coordinateuppers <- vector("double",n)
          coordinatehalfs <- vector("double",n)
          coordinatelowers <- vector("double",n)
          coordinateptuppers <- vector("double",n)
          coordinatepthalfs <- vector("double",n)
          coordinateptlowers <- vector("double",n)
          coordinatex <- vector("double",m)
          coordinates <- matrix("",m,n)
          for(j in 1:n)
          {
            ptuppers[j] <- private$OUPPassageTimeDensity(s,z[j],tuppers[j],k,omega,rho,mu,sigma,dt)
            pthalfs[j] <- private$OUPPassageTimeDensity(s,z[j],thalfs[j],k,omega,rho,mu,sigma,dt)
            ptlowers[j] <- private$OUPPassageTimeDensity(s,z[j],tlowers[j],k,omega,rho,mu,sigma,dt)
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
          }
          ptmeshuppers <- MeshWall(z,tuppers,ptuppers)
          ptmeshhalfs <- MeshWall(z,thalfs,pthalfs)
          ptmeshlowers <- MeshWall(z,tlowers,ptlowers)
          zgap <- 0.03*(z[n]-z[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(z[1]-zgap,z[n]+zgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
          else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(mindensity-0.03*(ptmax-mindensity),ptmax),tickmode="auto",nticks=5,mirror=TRUE) }
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=blu$e,font=list(color=blu$b))
          upperline <- list(color=blu$e,width=8)
          halfline <- list(color=blu$e,width=8)
          lowerline <- list(color=blu$e,width=8)
          ptupperline <- list(color=blu$d,width=8)
          pthalfline <- list(color=blu$d,width=8)
          ptlowerline <- list(color=blu$d,width=8)
          ptxline <- list(color=gry$d,width=6)
          ptline <- list(color=blu$d,width=6)
          gradientupper <- list(c(0,blu$c),c(1,blu$c))
          gradienthalf <- list(c(0,blu$b),c(1,blu$b))
          gradientlower <- list(c(0,blu$c),c(1,blu$c))
          rise <- list(x=0,y=-900,z=0)
          shine <- list(ambient=0.9,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimePercentiles3Ddensity")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=upperline,hoverinfo="text",text=coordinateuppers) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=halfline,hoverinfo="text",text=coordinatehalfs) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=lowerline,hoverinfo="text",text=coordinatelowers) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=ptuppers,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)"),mode="lines",line=ptupperline,hoverinfo="text",text=coordinateptuppers,legendgroup="ptuppers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshuppers$xvertex,y=ptmeshuppers$yvertex,z=ptmeshuppers$zvertex,i=ptmeshuppers$ivertex,j=ptmeshuppers$jvertex,k=ptmeshuppers$kvertex,intensity=ptmeshuppers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientupper,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateptuppers,legendgroup="ptuppers",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=pthalfs,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)"),mode="lines",line=pthalfline,hoverinfo="text",text=coordinatepthalfs,legendgroup="pthalfs",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshhalfs$xvertex,y=ptmeshhalfs$yvertex,z=ptmeshhalfs$zvertex,i=ptmeshhalfs$ivertex,j=ptmeshhalfs$jvertex,k=ptmeshhalfs$kvertex,intensity=ptmeshhalfs$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthalf,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatepthalfs,legendgroup="pthalfs",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=ptlowers,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)"),mode="lines",line=ptlowerline,hoverinfo="text",text=coordinateptlowers,legendgroup="ptlowers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=ptmeshlowers$xvertex,y=ptmeshlowers$yvertex,z=ptmeshlowers$zvertex,i=ptmeshlowers$ivertex,j=ptmeshlowers$jvertex,k=ptmeshlowers$kvertex,intensity=ptmeshlowers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientlower,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinateptlowers,legendgroup="ptlowers",visible="legendonly",showlegend=FALSE)
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="pt",visible="legendonly")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,hoverinfo="text",text=coordinates[,j],legendgroup="pt",visible="legendonly",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,visible="legendonly")
        }
        # OUP_A_PassageTimePercentiles3Dprobability
        else
        {
          if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
          if(k > mu) { spy <- list(x=0.85,y=-2.35,z=0.5) }
          else if(k == mu) { spy <- list(x=0,y=2.55,z=0.5) }
          else { spy <- list(x=-0.85,y=-2.35,z=0.5) }
          Ptuppers <- vector("double",n)
          Pthalfs <- vector("double",n)
          Ptlowers <- vector("double",n)
          coordinateuppers <- vector("double",n)
          coordinatehalfs <- vector("double",n)
          coordinatelowers <- vector("double",n)
          coordinatePtuppers <- vector("double",n)
          coordinatePthalfs <- vector("double",n)
          coordinatePtlowers <- vector("double",n)
          coordinatex <- vector("double",m)
          coordinates <- matrix("",m,n)
          kindex <- 0
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
          }
          Ptmeshuppers <- MeshWall(z,tuppers,Ptuppers)
          Ptmeshhalfs <- MeshWall(z,thalfs,Pthalfs)
          Ptmeshlowers <- MeshWall(z,tlowers,Ptlowers)
          zgap <- 0.03*(z[n]-z[1])
          tgap <- 0.03*(t[m]-t[1])
          xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(z[1]-zgap,z[n]+zgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,range=c(t[1]-tgap,t[m]+tgap),tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
          zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,range=c(0-0.03,1+0.03),tickmode="auto",nticks=5,mirror=TRUE)
          view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
          hover <- list(bgcolor=grn$e,font=list(color=grn$b))
          upperline <- list(color=grn$e,width=8)
          halfline <- list(color=grn$e,width=8)
          lowerline <- list(color=grn$e,width=8)
          Ptupperline <- list(color=grn$d,width=8)
          Pthalfline <- list(color=grn$d,width=8)
          Ptlowerline <- list(color=grn$d,width=8)
          Ptxline <- list(color=gry$d,width=6)
          Ptline <- list(color=grn$d,width=6)
          gradientupper <- list(c(0,grn$c),c(1,grn$c))
          gradienthalf <- list(c(0,grn$b),c(1,grn$b))
          gradientlower <- list(c(0,grn$c),c(1,grn$c))
          rise <- list(x=0,y=-800,z=0)
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimePercentiles3Dprobability")
          fig <- plot_ly() %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Pupper,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=upperline,hoverinfo="text",text=coordinateuppers) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Phalf,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=halfline,hoverinfo="text",text=coordinatehalfs) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=rep(0,n),name=paste(sep="","<i>t</i><sub>",format(Plower,digits=4),"</sub>(<i>z</i>)"),mode="lines",line=lowerline,hoverinfo="text",text=coordinatelowers) %>%
            add_trace(.,type="scatter3d",x=z,y=tuppers,z=Ptuppers,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Pupper,digits=4),"</sub>)"),mode="lines",line=Ptupperline,hoverinfo="text",text=coordinatePtuppers,legendgroup="Ptuppers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshuppers$xvertex,y=Ptmeshuppers$yvertex,z=Ptmeshuppers$zvertex,i=Ptmeshuppers$ivertex,j=Ptmeshuppers$jvertex,k=Ptmeshuppers$kvertex,intensity=Ptmeshuppers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientupper,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePtuppers,legendgroup="Ptuppers",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=thalfs,z=Pthalfs,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Phalf,digits=4),"</sub>)"),mode="lines",line=Pthalfline,hoverinfo="text",text=coordinatePthalfs,legendgroup="Pthalfs",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshhalfs$xvertex,y=Ptmeshhalfs$yvertex,z=Ptmeshhalfs$zvertex,i=Ptmeshhalfs$ivertex,j=Ptmeshhalfs$jvertex,k=Ptmeshhalfs$kvertex,intensity=Ptmeshhalfs$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradienthalf,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePthalfs,legendgroup="Pthalfs",visible="legendonly",showlegend=FALSE) %>%
            add_trace(.,type="scatter3d",x=z,y=tlowers,z=Ptlowers,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i><sub>",format(Plower,digits=4),"</sub>)"),mode="lines",line=Ptlowerline,hoverinfo="text",text=coordinatePtlowers,legendgroup="Ptlowers",visible="legendonly") %>%
            add_trace(.,type="mesh3d",x=Ptmeshlowers$xvertex,y=Ptmeshlowers$yvertex,z=Ptmeshlowers$zvertex,i=Ptmeshlowers$ivertex,j=Ptmeshlowers$jvertex,k=Ptmeshlowers$kvertex,intensity=Ptmeshlowers$zvertex,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradientlower,reversescale=reverse,opacity=0.5,hoverinfo="text",text=coordinatePtlowers,legendgroup="Ptlowers",visible="legendonly",showlegend=FALSE)
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(t|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",visible="legendonly")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",visible="legendonly",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Pt",visible="legendonly",showlegend=FALSE) }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(t|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,visible="legendonly")
        }
        fig <- config(fig,toImageButtonOptions=imageoptions) %>%
          layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
      }
      return(fig)
    },
    #' @description
    #' Plot passage time densities
    #' @param type  = 1 and 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type,NULL,ptmax)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      type <- private$plot_info$plottype$type
      ptmax <- private$plot_info$plottype$ptmax
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      blu <- private$plot_colors$blu
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(t)
      n <- length(z)
      ptx <- private$ptx
      pt <- private$pt
      if(is.null(ptx) | is.null(pt))
      {
        densities <- self$PassageTimeDensity(plotit=FALSE)
        ptx <- densities[[1]]
        pt <- densities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        ptx <- ptx[Ixbeg:Ixend]
        pt <- pt[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        pt <- pt[,Ixbeg:Ixend]
        n <- length(z)
      }
      # plot ----
      # 2D
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
      if(type < 2.5)
      {
        if(labels == TRUE)
        {
          if(is.null(xaxis)) { xaxis <- paste(sep="","<i>t</i><br>",syms) }
          else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
        }
        else if(is.null(title)) { title <- "" }
        if(is.null(yaxis)) { yaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        lookdown <- list(text=xaxis)
        lookleft <- list(text=yaxis)
        horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
        if(is.nan(ptmax)) { vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,rangemode="tozero",ticks="outside") }
        else
        {
          mindensity <- 9
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
        fig <- plot_ly()
        # OUP_A_PassageTimeDensity2DManyCurves
        if(type < 1.5)
        {
          legendpos <- list(orientation="h",x=1.0,y=1.0,xanchor="right")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeDensity2DManyCurves")
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          fig <- add_trace(fig,type="scatter",x=t,y=pt[,j],name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)"),mode="lines",line=ptline,opacity=lineopacity,legendgroup="pt",hoverinfo="x+y")
          while(j <= n)
          {
            fig <- add_trace(fig,type="scatter",x=t,y=pt[,j],name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|",z[j],")"),mode="lines",line=ptline,opacity=lineopacity,legendgroup="pt",showlegend=FALSE,hoverinfo="x+y")
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
          }
          fig <- add_trace(fig,type="scatter",x=t,y=ptx,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=ptxline,visible="legendonly",hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimeDensity2DOneCurve
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeDensity2DOneCurve")
          fig <- add_trace(fig,type="scatter",x=t,y=ptx,name=paste(sep="","<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=ptxline) %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>p<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        mindensity <- 9
        coordinatex <- vector("double",m)
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
        if( mindensity > 0) { mindensity <- 0 }
        hover <- list(bgcolor=blu$e,font=list(color=blu$b))
        if(k > mu) { spy <- list(x=2.35,y=0.85,z=0.5) }
        else if(k == mu) { spy <- list(x=0,y=2.35,z=0.5) }
        else { spy <- list(x=-2.35,y=0.85,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        if(is.nan(ptmax)) { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE) }
        else { zview <- list(title=zaxis,color=font$color,linecolor=blu$c,linewidth=3,gridcolor=blu$c,gridwidth=2,backgroundcolor=blu$b,showbackground=floor,range=c(mindensity-0.03*(ptmax-mindensity),ptmax),tickmode="auto",nticks=5,mirror=TRUE) }
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_PassageTimeDensity3DSurface
        if(type < 3.5)
        {
          gradient <- list(c(0,blu$c),c(1,blu$c))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeDensity3DSurface")
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=z,y=t,z=pt,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
        # OUP_A_PassageTimeDensity3DSurfaceScatter
        else if(type < 4.5)
        {
          ptxline <- list(color=font$color,width=8)
          ptline <- list(color=blu$d,width=8)
          gradient <- list(c(0,blu$c),c(1,blu$c))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeDensity3DSurfaceScatter")
          fig <- plot_ly()
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          j <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=ptline,hoverinfo="text",text=coordinates[,j],legendgroup="pt")
          j <- j+dx
          while(j <= n)
          {
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,hoverinfo="text",text=coordinates[,j],legendgroup="pt",showlegend=FALSE)
            j <- j+dx
          }
          fig <- add_trace(fig,type="surface",x=z,y=t,z=pt,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,legendgroup="pt",showlegend=FALSE)
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,visible="legendonly") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
        # OUP_A_PassageTimeDensity3DScatter
        else
        {
          ptxline <- list(color=gry$d,width=6)
          ptline <- list(color=blu$d,width=6)
          legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeDensity3DScatter")
          fig <- plot_ly()
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],name="<i>p<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="pt")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=pt[,j],mode="lines",line=ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="pt",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=ptx,name="<i>p<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=ptxline,hoverinfo="text",text=coordinatex,visible="legendonly") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
      }
      return(fig)
    },
    #' @description
    #' Plot passage time probabilities
    #' @param type  = 1 and 2 for 2D, 3, 4 and 5 for 3D
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
      self$set_plot_info(type)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      t <- private$t_stoch_args[[1]]
      k <- private$t_stoch_args[[2]]
      s <- private$t_stoch_args[[3]]
      x <- private$t_stoch_args[[4]]
      z <- private$t_stoch_args[[5]]
      omega <- private$t_stoch_args[[6]]
      type <- private$plot_info$plottype$type
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      walls <- private$plot_info$plot3D$walls
      floor <- private$plot_info$plot3D$floor
      labels <- private$plot_info$plotlabels
      grn <- private$plot_colors$grn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      reverse <- private$plot_colors$reverse
      m <- length(t)
      n <- length(z)
      Ptx <- private$Ptx
      Pt <- private$Pt
      if(is.null(Ptx) | is.null(Pt))
      {
        probabilities <- self$PassageTimeProbability(plotit=FALSE)
        Ptx <- probabilities[[1]]
        Pt <- probabilities[[2]]
      }
      Inx <- index(t,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < m)
      {
        t <- t[Ixbeg:Ixend]
        Ptx <- Ptx[Ixbeg:Ixend]
        Pt <- Pt[Ixbeg:Ixend,]
        m <- length(t)
      }
      Inx <- index(z,zbeg,zend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 | Ixend < n)
      {
        z <- z[Ixbeg:Ixend]
        Pt <- Pt[,Ixbeg:Ixend]
        n <- length(z)
      }
      # plot ----
      # 2D
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
      if(type < 2.5)
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
        if(Pt[m,1] > 0.8 | Pt[m,n] > 0.8)
        {
          if((Pt[m,1] < 0.5 | Pt[m,n] < 0.5) & sigma/(2*rho)^0.5 < 15) { legendy <- 0.9 }
          else { legendy <- 0.15 }
        }
        fig <- plot_ly()
        # OUP_A_PassageTimeProbability2DManyCurves
        if(type < 1.5)
        {
          kindex <- 0
          for(j in 1:n)  { if(z[j] == k) { kindex <- j } }
          legendpos <- list(orientation="h",x=1.0,y=legendy,xanchor="right")
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeProbability2DManyCurves")
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          fig <- add_trace(fig,type="scatter",x=t,y=Pt[,j],name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)"),mode="lines",line=Ptline,opacity=lineopacity,legendgroup="Pt",hoverinfo="x+y")
          while(j <= n)
          {
            fig <- add_trace(fig,type="scatter",x=t,y=Pt[,j],name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|",z[j],")"),mode="lines",line=Ptline,opacity=lineopacity,legendgroup="Pt",showlegend=FALSE,hoverinfo="x+y")
            j <- j+dx
            q <- q+1
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter",x=c(t[1],t[1]),y=c(0,Pt[1,kindex]),mode="lines",line=Ptline,legendgroup="Pt",showlegend=FALSE,hoverinfo="x+y") }
          fig <- add_trace(fig,type="scatter",x=t,y=Ptx,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=Ptxline,visible="legendonly",hoverinfo="x+y") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
        # OUP_A_PassageTimeProbability2DOneCurve
        else
        {
          imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_A_PassageTimeProbability2DOneCurve")
          fig <- add_trace(fig,type="scatter",x=t,y=Ptx,name=paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)"),mode="lines",line=Ptxline) %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))
        }
      }
      # 3D
      else
      {
        if(labels == TRUE) { lookdown <- list(text=syms,showarrow=FALSE,yref="container",y=0) }
        else { lookdown <- list(text="",showarrow=FALSE,yref="container",y=0) }
        if(is.null(xaxis)) { xaxis <- "<i>z</i>" }
        if(is.null(yaxis)) { yaxis <- "<i>t</i>" }
        if(is.null(zaxis)) { zaxis <- "<i>P<sub>t</sub></i>(<i>t</i>|<i>k,s,x</i>)" }
        kindex <- 0
        coordinatex <- vector("double",m)
        coordinates <- matrix("",m,n)
        for(i in 1:m)
        {
          coordinatex[i] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)=",format(Ptx[i],digits=4),"<br><i>t</i>=",t[i],"<br><i>x</i>=",x)
          for(j in 1:n) { coordinates[i,j] <- paste(sep="","<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)=",format(Pt[i,j],digits=4),"<br><i>t</i>=",t[i],"<br><i>z</i>=",z[j]) }
        }
        for(j in 1:n) { if(z[j] == k) { kindex <- j } }
        hover <- list(bgcolor=grn$e,font=list(color=grn$b))
        if(k > mu) { spy <- list(x=2.35,y=0.85,z=0.5) }
        else if(k == mu) { spy <- list(x=0,y=2.35,z=0.5) }
        else { spy <- list(x=-2.3,y=0.85,z=0.5) }
        xview <- list(title=xaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        yview <- list(title=yaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$a,showbackground=walls,tickmode="auto",nticks=5,zeroline=FALSE,mirror=TRUE)
        zview <- list(title=zaxis,color=font$color,linecolor=grn$c,linewidth=3,gridcolor=grn$c,gridwidth=2,backgroundcolor=grn$b,showbackground=floor,rangemode="tozero",tickmode="auto",nticks=5,mirror=TRUE)
        view <- list(camera=list(eye=spy),xaxis=xview,yaxis=yview,zaxis=zview,aspectratio=list(x=1,y=1,z=1))
        # OUP_A_PassageTimeProbability3DSurface
        if(type < 3.5)
        {
          gradient <- list(c(0,grn$c),c(1,grn$c))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeProbability3DSurface")
          fig <- plot_ly() %>%
            add_trace(.,type="surface",x=z,y=t,z=Pt,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates) %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,showlegend=FALSE,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
        # OUP_A_PassageTimeProbability3DSurfaceScatter
        else if(type < 4.5)
        {
          Ptxline <- list(color=gry$d,width=8)
          Ptline <- list(color=grn$e,width=8)
          gradient <- list(c(0,grn$c),c(1,grn$c))
          if(k > mu) { rise <- list(x=10,y=100,z=0) }
          else if(k == mu) { rise <- list(x=0,y=100,z=0) }
          else { rise <- list(x=-10,y=100,z=0) }
          shine <- list(ambient=0.7,diffuse=0.5,fresnel=0.2,roughness=0.5,specular=0.1)
          legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeProbability3DSurfaceScatter")
          fig <- plot_ly()
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          j <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=Ptline,hoverinfo="text",text=coordinates[,j],legendgroup="Pt")
          j <- j+dx
          while(j <= n)
          {
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",showlegend=FALSE)
            j <- j+dx
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Pt",showlegend=FALSE) }
          fig <-add_trace(fig,type="surface",x=z,y=t,z=Pt,showscale=FALSE,lighting=shine,lightposition=rise,colorscale=gradient,reversescale=reverse,hoverinfo="text",text=coordinates,legendgroup="Pt",showlegend=FALSE)
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,visible="legendonly") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
        # OUP_A_PassageTimeProbability3DScatter
        else
        {
          Ptxline <- list(color=gry$d,width=6)
          Ptline <- list(color=grn$d,width=6)
          legendpos <- list(x=1.0,y=0.5,xanchor="right",yanchor="center",tracegroupgap=0,itemsizing="constant")
          imageoptions <- list(format=file$format,width=file$width,height=file$width,filename="OUP_A_PassageTimeProbability3DScatter")
          fig <- plot_ly()
          xx <- vector("double",m)
          kk <- vector("double",m)
          dx <- as.integer((n-1)/10)
          if(dx < 1) { dx <- 1 }
          lineopacity <- 0.75
          j <- 1
          q <- 1
          for(i in 1:m) { xx[i] <- z[j] }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],name="<i>P<sub>t</sub></i>(<i>t</i>|<i>z</i>)",mode="lines",line=Ptline,opacity=lineopacity,hoverinfo="text",text=coordinates[,j],legendgroup="Pt")
          j <- j+dx
          q <- q+1
          while(j <= n)
          {
            if(q < 7) { lineopacity <- lineopacity+0.05 }
            else { lineopacity <- lineopacity-0.05 }
            for(i in 1:m) { xx[i] <- z[j] }
            fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Pt[,j],mode="lines",line=Ptline,hoverinfo="text",text=coordinates[,j],legendgroup="Pt",showlegend=FALSE)
            j <- j+dx
            q <- q+1
          }
          if(kindex > 0) { fig <- add_trace(fig,type="scatter3d",x=c(k,k),y=c(t[1],t[1]),z=c(0,Pt[1,kindex]),mode="lines",line=Ptline,hoverinfo="text",text=coordinates[1,kindex],legendgroup="Pt",showlegend=FALSE) }
          for(i in 1:m)
          {
            xx[i] <- x
            kk[i] <- k
          }
          fig <- add_trace(fig,type="scatter3d",x=xx,y=t,z=Ptx,name="<i>P<sub>t</sub></i>(<i>t</i>|<i>x</i>)",mode="lines",line=Ptxline,hoverinfo="text",text=coordinatex,visible="legendonly") %>%
            config(.,toImageButtonOptions=imageoptions) %>%
            layout(.,title=lookup,annotations=lookdown,scene=view,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,margin=list(t=0,r=0,b=0,l=0))
        }
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
    undo_args = NULL,
    plot_info = NULL,
    # private globals ----
    psiphi = NULL,
    undoIx = NULL,
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
    Oneg = NULL,
    Opos = NULL,
    Oshatneg = NULL,
    Oshatpos = NULL,
    Oscarfneg = NULL,
    Oscarfpos = NULL,
    sscarfneg = NULL,
    sscarfpos = NULL,
    KOneg = NULL,
    KOpos = NULL,
    BCneg = NULL,
    BCpos = NULL,
    tmode = NULL,
    tmodes = NULL,
    tmedian = NULL,
    tmedians = NULL,
    tmean = NULL,
    tmeans = NULL,
    tvariance = NULL,
    tvariances = NULL,
    tpercentile = NULL,
    tpercentiles = NULL,
    ptx = NULL,
    pt = NULL,
    Ptx = NULL,
    Pt = NULL,
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
    lists_equal = function(list1,list2)
    {
      allequal <- TRUE
      n1 <- length(list1)
      n2 <- length(list2)
      if(n1 == n2)
      {
        j <- 0
        while(j < n1 & allequal == TRUE )
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
    OUPDensity = function(s,x,t,y,rho,mu,sigma,dy)
    {
      mean <- mu+(x-mu)*exp(-rho*(t-s))
      if(rho < 0.0000000001) { variance <- sigma^2*(t-s) }
      else { variance <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
      if(variance^0.5 < dy)
      {
        lomean <- mean-0.5*dy
        himean <- mean+0.5*dy
        if(y >= lomean & y < himean)
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
    OUPDoubleIntegral = function(s,x,t,y,rho,mu,sigma,psi)
    {
      ltgt <- 1
      if(psi > 0) { ltgt <- -1 }
      if(t <= s) { doubleintegral <- ltgt*(y-x)*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,psi) }
      else if(sigma == 0) { doubleintegral <- ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,psi) }
      else if(rho<0.0000000001) { doubleintegral <- ltgt*(y-x)*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,psi)+sigma^2*(t-s)*private$OUPDnsty(s,x,t,y,rho,mu,sigma) }
      else { doubleintegral <- ltgt*(y-mu-(x-mu)*exp(-rho*(t-s)))*private$OUPPrbblty(s,x,t,y,rho,mu,sigma,psi)+sigma^2/(2*rho)*(1-exp(-2*rho*(t-s)))*private$OUPDnsty(s,x,t,y,rho,mu,sigma) }
      return(doubleintegral)
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
    OUPOptionMax = function(x,y,rho,mu,sigma,r,phi,b,c,tsguess,tsmax,tsmin,maxmin)
    {
      if(rho < 0.0000000001 & r < 0.0000000001)
      {
        if(maxmin > 0) { return(c(Inf,0)) }
        else
        {
          if(phi > 0)
          {
            if(x > y) { return(c(x-y+b,0)) }
            else { return(c(b,0)) }
          }
          else
          {
            if(y > x) { return(c(y-x-c,0)) }
            else { return(c(-c,0))  }
          }
        }
      }
      else
      {
        sigdig <- 12
        n <- 1000
        eps <- 0.000001
        Optn <- vector("double",5)
        ts <- vector("double",5)
        if(tsguess > 0)
        {
          ts[3] <- tsguess
          Optn[3] <- private$OUPOption(-ts[3],x,0,y,rho,mu,sigma,r,phi,b,c)
        }
        else
        {
          dt <- tsmax/11
          ts[3] <- tsmax
          ts[2] <- ts[3]-dt
          Optn[3] <- private$OUPOption(-ts[3],x,0,y,rho,mu,sigma,r,phi,b,c)
          Optn[2] <- private$OUPOption(-ts[2],x,0,y,rho,mu,sigma,r,phi,b,c)
          while(Optn[2] >= Optn[3] & ts[2] > tsmin & dt > eps)
          {
            i <- 10
            while(Optn[2] >= Optn[3] & ts[2] > tsmin & i > 1)
            {
              i <- i-1
              ts[3] <- ts[2]
              Optn[3] <- Optn[2]
              ts[2] <- ts[2]-dt
              Optn[2] <- private$OUPOption(-ts[2],x,0,y,rho,mu,sigma,r,phi,b,c)
            }
            dt <- dt/10
          }
        }
        dt <- (ts[3]-tsmin)/2
        ts[1] <- tsmin
        ts[2] <- ts[1]+dt
        ts[4] <- ts[3]+dt
        ts[5] <- ts[4]+dt
        Optn[1] <- private$OUPOption(-ts[1],x,0,y,rho,mu,sigma,r,phi,b,c)
        Optn[2] <- private$OUPOption(-ts[2],x,0,y,rho,mu,sigma,r,phi,b,c)
        Optn[4] <- private$OUPOption(-ts[4],x,0,y,rho,mu,sigma,r,phi,b,c)
        Optn[5] <- private$OUPOption(-ts[5],x,0,y,rho,mu,sigma,r,phi,b,c)
        m <- 1
        i <- 1
        while(i < 5)
        {
          i <- i+1
          if(maxmin*(Optn[i]-Optn[m]) > 0) { m <- i }
        }
        j <- 0
        while(10^sigdig*abs(Optn[4]-Optn[2]) >= abs(Optn[3]) & ts[3] < tsmax & ts[3] > tsmin & j < n)
        {
          j <- j+1
          if(m == 2 | m == 3 | m == 4)
          {
            if(m == 2)
            {
              ts[5] <- ts[3]
              Optn[5] <- Optn[3]
              ts[3] <- ts[2]
              Optn[3] <- Optn[2]
            }
            else if(m == 3)
            {
              ts[1] <- ts[2]
              Optn[1] <- Optn[2]
              ts[5] <- ts[4]
              Optn[5] <- Optn[4]
            }
            else
            {
              ts[1] <- ts[3]
              Optn[1] <- Optn[3]
              ts[3] <- ts[4]
              Optn[3] <- Optn[4]
            }
            ts[2] <- (ts[3]+ts[1]) / 2
            Optn[2] <- private$OUPOption(-ts[2],x,0,y,rho,mu,sigma,r,phi,b,c)
            ts[4] <- (ts[3]+ts[5]) / 2
            Optn[4] <- private$OUPOption(-ts[4],x,0,y,rho,mu,sigma,r,phi,b,c)
            if(maxmin*(Optn[2]-Optn[3]) > 0) { m <- 2 }
            else if(maxmin*(Optn[4]-Optn[3]) > 0) { m <- 4 }
            else { m <- 3 }
          }
          else if(m == 5)
          {
            i <- 1
            while(i < 5)
            {
              i <- i+1
              ts[i-1] <- ts[i]
              Optn[i-1] <- Optn[i]
            }
            ts[5] <- ts[5]+dt
            Optn[5] <- private$OUPOption(-ts[5],x,0,y,rho,mu,sigma,r,phi,b,c)
            if(maxmin*(Optn[5]-Optn[4]) > 0)
            {
              m <- 5
              dt <- 2*dt
            }
            else { m <- 4 }
          }
          else
          {
            i <- 6
            while(i > 2)
            {
              i <- i-1
              ts[i] <- ts[i-1]
              Optn[i] <- Optn[i-1]
            }
            ts[1] <- ts[1]-dt
            if(ts[1] < 0) { ts[1] <- 0 }
            Optn[1] <- private$OUPOption(-ts[1],x,0,y,rho,mu,sigma,r,phi,b,c)
            if(maxmin*(Optn[1]-Optn[2]) > 0)
            {
              m <- 1
              dt <- 2*dt
            }
            else { m <- 2 }
          }
        }
        if(maxmin*(Optn[3]-Optn[1]) <= 0) { return(c(Optn[1],ts[1])) }
        else { return(c(Optn[3],ts[3])) }
      }
    },
    OUPOptionDtZero = function()
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
      m <- length(s)
      n <- length(x)
      Oscarf <- matrix(0.0,2,n)
      sscarf <- matrix(0.0,2,n)
      if(phi > 0 )
      {
        j <- 1
        while(x[j] <= y & j < n+1)
        {
          Oscarf[1,j] <- NA
          sscarf[1,j] <- 0
          j <- j+1
        }
        tsguess <- 0.001
        tsmax <- 1.1*(t-s[m])
        tsmin <- 0.0
        while(tsguess < tsmax & j < n+1)
        {
          env <- private$OUPOptionMax(x[j],y,rho,mu,sigma,r,phi,b,c,tsguess[1],tsmax,tsmin,-1)
          if(env[2] > 0) { Oscarf[1,j] <- env[1] }
          else { Oscarf[1,j] <- NA }
          sscarf[1,j] <- env[2]
          tsguess <- env[2]+0.001
          j <- j+1
        }
        if(tsguess >= tsmax & j > 1) { j <- j-1 }
        jj <- j
        while(jj < n+1)
        {
          Oscarf[1,jj] <- NA
          sscarf[1,jj] <- tsmax
          jj <- jj+1
        }
        k <- 1
        tsguess <- 0
        env <- private$OUPOptionMax(x[k],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,tsmin,1)
        if(env[2] > 0) { Oscarf[2,k] <- env[1] }
        else { Oscarf[2,k] <- NA }
        sscarf[2,k] <- env[2]
        tsguess <- env[2]
        tsmin <- sscarf[1,k]
        k <- k+1
        while(tsguess > tsmin & k < n+1)
        {
          env <- private$OUPOptionMax(x[k],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,tsmin,1)
          if(env[2] > 0) { Oscarf[2,k] <- env[1] }
          else { Oscarf[2,k] <- NA }
          sscarf[2,k] <- env[2]
          tsguess <- env[2]
          tsmin <- sscarf[1,k]
          k <- k+1
        }
        if(tsguess <= tsmin & k > 1) { k <- k-1 }
        kk <- k
        while(kk < n+1)
        {
          Oscarf[2,kk] <- NA
          sscarf[2,kk] <- 0
          kk <- kk+1
        }
        if(j < n & j > 0 & k < n & k > 0)
        {
          if(j == k)
          {
            sscarf[1,j] <- (sscarf[1,j-1]+sscarf[2,k-1])/2
            Oscarf[1,j] <- private$OUPOption(-sscarf[1,j],x[j],0,y,rho,mu,sigma,r,phi,b,c)
            sscarf[2,k] <- sscarf[1,j]
            Oscarf[2,k] <- Oscarf[1,j]
          }
          else if(j == k-1)
          {
            sscarf[1,j] <- (sscarf[1,j-1]+sscarf[2,k-1])/2
            Oscarf[1,j] <- Oscarf[2,k-1]
          }
          else if(j == k+1)
          {
            sscarf[2,k] <- (sscarf[1,j-1]+sscarf[2,k-1])/2
            Oscarf[2,k] <- Oscarf[1,j-1]
          }
        }
        private$Oscarfpos <- Oscarf
        private$sscarfpos <- t-sscarf
      }
      else
      {
        j <- n
        while(x[j] >= y & j > 0)
        {
          Oscarf[1,j] <- NA
          sscarf[1,j] <- 0
          j <- j-1
        }
        tsguess <- 0.001
        tsmax <- 1.1*(t-s[m])
        tsmin <- 0.0
        while(tsguess < tsmax & j > 0)
        {
          env <- private$OUPOptionMax(x[j],y,rho,mu,sigma,r,phi,b,c,tsguess,tsmax,tsmin,-1)
          if(env[2] > 0) { Oscarf[1,j] <- env[1] }
          else { Oscarf[1,j] <- NA }
          sscarf[1,j] <- env[2]
          tsguess <- env[2]+0.001
          j <- j-1
        }
        if(tsguess >= tsmax  & j < n) { j <- j+1 }
        jj <- j
        while(jj > 0)
        {
          Oscarf[1,jj] <- NA
          sscarf[1,jj] <- tsmax
          jj <- jj-1
        }
        k <- n
        tsguess <- 0
        env <- private$OUPOptionMax(x[k],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,tsmin,1)
        if(env[2] > 0) { Oscarf[2,k] <- env[1] }
        else { Oscarf[2,k] <- NA }
        sscarf[2,k] <- env[2]
        tsguess <- env[2]
        tsmin <- sscarf[1,k]
        k <- k-1
        while(tsguess > tsmin & k > 0)
        {
          env <- private$OUPOptionMax(x[k],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,tsmin,1)
          if(env[2] > 0) { Oscarf[2,k] <- env[1] }
          else { Oscarf[2,k] <- NA }
          sscarf[2,k] <- env[2]
          tsguess <- env[2]
          tsmin <- sscarf[1,k]
          k <- k-1
        }
        if(tsguess <= tsmin & k < n) { k <- k+1 }
        kk <- k
        while(kk > 0)
        {
          Oscarf[2,kk] <- NA
          sscarf[2,kk] <- 0
          kk <- kk-1
        }
        if(j > 1 & j < n & k > 1 & k < n)
        {
          if(j == k)
          {
            sscarf[1,j] <- (sscarf[1,j+1]+sscarf[2,k+1])/2
            Oscarf[1,j] <- private$OUPOption(-sscarf[1,j],x[j],0,y,rho,mu,sigma,r,phi,b,c)
            sscarf[2,k] <- sscarf[1,j]
            Oscarf[2,k] <- Oscarf[1,j]
          }
          else if(j == k+1)
          {
            sscarf[1,j] <- (sscarf[1,j+1]+sscarf[2,k+1])/2
            Oscarf[1,j] <- Oscarf[2,k+1]
          }
          else if(j == k-1)
          {
            sscarf[2,k] <- (sscarf[1,j+1]+sscarf[2,k+1])/2
            Oscarf[2,k] <- Oscarf[1,j+1]
          }
        }
        private$Oscarfneg <- Oscarf
        private$sscarfneg <- sscarf
      }
      return(Oscarf)
    },
    OUPThresholdSearch = function(y,rho,mu,sigma,r,phi,b,c)
    {
      decision <- vector("double",2)
      if(rho < 0.0000000001 & r < 0.0000000001)
      {
        if(phi > 0) { decision[1] <- Inf}
        else { decision[1] <- -Inf}
        decision[2] <- Inf
      }
      else if(sigma^2 < 0.0000000001 & rho < 0.0000000001)
      {
        decision[1] <- y
        if(phi > 0) { decision[2] <- b }
        else { decision[2] <- -c }
      }
      else
      {
        sigdig <- 12
        n <- 1000
        x <- vector("double",5)
        t <- vector("double",5)
        Optn <- vector("double",5)
        tsguess <- 0
        ltgt <- 1
        if(phi > 0) { ltgt <- -1 }
        if(ltgt < 0)
        {
          x[1] <- y
          env <- private$OUPOptionMax(x[1],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,0,1)
          Optn[1] <- env[1]
          t[1] <- env[2]
          dx <- (Optn[1]-b)/2
          m <- 1
          i <- 1
          while(i < 5)
          {
            i <- i+1
            x[i] <- x[i-1]+dx
            env <- private$OUPOptionMax(x[i],y,rho,mu,sigma,r,phi,b,c,env[2],1100,0,1)
            Optn[i] <- env[1]
            t[i] <- env[2]
            if(Optn[i] > ltgt*(y-x[i])+b) { m <- i }
          }
        }
        else
        {
          x[5] <- y
          env <- private$OUPOptionMax(x[5],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,0,1)
          Optn[5] <- env[1]
          t[5] <- env[2]
          dx <- (Optn[5]+c)/2
          m <- 5
          i <- 5
          while(i > 1)
          {
            i <- i-1
            x[i] <- x[i+1]-dx
            env <- private$OUPOptionMax(x[i],y,rho,mu,sigma,r,phi,b,c,env[2],1100,0,1)
            Optn[i] <- env[1]
            t[i] <- env[2]
            if(Optn[i] > ltgt*(y-x[i])-c) { m <- i }
          }
        }
        j <- 0
        while(10^sigdig*(x[4]-x[2]) >= abs(x[3]) & j < n)
        {
          j <- j+1
          if(m == 2 | m == 3 | m == 4)
          {
            if(m == 2)
            {
              x[5] <- x[3]
              Optn[5] <- Optn[3]
              x[3] <- x[2]
              Optn[3] <- Optn[2]
            }
            else if(m == 3)
            {
              x[1] <- x[2]
              Optn[1] <- Optn[2]
              x[5] <- x[4]
              Optn[5] <- Optn[4]
            }
            else
            {
              x[1] <- x[3]
              Optn[1] <- Optn[3]
              x[3] <- x[4]
              Optn[3] <- Optn[4]
            }
            x[2] <- (x[3]+x[1])/2
            tsguess <- (t[3]+t[1])/2
            env <- private$OUPOptionMax(x[2],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,0,1)
            Optn[2] <- env[1]
            t[2] <- env[2]
            x[4] <- (x[3]+x[5])/2
            tsguess <- (t[3]+t[5])/2
            env <- private$OUPOptionMax(x[4],y,rho,mu,sigma,r,phi,b,c,tsguess,1100,0,1)
            Optn[4] <- env[1]
            t[4] <- env[2]
            if(ltgt < 0)
            {
              if(Optn[4] > ltgt*(y-x[4])+b) { m <- 4 }
              else if(Optn[3] > ltgt*(y-x[3])+b) { m <- 3 }
              else { m <- 2 }
            }
            else
            {
              if(Optn[2] > ltgt*(y-x[2])-c) { m <- 2 }
              else if(Optn[3] > ltgt*(y-x[3])-c) { m <- 3 }
              else { m <- 4 }
            }
          }
          else if(m == 5)
          {
            i <- 1
            while(i < 5)
            {
              i <- i+1
              x[i-1] <- x[i]
              Optn[i-1] <- Optn[i]
            }
            x[5] <- x[5]+dx
            env <- private$OUPOptionMax(x[5],y,rho,mu,sigma,r,phi,b,c,t[5],1100,0,1)
            Optn[5] <- env[1]
            t[5] <- env[2]
            if(ltgt < 0)
            {
              if(Optn[5] > ltgt*(y-x[5])+b)
              {
                m <- 5
                dx <- 2*dx
              }
              else { m <- 4 }
            }
            else
            {
              if(Optn[5] > ltgt*(y-x[5])-c) { m <- 4 }
              else
              {
                m <- 5
                dx <- 2*dx
              }
            }
          }
          else
          {
            i <- 6
            while(i > 2)
            {
              i <- i-1
              x[i] <- x[i-1]
              Optn[i] <- Optn[i-1]
            }
            x[1] <- x[1]-dx
            env <- private$OUPOptionMax(x[1],y,rho,mu,sigma,r,phi,b,c,t[1],1100,0,1)
            Optn[1] <- env[1]
            t[1] <- env[2]
            if(ltgt < 0)
            {
              if(Optn[1] > ltgt*(y-x[1])+b) { m <- 2 }
              else
              {
                m <- 1
                dx <- 2*dx
              }
            }
            else
            {
              if(Optn[1] > ltgt*(y-x[1])-c)
              {
                m <- 1
                dx <- 2*dx
              }
              else { m <- 2 }
            }
          }
          decision[1] <- x[3]
          if(ltgt < 0) { decision[2] <- ltgt*(y-x[3])+b }
          else { decision[2] <- ltgt*(y-x[3])-c }
        }
      }
      return(decision)
    },
    OUPPassageTimeModeSearch = function(s,x,k,omega,rho,mu,sigma,median)
    {
      if(median == Inf | abs(x - k) < 0.0000000001 | sigma^2 < 0.0000000001) { mode <- median }
      else
      {
        m <- vector("integer",2)
        dt <- vector("double",2)
        pass <- matrix(0.0,5,2)
        t <- matrix(0.0,5,2)
        sigdig <- 5
        n <- 1000
        dt[1] <- 0.01*(median-s)
        dt[2] <- 0.5*(median-s)
        t[1,1] <- 0
        t[1,2] <- 0
        pass[1,1] <- 0
        pass[1,2] <- 0
        m[1] <- 1
        m[2] <- 1
        i <- 1
        while(i < 5)
        {
          i <- i+1
          c <- 0
          while(c < 2)
          {
            c <- c+1
            t[i,c] <- t[i-1,c]+dt[c]
            pass[i,c] <- private$OUPPassageTimeDensity(0,x,t[i,c],k,omega,rho,mu,sigma,0.05)
            if(pass[i,c] > pass[m[c],c]) { m[c] <- i }
          }
        }
        c <- 0
        if(pass[m[1],1] < 0.000001) { c <- 1 }
        while(c < 2)
        {
          c <- c+1
          j <- 0
          while(10^sigdig*(t[4,c]-t[2,c]) >= t[3,c] & j < n)
          {
            j <- j+1
            if(m[c] == 2 | m[c] == 3 | m[c] == 4)
            {
              if(m[c] == 2)
              {
                t[5,c] <- t[3,c]
                pass[5,c] <- pass[3,c]
                t[3,c] <- t[2,c]
                pass[3,c] <- pass[2,c]
              }
              else if(m[c] == 3)
              {
                t[1,c] <- t[2,c]
                pass[1,c] <- pass[2,c]
                t[5,c] <- t[4,c]
                pass[5,c] <- pass[4,c]
              }
              else
              {
                t[1,c] <- t[3,c]
                pass[1,c] <- pass[3,c]
                t[3,c] <- t[4,c]
                pass[3,c] <- pass[4,c]
              }
              t[2,c] <- (t[3,c]+t[1,c])/2
              pass[2,c] <- private$OUPPassageTimeDensity(0,x,t[2,c],k,omega,rho,mu,sigma,0.05)
              t[4,c] <- (t[3,c]+t[5,c])/2
              pass[4,c] <- private$OUPPassageTimeDensity(0,x,t[4,c],k,omega,rho,mu,sigma,0.05)
              if(pass[2,c] > pass[3,c]) { m[c] <- 2 }
              else if(pass[4,c] > pass[3,c]) { m[c] <- 4 }
              else { m[c] <- 3 }
            }
            else if(m[c] == 5)
            {
              i <- 1
              while(i < 5)
              {
                i <- i+1
                t[i-1,c] <- t[i,c]
                pass[i-1,c] <- pass[i,c]
              }
              t[5,c] <- t[5,c]+dt[c]
              pass[5,c] <- private$OUPPassageTimeDensity(0,x,t[5,c],k,omega,rho,mu,sigma,0.05)
              if(pass[5,c] > pass[4,c])
              {
                m[c] <- 5
                dt[c] <- 2*dt[c]
              }
              else { m[c] <- 4 }
            }
            else
            {
              i <- 6
              while(i > 2)
              {
                i <- i-1
                t[i,c] <- t[i-1,c]
                pass[i,c] <- pass[i-1,c]
              }
              t[1,c] <- t[1,c]-dt[c]
              pass[1,c] <- private$OUPPassageTimeDensity(0,x,t[1,c],k,omega,rho,mu,sigma,0.05)
              if(pass[1,c] > pass[2,c])
              {
                m[c] <- 1
                dt[c] <- 2*dt[c]
              }
              else { m[c] <- 2 }
            }
          }
        }
        if(pass[3,1] > pass[3,2]) { mode <- s+t[3,1] }
        else { mode <- s+t[3,2] }
      }
      return(mode)
    },
    OUPPassageTimePctSearch = function(s,x,k,omega,rho,mu,sigma,Ppct)
    {
      if(k == Inf | k == -Inf) { median <- Inf }
      else if(abs(x - k) < 0.0000000001) { median <- s }
      else if(sigma^2 < 0.0000000001)
      {
        if((x < k & k < mu | x > k & k > mu) & rho > 0) { median <- s + log((x-mu)/(k-mu))/rho }
        else if(omega > 0) { median <- Inf }
        else { median <- s}
      }
      else
      {
        pass <- vector("double",5)
        t <- vector("double",5)
        sigdig <- 5
        n <- 1000
        pinf <- private$OUPPassageTimeProbabilityInf(x,k,omega,rho,mu,sigma)
        t[1] <- 0.0
        pass[1] <- 0.0
        dt <- 0.1
        i <- 1
        while(i < 5)
        {
          i <- i+1
          t[i] <- dt
          pass[i] <- private$OUPPassageTimeProbability(0,x,t[i],k,omega,rho,mu,sigma)
          dt <- dt*10.0
        }
        m <- 1
        while(m < 5 & pass[m] < Ppct*pinf) { m <- m+1 }
        if(m > 1)
        {
          j <- 0
          while(10^sigdig*(t[4]-t[2]) >= t[3] & j < n)
          {
            j <- j+1
            if(m == 2 | m == 3 | m == 4)
            {
              if(m == 2)
              {
                t[5] <- t[3]
                pass[5] <- pass[3]
                t[3] <- t[2]
                pass[3] <- pass[2]
              }
              else if(m == 3)
              {
                t[1] <- t[2]
                pass[1] <- pass[2]
                t[5] <- t[4]
                pass[5] <- pass[4]
              }
              else
              {
                t[1] <- t[3]
                pass[1] <- pass[3]
                t[3] <- t[4]
                pass[3] <- pass[4]
              }
              t[2] <- (t[3]+t[1]) / 2
              pass[2] <- private$OUPPassageTimeProbability(0,x,t[2],k,omega,rho,mu,sigma)
              t[4] <- (t[3]+t[5]) / 2
              pass[4] <- private$OUPPassageTimeProbability(0,x,t[4],k,omega,rho,mu,sigma)
              if(pass[2] > Ppct*pinf) { m <- 2 }
              else if(pass[4] < Ppct*pinf) { m <- 4 }
              else { m <- 3 }
            }
            else if(m == 5)
            {
              i <- 1
              while(i < 5)
              {
                i <- i+1
                t[i-1] <- t[i]
                pass[i-1] <- pass[i]
              }
              t[5] <- t[5]+dt
              pass[5] <- private$OUPPassageTimeProbability(0,x,t[5],k,omega,rho,mu,sigma)
              if(pass[5] <= Ppct*pinf)
              {
                m <- 5
                dt <- 2*dt
              }
              else
              {
                m <- 4
                dt <- 5
              }
            }
            else
            {
              i <- 6
              while(i > 2)
              {
                i <- i-1
                t[i] <- t[i-1]
                pass[i] <- pass[i-1]
              }
              t[1] <- t[1]-dt
              pass[1] <- private$OUPPassageTimeProbability(0,x,t[1],k,omega,rho,mu,sigma)
              if(pass[1] >= Ppct*pinf)
              {
                m <- 1
                dt <- 2*dt
              }
              else
              {
                m <- 2
                dt <- 5
              }
            }
          }
          median <- s+t[3]
        }
        else { median <- s }
      }
      return(median)
    },
    OUPPassageTimeMeanIntegrate = function(s,x,k,omega,rho,mu,sigma,median)
    {
      if(rho < 0.0000000001) { mean <- Inf }
      else if(median == Inf | abs(x - k) < 0.0000000001 | sigma^2 < 0.0000000001) { mean <- median }
      else
      {
        pinf <- private$OUPPassageTimeProbabilityInf(x,k,omega,rho,mu,sigma)
        dss <- (median-s)/100.0
        ss <- 0
        dnsty <- 0
        while(dnsty < 0.000001*pinf & ss < median-s)
        {
          ss <- ss+dss
          dnsty <- private$OUPPassageTimeDensity(0,x,ss,k,omega,rho,mu,sigma,0.05)
        }
        ss <- ss-dss
        t <- private$Laguerret()
        w <- private$Laguerrew()
        wghtfnct <- vector("double",41)
        density <- vector("double",41)
        tscale <- t[5]
        wscale <- 0
        mean <- 0
        i <- 0
        while(i < 41)
        {
          i <- i+1
          wghtfnct[i] <- exp(-t[i])
          t[i] <- t[i]*(median-s-ss)/tscale+ss
          w[i] <- w[i]*(median-s-ss)/tscale
          density[i] <- private$OUPPassageTimeDensity(0,x,t[i],k,omega,rho,mu,sigma,0.05)
          wscale <- wscale+w[i]*density[i]/wghtfnct[i]
          mean <- mean+w[i]*t[i]*density[i]/wghtfnct[i]
        }
        mean <- s+mean/wscale
      }
      return(mean)
    },
    OUPPassageTimeVarianceIntegrate = function(s,x,k,omega,rho,mu,sigma,median,mean)
    {
      if(rho < 0.000000093321 | median == Inf) { variance <- Inf }
      else if(abs(x - k) < 0.0000000001 | sigma^2 < 0.0000000001 | median <= s) { variance <- 0 }
      else
      {
        dss <- (median-s)/100.0
        ss <- 0
        pinf <- private$OUPPassageTimeProbabilityInf(x,k,omega,rho,mu,sigma)
        dnsty <- 0
        while(dnsty < 0.000001*pinf & ss < median)
        {
          ss <- ss+dss
          dnsty <- private$OUPPassageTimeDensity(0,x,ss,k,omega,rho,mu,sigma,0.05)
        }
        ss <- ss-dss
        t <- private$Laguerret()
        w <- private$Laguerrew()
        wghtfnct <- vector("double",41)
        density <- vector("double",41)
        tscale <- t[5]
        wscale <- 0
        variance <- 0
        i <- 0
        while(i < 41)
        {
          i <- i+1
          wghtfnct[i] <- exp(-t[i])
          t[i] <- t[i]*(median-s-ss)/tscale+ss
          w[i] <- w[i]*(median-s-ss)/tscale
          density[i] <- private$OUPPassageTimeDensity(0,x,t[i],k,omega,rho,mu,sigma,0.05)
          wscale <- wscale+w[i]*density[i]/wghtfnct[i]
          variance <- variance+w[i]*(t[i]-mean+s)^2*density[i]/wghtfnct[i]
        }
        variance <- variance/wscale
      }
      return(variance)
    },
    OUPPassageTimeDensity = function(s,x,t,k,omega,rho,mu,sigma,dt)
    {
      if(k == Inf | k == -Inf) { density <- 0 }
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
          if(k > lomeanx & k <= himeanx) { density <- 9 }
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
            if(k < mu & k < boundx) { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552) }
            else
            {
              if(lnlambda > 709.782712893384) { density <- u*exp(-vm2)/(2*2.50662827431001) }
              else { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*biggammab*dlnlambda*exp(lnlambda)/(2*1.77245385090552) }
            }
          }
          else
          {
            u <- ((k-mu)*exp(-2*rho*(t-s))-(x-mu)*exp(-rho*(t-s))+omega*((k-mu)*(exp(-rho*(t-s))-exp(-2*rho*(t-s)))-(x-k)*exp(-rho*(t-s))))*sigma^2/varx^1.5
            if(k > mu & k > boundx) { density <- u*exp(-vm2)/(2*2.50662827431001)+omega*dlnlambda*exp(lnlambda)*(1.77245385090552+smallgammab)/(2*1.77245385090552) }
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
    OUPPassageTimeProbability = function(s,x,t,k,omega,rho,mu,sigma)
    {
      if(k == Inf | k == -Inf) { probability <- 0 }
      else
      {
        if(x == k)
        {
          if(k == mu) { probability <- 0.5*(1+omega) }
          else
          {
            if(rho < 0.0000000001) { varx <- sigma^2*(t-s) }
            else { varx <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
            if(varx < 0.0000000001) { probability <- 0.5*(1+omega) }
            else
            {
              v2 <- 0.5*((k-mu)*(1-exp(-rho*(t-s))))^2/varx
              smallgamma <- private$GammaSmallOneHalf(v2)
              biggamma <- private$GammaBigOneHalf(v2)
              probability <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552)
            }
          }
        }
        else if(k == mu)
        {
          if(rho < 0.0000000001) { varx <- sigma^2*(t-s) }
          else { varx <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
          if(varx < 0.0000000001) { probability <- 0 }
          else
          {
            v2 <- 0.5*((k-x)*exp(-rho*(t-s)))^2/varx
            biggamma <- private$GammaBigOneHalf(v2)
            probability <- (1+omega)*biggamma/(2*1.77245385090552)
          }
        }
        else
        {
          meanx <- mu+(x-mu)*exp(-rho*(t-s))
          boundx <- mu+(k-mu)*exp(-rho*(t-s))-(x-k)*exp(-rho*(t-s))
          if(rho < 0.0000000001) { varx <- sigma^2*(t-s) }
          else { varx <- sigma^2*(1-exp(-2*rho*(t-s)))/(2*rho) }
          if(varx < 0.0000000001)
          {
            if(x > k & k > mu)
            {
              if(k < meanx) { probability <- 0 }
              else if(k == meanx) { probability <- 0.5}
              else { probability <- 1 }
            }
            else if(x > k & mu > k) { probability <- 0 }
            else if(x < k & k < mu)
            {
              if(k > meanx) { probability <- 0 }
              else if(k == meanx) { probability <- 0.5 }
              else { probability <- 1 }
            }
            else { probability <- 0 }
          }
          else
          {
            vm2 <- 0.5*(k-meanx)^2/varx
            vb2 <- 0.5*(k-boundx)^2/varx
            lnlambda <- 2*(k-mu)*(1-exp(-rho*(t-s)))*(x-k)*exp(-rho*(t-s))/varx
            smallgammam <- private$GammaSmallOneHalf(vm2)
            biggammam <- private$GammaBigOneHalf(vm2)
            smallgammab <- private$GammaSmallOneHalf(vb2)
            biggammab <- private$GammaBigOneHalf(vb2)
            if(x > k & k > mu)
            {
              if(k < meanx)
              {
                if(lnlambda > 709.782712893384) { probability <- biggammam/(2*1.77245385090552) }
                else { probability <- (biggammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552) }
              }
              else
              {
                if(lnlambda > 709.782712893384) { probability <- (1.77245385090552+smallgammam)/(2*1.77245385090552) }
                else { probability <- (1.77245385090552+smallgammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552) }
              }
            }
            else if(x > k & mu > k)
            {
              if(k > boundx) { probability <- (biggammam+omega*exp(lnlambda)*biggammab)/(2*1.77245385090552) }
              else { probability <- (biggammam+omega*exp(lnlambda)*(1.77245385090552+smallgammab))/(2*1.77245385090552) }
            }
            else if(x < k & k < mu)
            {
              if(k > meanx)
              {
                if(lnlambda > 709.782712893384) { probability <- biggammam/(2*1.77245385090552) }
                else { probability <- (biggammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552) }
              }
              else
              {
                if(lnlambda > 709.782712893384) { probability <- (1.77245385090552+smallgammam)/(2*1.77245385090552) }
                else { probability <- (1.77245385090552+smallgammam+omega*biggammab*exp(lnlambda))/(2*1.77245385090552) }
              }
            }
            else
            {
              if(k < boundx) { probability <- (biggammam+omega*exp(lnlambda)*biggammab)/(2*1.77245385090552) }
              else { probability <- (biggammam+omega*exp(lnlambda)*(1.77245385090552+smallgammab))/(2*1.77245385090552) }
            }
          }
        }
      }
      return(probability)
    },
    OUPPassageTimeProbabilityInf = function(x,k,omega,rho,mu,sigma)
    {
      if(k == Inf | k == -Inf) { pinf <- 0 }
      else
      {
        if(omega == 1) { pinf <- 1 }
        else if(k == mu) { pinf <- 0.5*(1+omega) }
        else if(sigma^2 < 0.0000000001)
        {
          if(x == k) { pinf <- 1 }
          else if(x > k)
          {
            if(k > mu) { pinf <- 1 }
            else { pinf <- omega }
          }
          else
          {
            if(k < mu) { pinf <- 1 }
            else { pinf <- omega }
          }
        }
        else
        {
          v2 <- rho*((k-mu)/sigma)^2
          smallgamma <- private$GammaSmallOneHalf(v2)
          biggamma <- private$GammaBigOneHalf(v2)
          if(x == k) { pinf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
          else if(x > k)
          {
            if(k > mu) { pinf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
            else { pinf <- (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552) }
          }
          else if(x < k)
          {
            if(k < mu) { pinf <- (1.77245385090552+smallgamma+omega*biggamma)/(2*1.77245385090552) }
            else { pinf <- (biggamma+omega*(1.77245385090552+smallgamma))/(2*1.77245385090552) }
          }
        }
      }
      return(pinf)
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
        while(10^sigdig*dgm >= gm & i < 2000)
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
        while(10^sigdig*dgm >= gm & i < n)
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
    Laguerret = function()
    {
      t <- vector("double",41)
      t[1] <- 0.034840064125239
      t[2] <- 0.183625087016639
      t[3] <- 0.451525058905516
      t[4] <- 0.838985610907935
      t[5] <- 1.34657594330787
      t[6] <- 1.97503995935016
      t[7] <- 2.72530671207309
      t[8] <- 3.59849894479794
      t[9] <- 4.59594284949542
      t[10] <- 5.71917969983466
      t[11] <- 6.96997971310658
      t[12] <- 8.3503584847076
      t[13] <- 9.86259639570966
      t[14] <- 11.5092614819208
      t[15] <- 13.2932363681025
      t[16] <- 15.2177500186414
      t[17] <- 17.2864152453923
      t[18] <- 19.5032731583782
      t[19] <- 21.872846064999
      t[20] <- 24.4002007458958
      t[21] <- 27.0910246000436
      t[22] <- 29.9517179152153
      t[23] <- 32.9895065670353
      t[24] <- 36.21258090697
      t[25] <- 39.6302686603258
      t[26] <- 43.2532526217178
      t[27] <- 47.0938482901412
      t[28] <- 51.1663631195495
      t[29] <- 55.4875691072365
      t[30] <- 60.0773363233172
      t[31] <- 64.9595008936776
      t[32] <- 70.1630847824147
      t[33] <- 75.7240620856209
      t[34] <- 81.6880101000193
      t[35] <- 88.1142662803787
      t[36] <- 95.0828121482822
      t[37] <- 102.706501661439
      t[38] <- 111.154921958725
      t[39] <- 120.707582561042
      t[40] <- 131.899754362277
      t[41] <- 146.110597447904

      return(t)
    },
    Laguerrew = function()
    {
      w <- vector("double",41)
      w[1] <- 8.63554570194484E-02
      w[2] <- 0.173332645592903
      w[3] <- 0.208568103938104
      w[4] <- 0.193350430207399
      w[5] <- 0.147724394333778
      w[6] <- 9.56297394958616E-02
      w[7] <- 5.31761202739599E-02
      w[8] <- 2.55882342423587E-02
      w[9] <- 1.06989176609361E-02
      w[10] <- 3.89523849800191E-03
      w[11] <- 1.2358867321078E-03
      w[12] <- 3.41684428543532E-04
      w[13] <- 8.22451802053959E-05
      w[14] <- 1.72110639995837E-05
      w[15] <- 3.12503224715398E-06
      w[16] <- 4.91094984483102E-07
      w[17] <- 6.65939575826955E-08
      w[18] <- 7.76479538935683E-09
      w[19] <- 7.75316169614965E-10
      w[20] <- 6.59856132000171E-11
      w[21] <- 4.76126275521546E-12
      w[22] <- 2.8950494277635E-13
      w[23] <- 1.47313265756315E-14
      w[24] <- 6.22367295054339E-16
      w[25] <- 2.1634252311101E-17
      w[26] <- 6.12366679384877E-19
      w[27] <- 1.39455738111393E-20
      w[28] <- 2.51965913863939E-22
      w[29] <- 3.5530244291388E-24
      w[30] <- 3.8348925676424E-26
      w[31] <- 3.09502899341052E-28
      w[32] <- 1.81545615594484E-30
      w[33] <- 7.47169937604288E-33
      w[34] <- 2.06348976691315E-35
      w[35] <- 3.60876174544819E-38
      w[36] <- 3.69668912238723E-41
      w[37] <- 1.98741198787653E-44
      w[38] <- 4.75752031803417E-48
      w[39] <- 3.87248620144704E-52
      w[40] <- 6.42131713584797E-57
      w[41] <- 5.91215183605371E-63

      return(w)
    },
    # private plot methods ----
    MeshWall = function(x,y,z)
    {
      m <- length(x)
      if(length(y) < m) { m <- length(y) }
      if(length(z) < m) { m <- length(z) }
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
        zv[i] <- 0
        xv[i+m] <- x[i]
        yv[i+m] <- y[i]
        zv[i+m] <- z[i]
      }
      i <- 0
      while(i < m-1)
      {
        i <- i+1
        xv[i+2*m] <- 0.5*(x[i]+x[i+1])
        yv[i+2*m] <- 0.5*(y[i]+y[i+1])
        zv[i+2*m] <- 0.25*(z[i]+z[i+1])
      }
      i <- 0
      while(i < m-1)
      {
        i <- i+1
        # right-hand side indexes are zero-based
        iv[i+(i-1)*(m-1)] <- i-1+2*m
        jv[i+(i-1)*(m-1)] <- i-1
        kv[i+(i-1)*(m-1)] <- i
        iv[i+1+(i-1)*(m-1)] <- i-1+2*m
        jv[i+1+(i-1)*(m-1)] <- i
        kv[i+1+(i-1)*(m-1)] <- i+m
        iv[i+2+(i-1)*(m-1)] <- i-1+2*m
        jv[i+2+(i-1)*(m-1)] <- i+m
        kv[i+2+(i-1)*(m-1)] <- i-1+m
        iv[i+3+(i-1)*(m-1)] <- i-1+2*m
        jv[i+3+(i-1)*(m-1)] <- i-1+m
        kv[i+3+(i-1)*(m-1)] <- i-1
      }
      return(list(xvertex=xv,yvertex=yv,zvertex=zv,ivertex=iv,jvertex=jv,kvertex=kv))
    },
    MeshVertical = function(x,y,ztop,zbottom)
    {
      m <- length(x)
      if(length(y) < m) { m <- length(y) }
      if(length(ztop) < m) { m <- length(ztop) }
      if(length(zbottom) < m) { m <- length(zbottom) }
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
        iv[i+(i-1)*(m-1)] <- i-1+2*m
        jv[i+(i-1)*(m-1)] <- i-1
        kv[i+(i-1)*(m-1)] <- i
        iv[i+1+(i-1)*(m-1)] <- i-1+2*m
        jv[i+1+(i-1)*(m-1)] <- i
        kv[i+1+(i-1)*(m-1)] <- i+m
        iv[i+2+(i-1)*(m-1)] <- i-1+2*m
        jv[i+2+(i-1)*(m-1)] <- i+m
        kv[i+2+(i-1)*(m-1)] <- i-1+m
        iv[i+3+(i-1)*(m-1)] <- i-1+2*m
        jv[i+3+(i-1)*(m-1)] <- i-1+m
        kv[i+3+(i-1)*(m-1)] <- i-1
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
