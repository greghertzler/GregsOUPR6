library(R6)
library(plotly)
library(stringr)

# roxygen ----
#' R6 Class implementing Maximum Likelihood estimation
#'
#' @description
#' Maximum Likelihood estimation of the Ornstein-Uhlenbeck Process, using a
#'  robust estimation algorithm and with testing for simple hypotheses.
#'
#' @details # Methods:
#'     rho, mu and sigma random
#'       LogLikelihood
#'       Estimates
#'       GoodnessOfFit
#'       LikelihoodRatioTest
#'
#' @details # Plots:
#'       PlotTimeSeries
#'       PlotEstimates
#'
#' @details # Arguments of functions:
#'       All arguments are optional in all functions.
#'       rho:    rate parameter 0<=rho<inf
#'       mu:     location parameter -inf<mu<inf
#'       sigma:  scale parameter -inf<sigma<inf
#'       tau:    vector of observed times -inf<tau<inf
#'       z:      vector of observed states -inf<z<inf
#'       df:     data frame containing columns for tau and z
#'       taucol: index of a column containing tau
#'       zcol:   index of a column containing z
#'       rhor:   constant to fix the rate parameter 0<=rhor<inf
#'       mur:    constant to fix the location parameter -inf<mur<inf
#'       sigmar: constant to fix scale parameter -inf<sigmar<inf
#'       rhos:   starting value for the rate parameter 0<=rhos<inf
#'       mus:    starting value for the location parameter -inf<mus<inf
#'       sigmas: starting value the scale parameter -inf<sigmas<inf
#'       lnLu:   unrestricted log likelihood -inf<lnLu<=0
#'       alphau: identifies distribution of lnLu 1<=alphau<=2
#'       lnLr:   restricted log likelihood -inf<lnLr<=lnLu
#'       alphar: identifies distribution of lnLr 1<=alphar<=2
#'       m1:     number of observations 0<m1=m-1
#'
#' @details # Using the methods:
#' Demonstration scripts are in files in the 'demo' directory. Identify a
#'  formula and in the console type:
#'
#'       demo(ML_FormulaName), or
#'       demo(ML_PlotFormulaName).
#'
#' @details # Discussion:
#' The Nelder-Mead algorithm is used to maximize a Log Likelihood derived from
#'  the Transition Density of the Ornstein-Uhlenbeck Process. The original
#'  algorithm is designed to minimize and has been called the amoeba algorithm.
#'  If spread across a bumpy surface, the amoeba pulls itself over bumps and out
#'  of hollows and flows down to the lowest level.  There it contracts around
#'  its center.  If spread across a flat surface, the amoeba shrinks without
#'  flowing.  Once the amoeba contracts or shrinks, it is spread across the
#'  surface opposite from where is came.  Again it flows, contracts and shrinks.
#'  Spreading, flowing, contracting and shrinking continue until the amoeba
#'  contracts around the same point or shrinks at the same level at least twice.
#'
#' In this implementation, the Nelder-Mead algorithm has been modified in five
#'  ways.
#'    1) It is modified to maximize;
#'    2) It can set parameters rho, mu and/or sigma to constants for simple
#'       hypothesis tests;
#'    3) It has a tie-breaking condition to prevent cycling or freezing;
#'    4) It can accelerate over long distances for searching flat likelihoods;
#'    5) It checks for log likelihoods greater than zero to prevent errors in
#'       case the data is faked or artificially smoothed.
#'
#' Once the unrestricted estimates are found, the arguments for restricting the
#'  parameters can be set for simple hypothesis tests. By default, starting
#'  values are calculated from the data and are then updated automatically. The
#'  arguments for starting values should seldom be needed.
#'
#' Data are entered as a data frame with columns for tau, a vector of times, and
#'  z, a vector of states. Arguments taucol and zcol are columns to extract from
#'  a data frame containing many columns.
#'
#' To read a csv file into a data frame, create a MaximumLikelihood object,
#'  estimate parameters and test the goodness-of-fit:
#'
#'       df <- read.csv("data/ML_OUP_SimulatedData.csv")
#'       ML <- MaximumLikelihood$new()
#'       ML$Estimates(df)
#'       ML$GoodnessOfFit()
#'
#' In this example, the 'data' directory is one level below the working
#'  directory. The taucol and zcol arguments are optional and will default to 1
#'  and 2. Continue on by restricting a parameter and performing a likelihood
#'  ratio test:
#'
#'       ML$Estimates(rhor=0.5)
#'       ML$LikelihoodRatioTest()
#'
#' The available data sets are documented. To see the documentation, type:
#'
#'       ?ML_
#'
#'  and then select from the drop-down list.
#'
#' Examples in R6 don't work the same way as other R modules.  There is only
#'  one example for an R6 object, not one for each function in the object.
#'  To run examples, the devtools::run_examples() works, but the R command
#'  example("MaximumLikelihood") doesn't.  You can copy commands to the clipboard,
#'  paste into the console and press Enter. Examples in this help and a simple
#'  example at the bottom can be run in this way.  A better alternative is demo().
#'
#' A better alternative is demo(), but this works sometimes, sometimes not.
#'  If you don't see the demos for this package, go to Files and demo.  You will
#'  see the demo names and then type something like demo(ML_Estimates).

# class ----
MaximumLikelihood <- R6::R6Class("MaximumLikelihood",
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
    #' Create a MaximumLikelihood object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new MaximumLikelihood object
    #' @examples
    #'   ML <- MaximumLikelihood$new()
    #'   ML$Estimates()
    #'   ML$GoodnessOfFit()
    #'   ML$Estimates(rhor=0.5,mur=-15,sigmar=15)
    #'   ML$LikelihoodRatioTest()
    initialize = function(OUP=NULL)
    {
      # pointer to container object ----
      if(!is.null(OUP) && class(OUP)[[1]] == "OUProcess") { private$OUP <- OUP }
      # arguments ----
      private$oup_params <- list(rho=0.459755421,mu=-22.8712176,sigma=13.7350886)
      private$oup_params_unrestr <- list(rhohat=0.459755421,muhat=-22.8712176,sigmahat=13.7350886)
      private$oup_params_restr <- list(rhor=0.5,mur=-15,sigmar=15)
      private$oup_params_start <- list(rhos=0.5,mus=-15,sigmas=15)
      # time series ----
      private$oup_stats <- list(lnLu=-447.788922,lnLr=-449.338874,alphar=0.9756147,m1=175)
      tau <- c(20,20.05,20.1,20.15,20.2,20.25,20.3,20.35,20.55,20.6,20.65,20.7,20.75,20.8,20.95,21,21.05,21.1,21.15,21.2,21.25,21.3,21.35,21.4,21.45,21.5,21.55,21.6,21.65,21.7,21.75,21.8,21.85,21.9,21.95,22,22.05,22.1,22.15,22.2,22.25,22.3,22.35,22.4,22.45,22.5,22.55,22.6,22.65,22.7,22.75,22.8,22.85,22.9,22.95,23,23.05,23.45,23.5,23.55,23.6,23.65,23.7,23.75,23.8,23.85,23.9,23.95,24.1,24.15,24.2,24.25,24.3,24.35,24.4,24.45,24.5,24.6,24.65,24.75,24.8,24.95,25,25.05,25.1,25.15,25.2,25.25,25.3,25.35,25.4,25.45,25.5,25.55,25.6,25.95,26,26.05,26.1,26.15,26.2,26.25,26.3,26.35,26.4,26.45,26.5,26.55,26.6,26.65,26.7,26.75,26.85,26.9,26.95,27,27.05,27.1,27.15,27.2,27.25,27.3,27.35,27.4,27.45,27.5,27.55,27.6,27.65,27.7,27.75,27.8,27.85,27.9,27.95,28,28.05,28.1,28.15,28.2,28.25,28.3,28.35,28.4,28.45,28.5,28.55,28.6,28.65,28.7,28.75,28.8,28.85,28.9,28.95,29,29.05,29.1,29.15,29.2,29.25,29.3,29.35,29.4,29.45,29.5,29.55,29.6,29.65,29.7,29.75,29.8,29.85,29.9,29.95,30)
      z <- c(30,27.49,24.08,23.45,24.72,21.5,22.15,21.05,27.9,27.71,33.01,30.55,32.42,26.43,24.37,23.22,17.02,17.89,20.04,17.9,19.52,16.78,16.13,16.91,13.92,15.02,13.04,8.19,10.62,12.04,10.39,11.94,14.62,6.61,1.7,1.71,-0.22,-2.69,1.69,-2.98,0.19,3.36,-0.3,-1.24,-4.35,-5.14,-14.63,-14.25,-19.83,-20.51,-20.98,-17.14,-18.94,-21.53,-25.73,-23.74,-24.57,-20.15,-22.54,-20.45,-17.59,-19.29,-21.85,-23.17,-23.8,-27.92,-26.87,-24.97,-17.32,-13.92,-13.49,-11.5,-9.23,-10.3,-9.54,-7.25,-4.04,-9.64,-9.42,-7.87,-6.86,-6.68,-9.03,-20.15,-19.21,-20.75,-21.07,-17.12,-9.99,-12.91,-12.19,-11.42,-8.96,-13.32,-17.8,-21.72,-26.66,-29.16,-28.77,-28.06,-22.83,-24.47,-29.79,-31.28,-34.05,-36.76,-38.97,-39.09,-32.89,-32.89,-31.25,-32.98,-36.19,-36.12,-40.74,-40.57,-42.75,-40.7,-38.5,-38.91,-40.66,-45.96,-48.73,-49.76,-48.19,-47.9,-47.73,-43.19,-41.62,-41.96,-42.98,-44.98,-41.98,-40.65,-36.21,-33.61,-33.73,-29.53,-34.07,-37.05,-33.84,-29.71,-29.24,-30.83,-26.19,-28.73,-30.27,-32.39,-30.3,-25.83,-29.39,-27.64,-25.63,-23.81,-23.34,-18.32,-16.53,-13.82,-13.74,-13.28,-16.98,-13,-13.03,-14.26,-16.94,-15.62,-13.54,-19.63,-14.06,-17.23,-17.6,-19.67,-19.87,-20.65,-14.42,-10.23)
      private$timeseries <- data.frame(tau=tau,z=z)
      private$timeseries_info <- list(tbeg=20,Ixbeg=1,tend=30,Ixend=176,dataname="Default",timename="time",statename="state",estimation="Unrestricted")
      # plot info ----
      plotfont <- list(family="Cambria",size=14)
      plotfile <- list(format="png",width=640,height=480)
      plottheme <- list(name="dark",opaque=1.0)
      private$plot_info <- list(plotfont=plotfont,plotfile=plotfile,plottheme=plottheme,plotlabels=TRUE)
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
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"ML")}
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
            private$loglikely <- NULL
            private$goodness <- NULL
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
            private$loglikely <- NULL
            private$goodness <- NULL
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
            private$loglikely <- NULL
            private$goodness <- NULL
          }
        }
        else { message("sigma not set.")}
      }
      return(private$oup_params)
    },
    #' Set OUP parameter restrictions
    #' @param rhor   rate parameter 0<=rhor<inf
    #' @param mur    location parameter -inf<mur<inf
    #' @param sigmar scale parameter -inf<sigmar<inf
    #' @return list(rhor,mur,sigmar)
    set_oup_params_restr = function(rhor=NULL,mur=NULL,sigmar=NULL)
    {
      rhosca <- private$extract_scalar(rhor)
      if(!is.null(rhor) & is.null(rhosca)) { message("rhor not set.")}
      musca <- private$extract_scalar(mur)
      if(!is.null(mur) & is.null(musca)) { message("mur not set.")}
      sigmasca <- private$extract_scalar(sigmar)
      if(!is.null(sigmar) & is.null(sigmasca)) { message("sigmar not set.")}
      rhoparam <- private$oup_params_restr[[1]]
      muparam <- private$oup_params_restr[[2]]
      sigmaparam <- private$oup_params_restr[[3]]
      hit <- FALSE
      if(!is.null(rhosca) & is.null(rhoparam) | is.null(rhosca) & !is.null(rhoparam)) { hit <- TRUE }
      else if(!is.null(musca) & is.null(muparam) | is.null(musca) & !is.null(muparam)) { hit <- TRUE }
      else if(!is.null(sigmasca) & is.null(sigmaparam) | is.null(sigmasca) & !is.null(sigmaparam)) { hit <- TRUE }
      else if(!is.null(rhosca) && (rhosca <= 0 & rhoparam > 0 | rhosca > 0 & rhosca != rhoparam)) { hit <- TRUE }
      else if(!is.null(musca) && musca != muparam) { hit <- TRUE }
      else if(!is.null(sigmasca) && sigmasca != sigmaparam) { hit <- TRUE }
      if(hit == TRUE)
      {
        private$oup_params_restr <- list(rhor=NULL,mur=NULL,sigmar=NULL)
        if(!is.null(rhosca))
        {
          if(rhosca < 0)
          {
            rhosca <- 0
            message("negative rhor set to zero.")
          }
          private$oup_params_restr$rhor <- rhosca
        }
        if(!is.null(musca)) { private$oup_params_restr$mur <- musca }
        if(!is.null(sigmasca)) { private$oup_params_restr$sigmar <- sigmasca }
        private$likelyratio <- NULL
        private$theta_r <- NULL
      }
      return(private$oup_params_restr)
    },
    #' Set OUP parameter starting values
    #' @param rhos   rate parameter 0<=rhos<inf
    #' @param mus    location parameter -inf<mus<inf
    #' @param sigmas scale parameter -inf<sigmas<inf
    #' @return list(rhos,mus,sigmas)
    set_oup_params_start = function(rhos=NULL,mus=NULL,sigmas=NULL)
    {
      if(!is.null(rhos))
      {
        sca <- private$extract_scalar(rhos)
        if(!is.null(sca))
        {
          if(sca < 0)
          {
            sca <- 0.0
            message("negative rhos set to zero.")
          }
          private$oup_params_start$rhos <- sca
        }
        else { message("rhos not set.")}
      }
      if(!is.null(mus))
      {
        sca <- private$extract_scalar(mus)
        if(!is.null(sca)) { private$oup_params_start$mus <- sca }
        else { message("mus not set.") }
      }
      if(!is.null(sigmas))
      {
        sca <- private$extract_scalar(sigmas)
        if(!is.null(sca)) { private$oup_params_start$sigmas <- sca }
        else { message("sigmas not set.") }
      }
      return(private$oup_params_start)
    },
    #' @description
    #' Set OUP statistics for estimation and hypothesis tests
    #' @param lnLu   unrestricted log likelihood -inf<lnLu<=0
    #' @param lnLr   restricted log likelihood -inf<lnLr<=lnLu
    #' @param alphar identifies probability of restricted log likelihood 0.5<=alphar<=1
    #' @param m1     number of observations 0<m1
    #' @return list(lnLu,lnLr,alphar,m1)
    set_oup_stats = function(lnLu=NULL,lnLr=NULL,alphar=NULL,m1=NULL)
    {
      if(!is.null(lnLu) & !is.null(lnLr))
      {
        lnLu <- private$extract_scalar(lnLu)
        lnLr <- private$extract_scalar(lnLr)
        if(!is.null(lnLu) & !is.null(lnLr))
        {
          if(lnLu < lnLr) { message("lnLu < lnLr.  lnLu and lnLr not set.") }
          else if(lnLu > 0) { message("lnLu > 0.  lnLu and lnLr not set.") }
          else if(lnLr > 0) { message("lnLr > 0.  lnLu and lnLr not set.") }
          else
          {
            private$oup_stats$lnLu <- lnLu
            private$oup_stats$lnLr <- lnLr
            private$likelyratio <- NULL
          }
        }
        else { message("lnLu and lnLr not set.") }
      }
      else if(!is.null(lnLu))
      {
        lnLu <- private$extract_scalar(lnLu)
        if(!is.null(lnLu))
        {
          lnL <- private$oup_stats[[2]]
          if(!is.null(lnL) && lnLu < lnL) { message("new lnLu < existing lnLr.  lnLu not set.") }
          else
          {
            private$oup_stats$lnLu <- lnLu
            private$likelyratio <- NULL
          }
        }
        else { message("lnLu not set.") }
      }
      else if(!is.null(lnLr))
      {
        lnLr <- private$extract_scalar(lnLr)
        if(!is.null(lnLr))
        {
          lnL <- private$oup_stats[[1]]
          if(!is.null(lnL) && lnLr > lnL) { message("new lnLr > existing lnLu.  lnLr not set.") }
          else
          {
            private$oup_stats$lnLr <- lnLr
            private$likelyratio <- NULL
          }
        }
        else { message("lnLr not set.") }
      }
      if(!is.null(alphar))
      {
        sca <- private$extract_scalar(alphar)
        if(!is.null(sca))
        {
          if(sca < 0.5)
          {
            message("alphar < 0.5 set to 0.5.")
            sca <- 0.5
          }
          if(sca > 1)
          {
            message("alphar > 1 set to 1.")
            sca <- 1
          }
          private$oup_stats$alphar <- sca
          private$likelyratio <- NULL
        }
        else { message("m1 not set.") }
      }
      if(!is.null(m1))
      {
        sca <- private$extract_scalar(m1)
        if(!is.null(sca))
        {
          if(sca < 1)
          {
            message("m1 < 1 set to 1.")
            sca <- 1
          }
          private$oup_stats$m1 <- sca
          private$likelyratio <- NULL
        }
        else { message("m1 not set.") }
      }
      return(private$oup_stats)
    },
    #' @description
    #' Set time series data for time tau and state z
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index a column containing z
    #' @return dataframe(tau,z)
    set_timeseries = function(df=NULL,taucol=NULL,zcol=NULL)
    {
      if(!is.null(df))
      {
        if(is.data.frame(df))
        {
          ncols <- ncol(df)
          if(ncols > 1)
          {
            nrows <- nrow(df)
            if(nrows > 2)
            {
              OK <- TRUE
              tc <- private$extract_scalar(taucol)
              zc <- private$extract_scalar(zcol)
              if(!is.null(tc) & !is.null(zc))
              {
                tc <- as.integer(tc)
                zc <- as.integer(zc)
                if(tc < 1)
                {
                  message("taucol too small,")
                  OK <- FALSE
                }
                else if(tc > ncols)
                {
                  message("taucol too large,")
                  OK <- FALSE
                }
                if(zc < 1)
                {
                  message("zcol too small,")
                  OK <- FALSE
                }
                else if(zc > ncols)
                {
                  message("zcol too large,")
                  OK <- FALSE
                }
                if(tc == zc)
                {
                  message("taucol equals zcol,")
                  OK <- FALSE
                }
              }
              else if(is.null(tc) & is.null(zc))
              {
                tc <- 1
                zc <- 2
                message("tau from column 1, z from column 2.")
              }
              else if(!is.null(tc))
              {
                tc <- as.integer(tc)
                if(tc < 1)
                {
                  message("taucol too small,")
                  OK <- FALSE
                }
                else if(tc > ncols)
                {
                  message("taucol too big,")
                  OK <- FALSE
                }
                else
                {
                  zc <- tc+1
                  if(zc > ncols) { zc <- tc-1}
                  message(paste(sep="","tau from column ",tc,", z from column ",zc,"."))
                }
              }
              else
              {
                zc <- as.integer(zc)
                if(zc < 1)
                {
                  message("zcol too small,")
                  OK <- FALSE
                }
                else if(zc > ncols)
                {
                  message("zcol too big,")
                  OK <- FALSE
                }
                else
                {
                  tc <- zc-1
                  if(tc < 1) { tc <- zc+1}
                  message(paste(sep="","tau from column ",tc,", z from column ",zc,"."))
                }
              }
              if(OK == TRUE)
              {
                df <- data.frame(df[,c(tc,zc)])
                df <- cleandata(df)
                nrows <- nrow(df)
                if(nrows > 2)
                {
                  i <- 1
                  while(i < nrows & OK == TRUE)
                  {
                    if(df[i+1,1] <= df[i,1])
                    {
                      message(paste(sep="","tau[",i+1,"]<=tau[",i,"] : ",df[i+1,1],"<=",df[i,1]))
                      OK <- FALSE
                    }
                    i <- i+1
                  }
                  if(OK == TRUE)
                  {
                    framenames <- colnames(df)
                    private$timeseries <- NULL
                    private$timeseries <- df[,c(1,2)]
                    private$timeseries_info <- list(tbeg=df[1,1],Ixbeg=1,tend=df[nrows,1],Ixend=nrows,dataname=framenames[2],timename=framenames[1],statename=framenames[2],estimation=NULL)
                    private$oup_params_unrestr <- list(rhohat=NULL,muhat=NULL,sigmahat=NULL)
                    private$oup_params_restr <- list(rhor=NULL,mur=NULL,sigmar=NULL)
                    private$oup_params_start <- list(rhos=NULL,mus=NULL,sigmas=NULL)
                    private$oup_stats <- list(lnLu=NULL,lnLr=NULL,alphar=NULL,m1=NULL)
                    private$loglikely <- NULL
                    private$theta_u <- NULL
                    private$theta_r <- NULL
                    private$goodness <- NULL
                  }
                  else { message("tau and z were not set.") }
                }
                else
                {
                  message("fewer than 3 rows remain after cleaning data,")
                  message("tau and z not set.")
                }
              }
              else { message("tau and z were not set.") }
            }
            else
            {
              message("fewer than 3 rows in data frame,")
              message("tau and z not set.")
            }
          }
          else
          {
            message("fewer than 2 columns in data frame,")
            message("tau and z not set.")
          }
        }
        else
        {
          message("data must be in a data.frame,")
          message("tau and z not set.")
        }
      }
      else
      {
        if(!is.null(taucol) | !is.null(zcol))
        {
          message("no data frame,")
          message("tau and z not set.")
        }
      }
      return(private$timeseries)
    },
    #' @description
    #' Set information for plotting times series and estimates
    #' @param tbeg       begin value for time axis
    #' @param tend       end value for time axis
    #' @param dataname   name for the data
    #' @param timename   name for time
    #' @param statename  name for state
    #' @param estimation description of estimation
    #' @return list(tbeg,Ixbeg,tend,Ixend,datename,timename,statename,estimation)
    set_timeseries_info = function(tbeg=NULL,tend=NULL,dataname=NULL,timename=NULL,statename=NULL,estimation=NULL)
    {
      tau <- private$timeseries[[1]]
      m <- length(tau)
      if(!is.null(tbeg))
      {
        if(tbeg == -Inf)
        {
          private$timeseries_info$tbeg <- tau[1]
          private$timeseries_info$Ixbeg <- 1
        }
        else if(tbeg == Inf)
        {
          private$timeseries_info$tbeg <- tau[m]
          private$timeseries_info$Ixbeg <- m
        }
        else
        {
          sca <- private$extract_scalar(tbeg)
          if(!is.null(sca))
          {
            if(sca < tau[1])
            {
              private$timeseries_info$tbeg <- tau[1]
              private$timeseries_info$Ixbeg <- 1
            }
            else if(sca > tau[m])
            {
              private$timeseries_info$tbeg <- tau[m]
              private$timeseries_info$Ixbeg <- m
            }
            else
            {
              hit <- FALSE
              i <- 0
              while(i < m & hit == FALSE)
              {
                i <- i+1
                if(sca <= tau[i])
                {
                  hit <- TRUE
                  private$timeseries_info$tbeg <- tau[i]
                  private$timeseries_info$Ixbeg <- i
                }
              }
            }
          }
          else { message("beg not set.") }
        }
      }
      if(!is.null(tend))
      {
        if(tend == Inf)
        {
          private$timeseries_info$tend <- tau[m]
          private$timeseries_info$Ixend <- m
        }
        else if(tend == -Inf)
        {
          private$timeseries_info$tend <- private$timeseries_info$tbeg
          private$timeseries_info$Ixend <- private$timeseries_info$Ixbeg
        }
        else
        {
          sca <- private$extract_scalar(tend)
          if(!is.null(sca))
          {
            if(sca > tau[m])
            {
              private$timeseries_info$tend <- tau[m]
              private$timeseries_info$Ixend <- m
            }
            else if(sca < private$timeseries_info$tbeg)
            {
              private$timeseries_info$tend <- private$timeseries_info$tbeg
              private$timeseries_info$Ixend <- private$timeseries_info$Ixbeg
            }
            else
            {
              hit <- FALSE
              i <- m+1
              while(i > 1 & hit == FALSE)
              {
                i <- i-1
                if(sca >= tau[i])
                {
                  hit <- TRUE
                  private$timeseries_info$tend <- tau[i]
                  private$timeseries_info$Ixend <- i
                }
              }
            }
          }
          else { message("end not set.") }
        }
      }
      if(!is.null(dataname))
      {
        chr <- private$extract_character(dataname)
        if(!is.null(chr)) { private$timeseries_info$dataname <- chr }
        else { message("data name not set.") }
      }
      if(!is.null(timename))
      {
        chr <- private$extract_character(timename)
        if(!is.null(chr)) { private$timeseries_info$timename <- chr }
        else { message("time name not set.") }
      }
      if(!is.null(statename))
      {
        chr <- private$extract_character(statename)
        if(!is.null(chr)) { private$timeseries_info$statename <- chr }
        else { message("state name not set.") }
      }
      if(!is.null(estimation))
      {
        chr <- private$extract_character(estimation)
        if(!is.null(chr)) { private$timeseries_info$estimation <- chr }
        else { message("estimation description not set.") }
      }
      return(private$timeseries_info)
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
    #' @param labels     title and parameters TRUE or FALSE
    #' @param who        identifier for object sending the parameters
    #' @return list(font,file,theme)
    set_plot_info = function(fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) & !is.null(private$OUP)) { private$OUP$send_plot_info(NULL,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,NULL,NULL,labels,"ML")}
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
    #' @return list(oup_params,oup_params_unrestr,oup_params_restr,oup_params_start,oup_stats,timeseries,timeseries_info,plot_info)
    get_all = function()
    {
      all <- list(oup_params = private$oup_params,
        oup_params_restr = private$oup_params_restr,
        oup_params_unrestr = private$oup_params_unrestr,
        oup_params_start = private$oup_params_start,
        oup_stats = private$oup_stats,
        timeseries = private$timeseries,
        timeseries_info = private$timeseries_info,
        plot_info = private$plot_info)
      return(all)
    },
    #' @description
    #' Get OUP parameters
    #' @return list(rho,mu,sigma)
    get_oup_params = function() { return(private$oup_params) },
    #' @description
    #' Get OUP parameters unrestricted estimates
    #' @return list(rhohat,muhat,sigmahat)
    get_oup_params_unrestr = function() { return(private$oup_params_unrestr) },
    #' @description
    #' Get OUP parameter restrictions
    #' @return list(rhor,mur,sigmar)
    get_oup_params_restr = function() { return(private$oup_params_restr) },
    #' @description
    #' Get OUP parameter starting values
    #' @return list(rhos,mus,sigmas)
    get_oup_params_start = function() { return(private$oup_params_start) },
    #' @description
    #' Get OUP statistics for estimation and hypothesis tests
    #' @return list(lnLu,lnLr,alphar,m1)
    get_oup_stats = function() { return(private$oup_stats) },
    #' @description
    #' Get time series data for time tau and state z
    #' @return dataframe(tau,z)
    get_timeseries = function() { return(private$timeseries) },
    #' @description
    #' Get information for times series
    #' @return list(tbeg,Ixbeg,tend,Ixend,datename,timename,statename,estimation)
    get_timeseries_info = function() { return(private$timeseries_info) },
    #' @description
    #' Get information for plotting
    #' @return list(font,file,theme,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)
    get_plot_colors = function() { return(private$plot_colors) },
    # public calculate methods ----
    #' @description
    #' Calculate the log likelihood of estimates
    #' @param rho    rate parameter 0<=rho<inf
    #' @param mu     location parameter -inf<mu<inf
    #' @param sigma  scale parameter -inf<sigma<inf
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index of a column containing z
    #' @return list(lnL)
    LogLikelihood = function(rho=NULL,mu=NULL,sigma=NULL,df=NULL,taucol=NULL,zcol=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_timeseries(df,taucol,zcol)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      # calculate ----
      loglikely <- private$loglikely
      if(is.null(loglikely))
      {
        loglikely <- private$OUPLogLikelihood(tau,z,rho,mu,sigma)
        private$loglikely <- loglikely
      }
      return(list(lnL=loglikely))
    },
    #' @description
    #' Calculate unrestricted or restricted maximum likelihood estimates
    #' @param df      data frame containing columns for tau and z
    #' @param taucol  index of a column containing tau
    #' @param zcol    index of a column containing z
    #' @param rhor    constant to fix the rate parameter 0<=rhor<inf
    #' @param mur     constant to fix the location parameter -inf<mur<inf
    #' @param sigmar  constant to fix scale parameter -inf<sigmar<inf
    #' @param rhos    starting value for the rate parameter 0<=rhos<inf
    #' @param mus     starting value for the location parameter -inf<mus<inf
    #' @param sigmas  starting value the scale parameter -inf<sigmas<inf
    #' @param plotit  TRUE or FALSE
    #' @return list(rhohat,muhat,sigmahat,lnLu,ku,alphau,m1) or list(rhor,mur,sigmar,lnLr,kr,alphar,m1)
    Estimates = function(df=NULL,taucol=NULL,zcol=NULL,rhor=NULL,mur=NULL,sigmar=NULL,rhos=NULL,mus=NULL,sigmas=NULL,plotit=TRUE)
    {
      # set / get ----
      self$set_timeseries(df,taucol,zcol)
      self$set_oup_params_start(rhos,mus,sigmas)
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      rhos <- private$oup_params_start[[1]]
      mus <- private$oup_params_start[[2]]
      sigmas <- private$oup_params_start[[3]]
      if(is.null(rhor) & is.null(mur) & is.null(sigmar)) { theta <- private$theta_u }
      else
      {
        self$set_oup_params_restr(rhor,mur,sigmar)
        rhor <- private$oup_params_restr[[1]]
        mur <- private$oup_params_restr[[2]]
        sigmar <- private$oup_params_restr[[3]]
        theta <- private$theta_r
      }
      # calculate ----
      if(is.null(theta))
      {
        nmstart <- private$OUPNMStart(rhor,mur,sigmar,rhos,mus,sigmas,tau,z)
        nm <- private$NelderMead(nmstart[[1]],nmstart[[2]],tau,z,nmstart[[4]])
        # reset and set ----
        rho <- nm[[1]][1]
        mu <- nm[[1]][2]
        sigma <- nm[[1]][3]
        lnL <- nm[[2]]
        nk <- nmstart[[3]]
        m1 <- nmstart[[4]]-1
        alpha <- 0
        i <- 0
        while(i < m1)
        {
          i <- i+1
          alpha <- alpha+1+exp(-2*rho*(tau[i+1]-tau[i]))
        }
        alpha <- 0.5*alpha/m1
        if(nmstart[[2]][1] != 0 & nmstart[[2]][2] != 0 & nmstart[[2]][3] != 0)
        {
          self$set_oup_stats(lnLu=lnL,NULL,NULL,m1=m1)
          self$set_timeseries_info(NULL,NULL,NULL,NULL,NULL,"Unrestricted")
          private$oup_params_start$rhos <- rho
          private$oup_params_start$mus <- mu
          private$oup_params_start$sigmas <- sigma
          theta <- list(rhohat=rho,muhat=mu,sigmahat=sigma,lnLu=lnL,ku=nk,alphau=alpha,m1=m1)
          private$theta_u <- theta
          private$oup_params_unrestr <- list(rhohat=rho,muhat=mu,sigmahat=sigma)
        }
        else
        {
          self$set_oup_stats(NULL,lnLr=lnL,alphar=alpha,m1=m1)
          self$set_timeseries_info(NULL,NULL,NULL,NULL,NULL,"Restricted")
          if(nmstart[[2]][1] != 0) { theta <- list(rhobar=rho) }
          else { theta <- list(rhor=rho) }
          if(nmstart[[2]][2] != 0) { theta <- append(theta,list(mubar=mu)) }
          else { theta <- append(theta,list(mur=mu)) }
          if(nmstart[[2]][3] != 0) { theta <- append(theta,list(sigmabar=sigma)) }
          else { theta <- append(theta,list(sigmar=sigma)) }
          theta <- append(theta,list(lnLr=lnL))
          theta <- append(theta,list(kr=nk))
          theta <- append(theta,list(alphar=alpha))
          theta <- append(theta,list(m1=m1))
          private$theta_r <- theta
        }
      }
      self$set_oup_params(rho=theta[[1]],mu=theta[[2]],sigma=theta[[3]])
      private$loglikely <- NULL
      private$goodness <- NULL
      # plot ----
      if(plotit == TRUE) { print(self$PlotEstimates(NULL,NULL,NULL)) }

      return(theta)
    },
    #' @description
    #' Calculate the goodness of fit compared with invariant and scaled brownian motion estimates
    #' @param rho    rate parameter 0<=rho<inf
    #' @param mu     location parameter -inf<mu<inf
    #' @param sigma  scale parameter -inf<sigma<inf
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index a column containing z
    #' @return list(Inv,SBM) with Inv=list(R2,Pval) and SBM=list(R2,Pval)
    GoodnessOfFit = function(rho=NULL,mu=NULL,sigma=NULL,df=NULL,taucol=NULL,zcol=NULL)
    {
      # set / get ----
      self$set_oup_params(rho,mu,sigma)
      self$set_timeseries(df,taucol,zcol)
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      # calculate ----
      goodness <- private$goodness
      if(is.null(goodness))
      {
        # unrestricted
        lnLu <- private$OUPLogLikelihood(tau,z,rho,mu,sigma)
        # invariant
        m <- length(tau)
        ave <- 0
        for(i in 1:(m-1)) { ave <- ave+z[i+1] }
        ave <- ave/(m-1)
        sumsq <- 0
        for(i in 1:(m-1)) { sumsq <- sumsq+(z[i+1]-ave)^2 }
        sumsq <- sumsq/(m-1)
        lnLr <- -0.5*(m-1)*(log(2*3.14159265358979*sumsq)+1)
        Upsilon2 <- 2*(lnLu-lnLr)
        if(Upsilon2 < 0)
        {
          message("parameters do not fit the data.")
          message("Inv R2 is negative.")
        }
        R2 <- 1-exp(-log(2)*Upsilon2/(m-1))
        Pval <- private$GammaBigRatio(0.5*(m-1),0.5*Upsilon2)
        Inv <- list(R2=R2,Pval=Pval)
        # scaled brownian motion
        sumsq <- 0
        sumlntau <- 0
        for(i in 1:(m-1))
        {
          sumsq <- sumsq+(z[i+1]-z[i])^2/(tau[i+1]-tau[i])
          sumlntau <- sumlntau + log(tau[i+1]-tau[i])
        }
        sumsq <- sumsq/(m-1)
        lnLr <- -0.5*(m-1)*(log(2*3.14159265358979*sumsq)+1)-0.5*sumlntau
        Upsilon2 <- 2*(lnLu-lnLr)
        if(Upsilon2 < 0)
        {
          message("parameters do not fit the data.")
          message("SBM R2 is negative.")
        }
        R2 <- 1-exp(-log(2)*0.5*Upsilon2/(m-1))
        Pval <- private$GammaBigRatio(m-1,0.5*Upsilon2)
        SBM <- list(R2=R2,Pval=Pval)
        goodness <- list(Inv=Inv,SBM=SBM)
        private$goodness <- goodness
      }
      return(goodness)
    },
    #' @description
    #' Test for significant difference between two likelihoods
    #' @param lnLu   unrestricted log likelihood -inf<lnLu<=0
    #' @param lnLr   restricted log likelihood -inf<lnLr<=lnLu
    #' @param alphar identifies distribution of restricted log likelihood 0.5<=alphar<=1
    #' @param m1     number of observations 0<m1=m-1
    #' @return list(R2,Pval)
    LikelihoodRatioTest = function(lnLu=NULL,lnLr=NULL,alphar=NULL,m1=NULL)
    {
      # set / get ----
      self$set_oup_stats(lnLu,lnLr,alphar,m1)
      lnLu <- private$oup_stats[[1]]
      lnLr <- private$oup_stats[[2]]
      alphar <- private$oup_stats[[3]]
      m1 <- private$oup_stats[[4]]
      likelyratio <- private$likelyratio
      # calculate ----
      if(is.null(likelyratio))
      {
        if(is.null(lnLu) | is.null(m1))
        {
          message("no unrestricted log likelihood.")
          R2 <- 1
          Pval <- 0
        }
        else if(is.null(lnLr) | is.null(alphar) | is.null(m1))
        {
          message("no restricted log likelihood.")
          R2 <- 0
          Pval <- 1
        }
        else if(lnLu-lnLr < 0)
        {
          message("lnLu is less than lnLr.")
          R2 <- 0
          Pval <- 0
        }
        else
        {
          Upsilon2 <- 2*(lnLu-lnLr)
          R2 <- 1-exp(-log(2)*0.5*Upsilon2/(alphar*m1))
          Pval <- private$GammaBigRatio(alphar*m1,0.5*Upsilon2)
        }
        likelyratio <- list(R2=R2,Pval=Pval)
        private$likelyratio <- likelyratio
      }
      return(likelyratio)
    },
    # public plot methods ----
    #' @description
    #' Plot time series
    #' @param df      data frame containing columns for tau and z
    #' @param taucol  index of a column containing tau
    #' @param zcol    index a column containing z
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @return plot
    PlotTimeSeries = function(df=NULL,taucol=NULL,zcol=NULL,tbeg=NULL,tend=NULL,title=NULL,xaxis=NULL,yaxis=NULL)
    {
      # set / get ----
      self$set_timeseries(df,taucol,zcol)
      self$set_timeseries_info(tbeg,tend,NULL,NULL,NULL,NULL)
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      Ixbeg <- private$timeseries_info[[2]]
      Ixend <- private$timeseries_info[[4]]
      dataname <- private$timeseries_info[[5]]
      timename <- private$timeseries_info[[6]]
      statename <- private$timeseries_info[[7]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      tauplot <- tau[Ixbeg:Ixend]
      zplot <- z[Ixbeg:Ixend]
      # plot ----
      # OUP_ML_TimeSeries2D
      if(labels == TRUE)
      {
        if(is.null(title)) { title <- dataname }
      }
      else if(is.null(title)) { title <- "" }
      if(is.null(xaxis))
      {
        bsym <- "<span style='font-family: Symbol;'>"
        esym <- "</span>"
        xaxis <- paste(sep="",bsym,"<i>t</i>",esym," (",timename,")")
      }
      if(is.null(yaxis)) { yaxis <- paste(sep="","<i>z</i> (",statename,")") }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      zline <- list(color=gry$e,width=2,dash="dot")
      zmarker <- list(color=gry$e,size=8)
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_ML_TimeSeries2D")
      fig <- plot_ly()  %>%
        add_trace(.,type="scatter",x=tauplot,y=zplot,name="<i>z</i>",hoverinfo="x+y+name",mode="lines+markers",line=zline,marker=zmarker) %>%
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    },
    #' @description
    #' Plot time series with estimated means and variances
    #' @param tbeg    begin value for time axis
    #' @param tend    end value for time axis
    #' @param title   text for plot title
    #' @param xaxis   text for x-axis label
    #' @param yaxis   text for y-axis label
    #' @return plot
    PlotEstimates = function(tbeg=NULL,tend=NULL,title=NULL,xaxis=NULL,yaxis=NULL)
    {
      # set / get ----
      self$set_timeseries_info(tbeg,tend,NULL,NULL,NULL,NULL)
      Ixbeg <- private$timeseries_info[[2]]
      Ixend <- private$timeseries_info[[4]]
      dataname <- private$timeseries_info[[5]]
      timename <- private$timeseries_info[[6]]
      statename <- private$timeseries_info[[7]]
      estimation <- private$timeseries_info[[8]]
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      tauplot <- tau[Ixbeg:Ixend]
      taumv <- tau[Ixbeg:(Ixend-1)]
      zplot <- z[Ixbeg:Ixend]
      m <- length(tauplot)
      # calculate ----
      mean <- vector("double",m-1)
      variance <- vector("double",m-1)
      resid <- vector("double",m-1)
      if(rho < 0.0000000001)
      {
        for(i in 1:(m-1))
        {
          mean[i] <- zplot[i]
          variance[i] <- sigma^2*(tauplot[i+1]-tauplot[i])
          if(abs(sigma) < 0.0000000001) { resid[i] <- 0 }
          else { resid[i] <- (zplot[i+1]-mean[i])/variance[i]^0.5 }
        }
      }
      else
      {
        for(i in 1:(m-1))
        {
          mean[i] <- mu+(zplot[i]-mu)*exp(-rho*(tauplot[i+1]-tauplot[i]))
          variance[i] <- sigma^2/(2*rho)*(1-exp(-2*rho*(tauplot[i+1]-tauplot[i])))
          if(abs(sigma) < 0.0000000001) { resid[i] <- 0 }
          else { resid[i] <- (zplot[i+1]-mean[i])/variance[i]^0.5 }
        }
      }
      # plot ----
      # OUP_ML_Estimates2D
      bsml <- "<span style='font-size: 10pt;'>"
      esml <- "</span>"
      bsym <- "<span style='font-family: Symbol;'>"
      esym <- "</span>"
      if(labels == TRUE)
      {
        syms <- paste(sep="",bsml,"(",bsym,"<i>r</i>=",esym,format(rho,digits=4),",",bsym,"<i>m</i>=",esym,format(mu,digits=4),",",bsym,"<i>s</i>=",esym,format(sigma,digits=4),")",esml)
        if(is.null(title)) { title <- estimation }
        if(is.null(xaxis)) { xaxis <- paste(sep="",bsym,"<i>t</i>",esym," (",timename,")<br>",syms) }
        else{ xaxis <- paste(sep="",xaxis,"<br>",syms) }
      }
      else
      {
        if(is.null(title)) { title <- "" }
        if(is.null(xaxis)) { xaxis <- paste(sep="",bsym,"<i>t</i>",esym," (",timename,")") }
      }
      if(is.null(yaxis)) { yaxis <- paste(sep="","<i>z</i> (",statename,"), <i>G</i>, <i>H</i><sup>2</sup>, <i>v</i>") }
      lookup <- list(text=title,yref="container",y=0.95)
      lookdown <- list(text=xaxis)
      lookleft <- list(text=yaxis)
      horz <- list(title=lookdown,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      vert <- list(title=lookleft,color=font$color,linewidth=1,showgrid=FALSE,ticks="outside",zeroline=FALSE)
      zmarker <- list(color=gry$e,size=4)
      meanline <- list(color=cyn$e,width=2,dash="dot")
      meanmarker <- list(color=cyn$b,size=6,line=list(width=2,color=cyn$e))
      varianceline <- list(color=mgn$e,width=2,dash="dot")
      variancemarker <- list(color=mgn$d,size=6,symbol="cross")
      residline <- list(color=red$e,width=2,dash="dot")
      residmarker <- list(color=red$d,size=6,symbol="x")
      imageoptions <- list(format=file$format,width=file$width,height=file$height,filename="OUP_ML_Estimates2D")
      legendpos <- list(orientation="h",x=1.0,y=1.05,xanchor="right")
      fig <- plot_ly()  %>%
        add_trace(.,type="scatter",x=tauplot,y=zplot,name="<i>z</i>",mode="markers",marker=zmarker) %>%
        add_trace(.,type="scatter",x=taumv,y=mean,name="<i>G</i>",mode="lines+markers",line=meanline,marker=meanmarker) %>%
        add_trace(.,type="scatter",x=taumv,y=variance,name="<i>H</i><sup>2</sup>",mode="lines+markers",line=varianceline,marker=variancemarker,visible="legendonly") %>%
        add_trace(.,type="scatter",x=taumv,y=resid,name="<i>v</i>",mode="lines+markers",line=residline,marker=residmarker,visible="legendonly") %>%
        config(.,toImageButtonOptions=imageoptions) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    OUP = NULL,
    # private input fields ----
    oup_params = NULL,
    oup_params_unrestr = NULL,
    oup_params_restr = NULL,
    oup_params_start = NULL,
    oup_stats = NULL,
    timeseries = NULL,
    timeseries_info = NULL,
    plot_info = NULL,
    # private output fields ----
    loglikely = NULL,
    theta_u = NULL,
    theta_r = NULL,
    goodness = NULL,
    likelyratio = NULL,
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
      }
      return(list(red=red,ylw=ylw,grn=grn,cyn=cyn,blu=blu,mgn=mgn,gry=gry,background=background,font=font))
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
        # if(is.list(input)) { input <- input[[1]] }
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
    cleandata = function(df)
    {
      if(!is.null(df))
      {
        if(is.data.frame(df))
        {
          m <- nrow(df)
          if(m > 2)
          {
            n <- ncol(df)
            if(n == 2)
            {
              framenames <- colnames(df)
              col1 <- vector("double",m)
              col2 <- vector("double",m)
              if(!is.numeric(df[,1]))
              {
                message(paste("column",framenames[1], "contains text entries."))
                suppressWarnings(df[,1] <- as.numeric(df[,1]))
                }
              if(!is.numeric(df[,2]))
              {
                message(paste("column",framenames[2], "contains text entries."))
                suppressWarnings(df[,2] <- as.numeric(df[,2]))
              }
              ix <- 0
              i <- 0
              while(i < m)
              {
                i <- i+1
                if(is.finite(df[i,1]) & !is.na(df[i,1]) & is.finite(df[i,2]) & !is.na(df[i,2]))
                {
                  ix <- ix+1
                  col1[ix] <- df[i,1]
                  col2[ix] <- df[i,2]
                }
              }
              col1 <- col1[1:ix]
              col2 <- col2[1:ix]
              df <- data.frame(col1,col2)
              colnames(df) <- c(framenames[1],framenames[2])
              return(df[order(df[[1]]),])
            }
            else
            {
              message("dataframe must have exactly two columns.")
              return(data.frame(tau=NULL,z=NULL))
            }
          }
          else
          {
            message("dataframe must have at least 3 rows.")
            return(data.frame(tau=NULL,z=NULL))
          }
        }
        else
        {
          message("argument is not a data frame")
          return(data.frame(tau=NULL,z=NULL))
        }
      }
      else
      {
        message("data frame is NULL.")
        return(data.frame(tau=NULL,z=NULL))
      }
    },
    # private calculate methods ----
    OUPNMStart = function(rhor,mur,sigmar,rhos,mus,sigmas,tau,z)
    {
      theta <- vector("double",3)
      steps <- vector("double",3)
      m <- length(tau)
      nk <- 0
      # sigma
      if(!is.null(sigmar))
      {
        steps[3] <- 0
        theta[3] <- sigmar
      }
      else
      {
        steps[3] <- 1
        nk <- nk+1
        if(!is.null(sigmas))
        {
          theta[3] <- sigmas
        }
        else
        {
          theta[3] <- 0
          for(i in 1:(m-1))
          {
            theta[3] <- theta[3]+(z[i+1]-z[i])^2/(tau[i+1]-tau[i])
          }
          theta[3] <- (theta[3]/(m-1))^0.5
        }
      }
      # mu
      if(!is.null(mur))
      {
        steps[2] <- 0
        theta[2] <- mur
      }
      else
      {
        steps[2] <- 1
        nk <- nk+1
        if(!is.null(mus))
        {
          theta[2] <- mus
        }
        else
        {
          theta[2] <- 0
          for(i in 1:(m-1))
          {
            theta[2] <- theta[2]+z[i+1]
          }
          theta[2] <- theta[2]/(m-1)
        }
      }
      # rho
      if(!is.null(rhor))
      {
        steps[1] <- 0
        if(rhor <= 0)
        {
          theta[1] <- 0
          if(steps[2] == 1) { nk <- nk-1 }
          steps[2] <- 0
          theta[2] <- 0
        }
        else { theta[1] <- rhor }
      }
      else
      {
        steps[1] <- 1
        nk <- nk+1
        if(!is.null(rhos))
        {
          if(rhos <= 0) { theta[1] <- 0 }
          else { theta[1] <- rhos }
        }
        else
        {
          scalrho <- 0
          for(i in 1:(m-1))
          {
            scalrho <- scalrho+(z[i+1]-theta[2])^2
          }
          scalrho <- scalrho/(m-1)
          if(scalrho > 0) { theta[1] <- theta[3]^2/scalrho/2 }
          else { theta[1] <- 0 }
        }
      }
      return(list(theta=theta,steps=steps,nk=nk,m=m))
    },
    NelderMead = function(theta,steps,tau,z,m)
    {
      # algorithm parameters
      rho <- 1         # reflect
      epsilon <- 2     # expand
      mu <- 1.1        # move (not in the usual Nelder-Mead algorithm)
      chi <- 0.5       # contract
      sigma <- 0.5     # shrink
      iota <- 5        # steps increment (not in the usual Nelder-Mead algorithm)
      # iteration parameters
      sigdig <- 12
      cmax <- 9999
      tensig <- 1/10^sigdig
      message("Nelder-Mead,    0:0")
      # cement constant thetas and create index of thetas in the simplex
      k <- length(theta)
      Ix <- vector("integer",k)
      tj <- vector("double",k)
      tbar <- vector("double",k)
      tr <- vector("double",k)
      te <- vector("double",k)
      tm <- vector("double",k)
      tc <- vector("double",k)
      ts <- vector("double",k)
      nk <- 0
      for(i in 1:k)
      {
        if(steps[i] == 0)
        {
          tj[i] <- theta[i]
          tbar[i] <- theta[i]
          tr[i] <- theta[i]
          te[i] <- theta[i]
          tm[i] <- theta[i]
          tc[i] <- theta[i]
          ts[i] <- theta[i]
        }
        else
        {
          steps[i] <- theta[i]/iota
          if(steps[i] < 0.1) {steps[i] <- 0.1}
          nk <- nk+1
          Ix[nk] <- i
        }
      }
      # outer loop
      LnL <- -private$OUPNMLnL(tau,z,m,theta)
      if(nk > 0)
      {
        tplex <- matrix(0.0,nk,nk+1)
        Lplex <- vector("double",nk+1)
        Lprev <- 1.79769313486231E+308
        sign <- -1
        cnt <- 0
        starts <- 0
        while(abs(Lprev-LnL) > abs(LnL*tensig) & cnt < cmax)
        {
          starts <- starts+1
          sign <- -1*sign
          # initial simplex
          for(i in 1:nk)
          {
            tplex[i,1] <- theta[Ix[i]]
          }
          Lplex[1] <- LnL
          for(j in 2:(nk+1))
          {
            i <- 0
            while(i < j-2)
            {
              i <- i+1
              tplex[i,j] <- theta[Ix[i]]
              tj[Ix[i]] <- tplex[i,j]
            }
            tplex[j-1,j] <- theta[Ix[j-1]]+sign*steps[Ix[j-1]]
            tj[Ix[j-1]] <- tplex[j-1,j]
            i <- j-1
            while(i < nk)
            {
              i <- i+1
              tplex[i,j] <- theta[Ix[i]]
              tj[Ix[i]] <- tplex[i,j]
            }
            Lplex[j] <- -private$OUPNMLnL(tau,z,m,tj)
          }
          # inner loop
          tdev <- 1.79769313486231E+308
          Ldev <- 1.79769313486231E+308
          while(tdev > tensig & Ldev > tensig & cnt < cmax)
          {
            #   minimum and maximum function values
            Lmin <- Lplex[1]
            jmin <- 1
            Lmax <- Lplex[1]
            jmax <- 1
            for(j in 2:(nk+1))
            {
              if(Lplex[j] < Lmin)
              {
                Lmin <- Lplex[j]
                jmin <- j
              }
              else if(Lplex[j] > Lmax)
              {
                Lmax <- Lplex[j]
                jmax <- j
              }
            }
            #  penultimate function value
            Lpenult <- Lmin
            for(j in 1:(nk+1))
            {
              if(Lpenult < Lplex[j] & Lplex[j] < Lmax)
              {
                Lpenult <- Lplex[j]
              }
            }
            #   calculate centroid of all but Lmax
            for(i in 1:nk)
            {
              tbar[Ix[i]] <- 0
              j <- 0
              while(j < jmax-1)
              {
                j <- j+1
                tbar[Ix[i]] <- tbar[Ix[i]]+tplex[i,j]
              }
              j <- jmax
              while(j < nk+1)
              {
                j <- j+1
                tbar[Ix[i]] <- tbar[Ix[i]]+tplex[i,j]
              }
              tbar[Ix[i]] <- tbar[Ix[i]]/nk
            }
            Lbar <- -private$OUPNMLnL(tau,z,m,tbar)
            #    calculate reflection of Lmax
            for(i in 1:nk)
            {
              tr[Ix[i]] <- tbar[Ix[i]]+rho*(tbar[Ix[i]]-tplex[i,jmax])
            }
            Lr <- -private$OUPNMLnL(tau,z,m,tr)
            #    expansion
            if(Lr < Lmin)
            {
              for (i in 1:nk)
              {
                te[Ix[i]] <- tbar[Ix[i]]+epsilon*(tr[Ix[i]]-tbar[Ix[i]])
              }
              Le <- -private$OUPNMLnL(tau,z,m,te)
              if(Le < Lr)
              {
                for(i in 1:nk)
                {
                  tplex[i,jmax] <- te[Ix[i]]
                }
                Lplex[jmax] <- Le
              }
              else
              {
                for(i in 1:nk)
                {
                  tplex[i,jmax] <- tr[Ix[i]]
                }
                Lplex[jmax] <- Lr
              }
            }
            #    reflection
            else if(Lmin <= Lr & Lr < Lpenult)
            {
              for(i in 1:nk)
              {
                tplex[i,jmax] <- tr[Ix[i]]
                tm[Ix[i]] <- tbar[Ix[i]]+mu*(tplex[i,jmin]-tbar[Ix[i]])
              }
              Lplex[jmax] <- Lr
              #    move minimum
              Lm <- -private$OUPNMLnL(tau,z,m,tm)
              if(Lm < Lplex[jmin])
              {
                for(i in 1:nk)
                {
                  tplex[i,jmin] <- tm[Ix[i]]
                }
                Lplex[jmin] <- Lm
              }
            }
            #    outside contraction
            else if(Lr < Lmax)
            {
              for(i in 1:nk)
              {
                tc[Ix[i]] <- tbar[Ix[i]]+chi*(tr[Ix[i]]-tbar[Ix[i]])
              }
              Lc <- -private$OUPNMLnL(tau,z,m,tc)
              if(Lc < Lr)
              {
                for(i in 1:nk)
                {
                  tplex[i,jmax] <- tc[Ix[i]]
                }
                Lplex[jmax] <- Lc
              }
              #   shrink
              else
              {
                j <- 0
                while(j < jmin-1)
                {
                  j <- j+1
                  for(i in 1:nk)
                  {
                    ts[Ix[i]] <- tplex[i,jmin]+sigma*(tplex[i,j]-tplex[i,jmin])
                    tplex[i,j] <- ts[Ix[i]]
                  }
                  Lplex[j] <- -private$OUPNMLnL(tau,z,m,ts)
                }
                j <- jmin
                while(j < nk+1)
                {
                  j <- j+1
                  for(i in 1:nk)
                  {
                    ts[Ix[i]] <- tplex[i,jmin]+sigma*(tplex[i,j]-tplex[i,jmin])
                    tplex[i,j] <- ts[Ix[i]]
                  }
                  Lplex[j] <- -private$OUPNMLnL(tau,z,m,ts)
                }
              }
            }
            #    inside contraction
            else
            {
              for(i in 1:nk)
              {
                tc[Ix[i]] <- tbar[Ix[i]]+chi*(tplex[i,jmax]-tbar[Ix[i]])
              }
              Lc <- -private$OUPNMLnL(tau,z,m,tc)
              if(Lc < Lmax)
              {
                for(i in 1:nk)
                {
                  tplex[i,jmax] <- tc[Ix[i]]
                }
                Lplex[jmax] <- Lc
              }
              #   shrink
              else
              {
                j <- 0
                while(j < jmin-1)
                {
                  j <- j+1
                  for(i in 1:nk)
                  {
                    ts[Ix[i]] <- tplex[i,jmin]+sigma*(tplex[i,j]-tplex[i,jmin])
                    tplex[i,j] <- ts[Ix[i]]
                  }
                  Lplex[j] <- -private$OUPNMLnL(tau,z,m,ts)
                }
                j <- jmin
                while(j < nk+1)
                {
                  j <- j+1
                  for(i in 1:nk)
                  {
                    ts[Ix[i]] <- tplex[i,jmin]+sigma*(tplex[i,j]-tplex[i,jmin])
                    tplex[i,j] <- ts[Ix[i]]
                  }
                  Lplex[j] <- -private$OUPNMLnL(tau,z,m,ts)
                }
              }
            }
            #   deviations
            tdev <- 0
            Ldev <- 0
            for(j in 1:(nk+1))
            {
              for(i in 1:nk)
              {
                if(abs((tplex[i,j]-tbar[Ix[i]])/steps[Ix[i]]) > tdev)
                {
                  tdev <- abs((tplex[i,j]-tbar[Ix[i]])/steps[Ix[i]])
                }
              }
              if(abs(Lplex[j]/Lbar-1) > Ldev)
              {
                Ldev <- abs(Lplex[j]/Lbar-1)
              }
            }
            nstarts <- str_length(as.integer(starts))
            ncnt <- str_length(as.integer(cnt))
            back <- strrep("\b",nstarts+ncnt+2)
            cnt <- cnt+1
            message(paste(sep="",back,starts,":",cnt))
          }
          # new minimum
          Lprev <- LnL
          for(i in 1:nk)
          {
            theta[Ix[i]] <- tplex[i,jmin]
            steps[Ix[i]] <- theta[Ix[i]]/iota
            if(steps[i]< 0.1) {steps[i] <- 0.1}
          }
          LnL <- Lplex[jmin]
        }
      }
      LnL <- -LnL
      return(list(theta=theta,lnL=LnL))
    },
    OUPNMLnL = function(tau,z,m,theta)
    {
      LnL <- 0
      if(abs(theta[3]) < 0.00001) { theta[3] <- 0.000001 }
      if(abs(theta[1]) < 0.0000000001)
      {
        for(i in 1:(m-1))
        {
          mean <- z[i]
          variance <- theta[3]^2*(tau[i+1]-tau[i])
          lnu <- -0.5*log(2*3.14159265358979*variance)
          v <- 0.5*(z[i+1]-mean)^2/variance
          LnL <- LnL+lnu-v
        }
      }
      else
      {
        for(i in 1:(m-1))
        {
          mean <- theta[2]+(z[i]-theta[2])*exp(-theta[1]*(tau[i+1]-tau[i]))
          variance <- theta[3]^2/(2*theta[1])*(1-exp(-2*theta[1]*(tau[i+1]-tau[i])))
          lnu <- -0.5*log(2*3.14159265358979*variance)
          v <- 0.5*(z[i+1]-mean)^2/variance
          LnL <- LnL+lnu-v
        }
      }
      return(LnL)
    },
    OUPLogLikelihood = function(tau,z,rho,mu,sigma)
    {
      if(abs(sigma) < 0.000001) { sigma <- 0.000001 }
      if(abs(rho) < 0.0000000001)
      {
        loglikelihood <- 0
        m <- length(tau)
        for(i in 1:(m-1))
        {
          mean <- z[i]
          variance <- sigma^2*(tau[i+1]-tau[i])
          lnu <- -0.5*log(2*3.14159265358979*variance)
          v <- 0.5*(z[i+1]-mean)^2/variance
          loglikelihood <- loglikelihood+lnu-v
        }
      }
      else
      {
        loglikelihood <- 0
        m <- length(tau)
        for(i in 1:(m-1))
        {
          mean <- mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]))
          variance <- sigma^2/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])))
          lnu <- -0.5*log(2*3.14159265358979*variance)
          v <- 0.5*(z[i+1]-mean)^2/variance
          loglikelihood <- loglikelihood+lnu-v
        }
      }
      return(loglikelihood)
    },
    GammaComplete = function(a)
    {
      # small to medium a
      if(a < 50)
      {
        # split a into integer portion, p, and remainder, r
        p <- floor(a)
        r <- a - p
        # check for a equal to zero or negative integer
        if(r == 0 & p < 1) { gamma <- NaN }
        else
        {
          # special cases
          if(2*a-floor(2*a) == 0)
          {
            if(r == 0)
            {
              gamma <- 1
              j <- 0
              while(p-1-j > 1e-15)
              {
                j <- j+1
                gamma <- gamma*j
              }
            }
            else
            {
              gamma <- 1.77245385090552
              if(a > 0)
              {
                j <- -0.5
                while(p-0.5-j > 1e-15)
                {
                  j <- j+1
                  gamma <- gamma*j
                }
              }
              else
              {
                j <- 0.5
                while(j-(p+0.5) > 1e-15)
                {
                  j <- j-1
                  gamma <- gamma/j
                }
              }
            }
          }
          # other cases
          else
          {
            sigdig <- 15
            cmax <- 2000
            # confluent hypergeometric in well-behaved region 1 < r < 2
            dgm <- 1
            gm <- dgm
            cnt <- 0
            while(cnt < cmax & 10^sigdig*dgm >= gm)
            {
              cnt <- cnt+1
              dgm <- dgm*199/(r+1+cnt)
              gm <- gm+dgm
            }
            gm <- 199^(r+1)/(r+1)*exp(-199)*gm
            # factorial calculations (gm is gamma of r+1)
            j <- r
            while(p-1+r-j > 1e-15)
            {
              j <- j+1
              gm <- gm*j
            }
            j <- 1+r
            while(j-(p+r) > 1e-15)
            {
              j <- j-1
              gm <- gm/j
            }
            gamma <- gm
          }
        }
      }
      else
      {
        gamma <- exp(GammaLn(a))
      }
      return(gamma)
    },
    GammaLn = function(a)
    {
      # positive a
      if(a <= 0) { gammaln <- NaN }
      # ln of GammaComplete for small to medium a
      else if(a < 50) { gammaln <- log(GammaComplete(a)) }
      # Poincaire expansion
      else
      {
        B2k <- vector("double",16)
        # absolute value of Bernoulli numbers (they alternate in sign)
        B2k[1] <- 0.166666666666667
        B2k[2] <- 3.33333333333333E-02
        B2k[3] <- 2.38095238095238E-02
        B2k[4] <- 3.33333333333333E-02
        B2k[5] <- 7.57575757575758E-02
        B2k[6] <- 0.253113553113553
        B2k[7] <- 1.16666666666667
        B2k[8] <- 7.0921568627451
        B2k[9] <- 54.9711779448622
        B2k[10] <- 529.124242424242
        B2k[11] <- 6192.1231884058
        B2k[12] <- 86580.2531135531
        B2k[13] <- 1425517.16666667
        B2k[14] <- 27298231.0678161
        B2k[15] <- 601580873.900642
        B2k[16] <- 15116315767.0922
        gammaln <- (a-0.5)*log(a)-a+0.918938533204673
        k <- 0
        while(k < 15)
        {
          k <- k+1
          dgm <- log(B2k[k])-log(2*k*(2*k-1))-(2*k-1)*log(a)
          gammaln <- gammaln+exp(dgm)
          k <- k+1
          dgm <- log(B2k[k])-log(2*k*(2*k-1))-(2*k-1)*log(a)
          gammaln <- gammaln-exp(dgm)
        }
      }
      return(gammaln)
    },
    GammaSmallRatio = function(a,x)
    {
      # zero or negative a
      if(a <= 0) { gammaratio <- NaN }
      # degenerate solution
      else if(x <= 0) { gammaratio <- 0 }
      # series expansion using logs
      else
      {
        s <- 15
        n <- 30000
        lnx <- log(x)
        sumlnxlna <- 0
        dgm <- 1
        gm <- dgm
        i <- 1
        while(10^s*dgm >= gm & i <= n)
        {
          lnxlna <- lnx-log(i+a)
          sumlnxlna <- sumlnxlna+lnxlna
          dgm <- exp(sumlnxlna)
          gm <- gm+dgm
          i <- i+1
        }
        lngm <- -log(a)+a*lnx-x+log(gm)-GammaLn(a)
        if(lngm <= -27.6310211159285) { gammaratio <- 0 }
        else if(lngm >= -1e-12) { gammaratio <- 1 }
        else { gammaratio <- exp(lngm) }
      }
      return(gammaratio)
    },
    GammaBigRatio = function(a,x)
    {
      # zero or negative a
      if(a <= 0) { gammaratio <- NaN }
      # degenerate solution
      else if(x <= 0) { gammaratio <- 1 }
      # by subtraction
      else { gammaratio <- 1-GammaSmallRatio(a, x) }
      return(gammaratio)
    }
 )
  # class end ----
)
# data roxygen ----

#' My data for the Ornstein-Uhlenbeck Process
#'
#' Data to estimate parameters, rho, mu and sigma, where rho is the rate of convergence,
#'  mu is the location and sigma is the scale.
#'
#' \itemize{
#'   \item tau: time variable
#'   \item z: state variable
#' }
#'
#' The data must be in a .csv (Comma Separated Values) file.  The first column should
#'  be times and the second column should be states of nature.  There can be more columns
#'  for times and states if you wish.  There can be blank entries.  The data will be
#'  cleaned and sorted by time before it is used.
#'
#' @docType data
#' @keywords datasets
#' @name MyData
#' @format csv file with at least 3 rows and 2 columns
NULL

#' Rates of convergence for the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate different rates of convergence, rho.
#'
#' \itemize{
#'   \item year: time variable in annual increments for all sample paths
#'   \item z1-z5: sample paths in sets of three, each set with the same pseudo-random shocks but different rates of convergence
#' }
#'
#' The rate of convergence, rho, determines the probability distribution of the
#'  estimated parameters and the correlation between two sets of parameters in
#'  hypothesis tests.  Small rho tends toward Brownian Motion, which does not
#'  converge.  Large rho tends toward a stationary or ergodic process which has
#'  converged every time it is observed.  In between is an Ornstein-Uhlenbeck
#'  Process which converges but has not yet converged.
#'
#' Parameters for Browian Motion have an Erlang distribution.  Parameters for
#'  a stationary or ergodic process have a Chi^2 distribution.  These distributions
#'  are special cases of a Gamma distribution.  In general, parameters for the
#'  Ornstein-Uhlenbeck Process have a Gamma distribution.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_Convergence
#' @format csv file with 177 rows and 16 columns
NULL

#' Observation intervals and the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate the effect of the observation interval on the
#'  rate of convergence, rho, and the scale, sigma.
#'
#' \itemize{
#'   \item year 1: time variable in annual increments
#'   \item z rho 0.5: sample paths with annual time increments and a rate of convergence of 0.5
#'   \item year 0.05: time variable in increments of 1/20 year
#'   \item z rho 10: sample paths with time increments of 1/20 year and a rate of convergence of 10
#' }
#'
#' This data has two sample paths.  Within a sample path, observation intervals are
#'  of equal length.  Between sample paths the observation intervals are of unequal
#'  lengths.  Otherwise, the sample paths have the same shocks and give equivalent
#'  estimates of the parameters.
#'
#'  Both sets of estimates have the same rho(t-s), where t-s is the observation interval.
#'  A common example is annual versus daily intervals; t-s=1 and rho=1/2 is equivalent to
#'  t-s=1/365 and rho=365/2.  Location parameter, mu, is the same in both sets of estimates,
#'  but scale parameter, sigma, is larger if rho is larger because the asymptotic variance,
#'  sigma^2/2rho, is the same in both sets of estimates.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_ObservationInterval
#' @format csv file with 201 rows and 4 columns
NULL

#' Sample sizes for the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate the effect of sample size on hypothesis tests
#'  for the Ornstein-Uhlenbeck Process.
#'
#' \itemize{
#'   \item year: time variable in annual increments for all sample paths
#'   \item xsmall: sample with 37 measurements or 36 observations
#'   \item small: sample with 2 x 37 measurements or 73 observations
#'   \item medium: sample with 8 x 37 measurements or 295 observations
#'   \item large: sample with 32 x 37 measurements or 1183 observations
#'   \item xlarge: sample with 128 x 37 measurements or 4735 observations
#' }
#'
#' This data starts with a very small sample and then replicates it.  Larger
#'  samples contain no more information than smaller samples and hypothesis
#'  tests for the Ornstein-Uhlenbeck Process do not change with the sample size.
#'  The usual Chi^2 tests, however, will reject all null hypotheses for medium
#'  to large sample sizes.
#'
#' If larger samples contain more information, hypothesis tests for the
#'  Ornstein-Uhlenbeck Process will reject more null hypotheses, but not all
#'  hypotheses as do Chi^2 tests.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_SampleSize
#' @format csv file with 4736 rows and 6 columns
NULL

#' Smoothed sample paths for the Ornstein-Uhlenbeck Process
#'
#' Simulated and smoothed data to demonstrate the effect of smoothing on
#'  parameter estimates.
#'
#' \itemize{
#'   \item year: time variable in annual increments
#'   \item Data: simulated sample path
#'   \item G uno-G diez: ten successively smoother paths
#' }
#'
#' Simulated data is smoothed ten times.  First, the parameters are estimated
#'  and the means calculated for each observation. Then the calculated means
#'  are used, as if they are data, in a subsequent estimation. Means are
#'  calculated again and used in the next estimation and so on ten times.
#'  The true rate of convergence, rho, and location, mu, are recovered from
#'  each estimation, but the scale, sigma, goes toward zero.
#'
#' By the ninth smoothing, the log of the likelihood becomes positive. Hence,
#'  the likelihood, as the anti-log, becomes greater than one. In other words,
#'  if the ninth and tenth smoothings were real samples, the probability of
#'  observing them would be greater than 100%. Further smoothings would increase
#'  this probability.  A small sigma is a tell-tale sign the data has been smoothed.
#'  A positive log likelihood is a sure sign.
#'
#' This is the best possible smoothing method. The model used for the smoothing
#'  is consistent with the data. Surprisingly, hypothesis tests and decision
#'  thresholds are not greatly affected.  However, passage times are calculated
#'  to be much larger and are completely unreliable.
#'
#' The best possible smoothing is unlikely. Any model used to smooth the data is
#'  probably not the Ornstein-Uhlenbeck Process.  The estimates will be wrong
#'  and the actual system will be much more uncertain. How much more uncertain
#'  is uncertain.
#'
#' Always use the raw data.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_SmoothedData
#' @format csv file with 177 rows and 11 columns
NULL

#' Unequal observation intervals and the Ornstein-Uhlenbeck Process
#'
#' Simulated data to demonstrate that observations are never 'missing'.
#'
#' \itemize{
#'   \item year 1: time variable in annual increments
#'   \item z rho 0.5: sample paths with annual time increments or longer and a rate of convergence of 0.5
#'   \item year 0.05: time variable in increments of 1/20 year
#'   \item z rho 10: sample paths with time increments of 1/20 year or longer and a rate of convergence of 10
#' }
#'
#' An observation consists of the initial state, the initial time, the terminal state
#'  and the terminal time.  The terminal time minus the initial time can differ among
#'  observations because each observation has its own mean and variance.  Longer
#'  observation intervals have smaller means and larger variances but, otherwise,
#'  parameters estimated from unequal observation intervals are statistically
#'  equivalent to parameters estimated from equal intervals.
#'
#' The usual time-series methods are special cases of estimating the Ornstein-Uhlenbeck
#'  Process.  They estimate location, mu, and scale, sigma, but not the rate of
#'  convergence, rho.  They eliminate rho by assuming weak stationarity.  The mean of
#'  each observation is assumed to have converged and lost its connection to its
#'  initial state.  The variance of each observation does not depend upon an initial
#'  state and may still be converging.  Observations will have equal variances if
#'  observation intervals are equal.
#'
#' For the Ornstein-Uhlenbeck Process, variances converge twice as fast as means.
#'  If variances are still converging, so are means.  Observations over equal time
#'  intervals will have equal variances but different means.  In other words, weak
#'  stationarity does not exist.  Any assumption of stationarity is strong stationarity
#'  in which both the means and variances have converged.
#'
#' Stationarity is an hypothesis to test, not an assumption to impose.  Goodness
#'  of fit for the Ornstein-Uhlenbeck Process tests for stationarity, and also
#'  for the other extreme of Brownian Motion.
#'
#' If you have a time series with measurements at sporadic times, don't fill in for
#'  observations that aren't missing.  If you must, first estimate the parameters.
#'  Then use the parameters to calculate means for the times of your choosing.
#'  Means are the maximum likelihood estimates of unobserved states of nature.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_UnequalIntervals
#' @format csv file with 151 rows and 4 columns
NULL

#' Water in Farm Dams in the Riverina of New South Wales
#'
#' Matlab generated water volumes in a 4000 cubic metre dam, with and without
#'  tree windbreaks to control evaporation, and a merino sheep flock size of
#'  1500 dry sheep equivalents.
#'
#' \itemize{
#'   \item Day: time variable in daily increments
#'   \item Without: water volumes in cubic metres without a windbreak
#'   \item With: water volumes in cubic metres with a windbreak
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_FarmDamsRiverina
#' @format csv file with 19313 rows and 3 columns
#' @author { Tim Capon \email{tim.capon.csiro.au}, Helena Clayton \email{helena.clayton@anu.edu.au}, Sally Thompson \email{sally.thompson@uwa.edu.au}, Greg Hertzler \email{ghertzlerau@gmail.com}, Philip Graham \email{phil@gramadvisory.com.au}, David Lindemayer \email{david.lindemayer@anu.edu} }
NULL

#' Farm Adaptation at Cootamundra, New South Wales
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S)
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMCootamundra
#' @format csv file with 50 rows and 5 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Narrendera, New South Wales
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S)
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMNarrendera
#' @format csv file with 50 rows and 5 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Temora, New South Wales
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S)
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMTemora
#' @format csv file with 50 rows and 5 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Soil Mineral Nitrogen with Stubble Management
#'
#' CSIRO Harden Long-Term Tillage Experiment
#'
#' \itemize{
#'   \item Year: time variable in sporadic increments
#'   \item Burn: kilograms per hectare of mineral soil nitrogen in the soil profile from 0 to 160 millimetres with stubble burning
#'   \item Bash: kilograms per hectare of mineral soil nitrogen in the soil profile from 0 to 160 millimetres with stubble bashing
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_SoilNitrogenHarden
#' @format csv file with 27 rows and 3 columns
#' @source https://data.csiro.au/collection/csiro:61293?q=population%20data&_st=keyword&_str=210&_si=7
NULL

#' Soil Water with Stubble Management
#'
#' CSIRO Harden Long-Term Tillage Experiment
#'
#' \itemize{
#'   \item Year: time variable in sporadic increments
#'   \item Burn: millimetres of water in the soil profile from 0 to 160 millimetres with stubble burning
#'   \item Bash: millimetres of water in the soil profile from 0 to 160 millimetres with stubble bashing
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_SoilWaterHarden
#' @format csv file with 21 rows and 3 columns
#' @source https://data.csiro.au/collection/csiro:61293?q=population%20data&_st=keyword&_str=210&_si=7
NULL

#' Farm Adaptation at Clare, South Australia
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMClare
#' @format csv file with 108 rows and 7 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Hawker, South Australia
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMHawker
#' @format csv file with 108 rows and 7 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Orroroo, South Australia
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMOrroroo
#' @format csv file with 108 rows and 7 columns
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Tree shelter belts in Tasmania
#'
#' CSIRO Perennial Prosperity Project
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item With trees: sheep gross margins per paddock with tree shelter belt
#'   \item Without trees: sheep gross margins per paddock without tree shelter belt
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Agric_Tas_TreeShelterBelts
#' @format csv file with 52 rows and 3 columns
#' @author { Tim Capon \email{tim.capon@csiro.au}, Daniel Mendham \email{daniel.mendham@csiro.au} }
NULL

#' Sea Level at Port Kembla
#'
#' Bureau of Meteorology Australian Baseline Sea Level Monitoring Project, with
#' measurements smoothed to eliminate variability.
#'
#' \itemize{
#'   \item Day in 2023: time variable in hourly increments
#'   \item Sea Level: metres above Tide Gauge Zero
#'   \item Water Temperature: degrees Celsius
#'   \item Air Temperature: degrees Celsius
#'   \item Barometric Pressure: hPa
#' }
#'
#' This is radically smoothed data.  Do not use this data for decisions under uncertainty.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_SeaLevel_PortKembla
#' @format csv file with 7985 rows and 5 columns
#' @source http://www.bom.gov.au/oceanography/projects/abslmp/data/data.shtml
NULL

#' Sunspot Numbers
#'
#' World Data Center SILSO, Royal Observatory of Belgium
#'
#' \itemize{
#'   \item DecimalYear: time variable in daily increments
#'   \item Average: Simple average of the daily total sunspot number over all days of each calendar month.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Climate_Sunspots
#' @format csv file with 3316 rows and 2 columns
#' @source https://doi.org/10.24414/qnza-ac80
NULL

#' Maximum daily temperatures at Cape Otway
#'
#' Bureau of Meteorology, raw (unhomogenized) data.
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_CapeOtway
#' @format csv file with 8888 rows and 2 columns
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 090015
NULL

#' Maximum daily temperatures at Darwin
#'
#' Bureau of Meteorology, raw (Unhomogenized) data.
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_Darwin
#' @format csv file with 9127 rows and 2 columns
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 014040
NULL

#' Maximum daily temperatures at Tennant Creek
#'
#' Bureau of Meteorology, raw (Unhomogenized) data.
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_TennantCreek
#' @format csv file with 9115 rows and 2 columns
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 015135
NULL

#' Water Supply for Irrigated Agriculture in Eastern Australia
#'
#' CSIRO Agricultural water supply data for national ecosystem services accounts.
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item WaterSupply: irrigation water supplied in Australia in megalitres
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_IrrigationWaterSupply
#' @format csv file with 13 rows and 2 columns
#' @source https://doi.org/10.25919/7jj7-8826
NULL

#' Kangaroo Population and Harvest in South Australia
#'
#' Government of South Australia, Department of Environment, population surveys and
#'  meat processor records of harvest aggregated across all regions of South Australia.
#'
#' \itemize{
#'   \item Year: time variable
#'   \item Red Pop: Estimated population of Red Kangaroos
#'   \item Red Harvest: Number of Red Kangaroos commercially harvested
#'   \item Grey Pop: Estimated population of Western Grey Kangaroos
#'   \item Grey Harvest: Number of Western Grey Kangaroos commercially harvested
#'   \item Euro Pop: Estimated population of European Kangaroos
#'   \item Euro Harvest: Number of European Kangaroos commercially harvested
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_Kangaroos
#' @format csv file with 46 rows and 7 columns
#' @source https://www.environment.sa.gov.au/topics/animals-and-plants/sustainable-use-of-animals-and-plants/kangaroo-conservation-and-management/quotas-harvest-data
NULL

#' Southern Bluefin Tuna in Australian Waters
#'
#' CSIRO Fisheries Biomass Provisioning Services
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item Catch: kilograms caught per year
#'   \item GVP: Gross value product in dollars as market price times catch
#'   \item EV: Exchange value in dollars as the value of ecosystem services as if a market existed
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_SouthernBluefinTuna
#' @format csv file with 22 rows and 4 columns
#' @source https://doi.org/10.25919/nzrp-6702
NULL

#' Tropical Rock Lobster in Australian Waters
#'
#' CSIRO Fisheries Biomass Provisioning Services
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item Catch: kilograms caught per year
#'   \item GVP: Gross value product in dollars as market price times catch
#'   \item EV: Exchange value in dollars as the value of ecosystem services as if a market existed
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_TropicalRockLobsters
#' @format csv file with 22 rows and 4 columns
#' @source https://doi.org/10.25919/nzrp-6702
NULL

#' Sydney Drinking Water Catchment
#'
#' WaterNSW WaterInsight, Water Storage from August 2015 to July 2025
#'
#' \itemize{
#'   \item Day: time variable for days from August 2015
#'   \item All: All Storage volume in gigalitres
#'   \item Blue Mtns: Blue Mountains Dams volume in gigalitres
#'   \item Nepean: Nepean Dam volume in gigalitres
#'   \item Avon: Avon Dam volume in gigalitres
#'   \item Wingecarribe: Wingecarribe Reservoir volume in gigalitres
#'   \item Cordeaux: Cordeaux Dam volume in gigalitres
#'   \item Cataract: Cataract Dam volume in gigalitres
#'   \item Warragamba: Warragamba Dam volume in gigalitres
#'   \item Woronora: Woronora Dam volume in gigalitres
#'   \item Prospect: Prospect Reservoir volume in gigalitres
#'   \item Tallowa: Tallowa Dam volume in gigalitres
#'   \item Fitzroy: Fitzroy Falls Dam volume in gigalitres
#'   \item Year: time variable as decimal year
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_SydneyWater
#' @format csv file with 120 rows and 14 columns
#' @source https://waterinsights.waternsw.com.au/12964-sydney-drinking-water-catchment/storage
NULL

#' Metals and Energy Commodities
#'
#' Closing prices on the first day of each month from 1 January 2009 to 1 May 2025.
#'
#' \itemize{
#'   \item Day: time variable in days since 1 January 2009
#'   \item Gold: London Bullion Market fix at the Perth (Aus) Mint in USD/troy ounce
#'   \item Silver: London Bullion Market fix in USD/troy ounce
#'   \item Copper: London Metals Exchange spot price, Copper, Grade A, USD/tonne
#'   \item Iron Ore: Cleared for Export Tianjin port, Iron Ore Fines 62% FE spot price in USD/tonne
#'   \item WTI: West Texas Intermediate in USD/bbl
#'   \item Brent: Brent Crude in USD/bbl
#'   \item Year: time variable as decimal year
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Finance_Commodities
#' @format csv file with 197 rows and 8 columns
#' @source https://files.marketindex.com.au/files/workbooks/commodities-workbook.xlsx
NULL

#' Kansas City Hard Red Wheat Futures
#'
#' Daily futures for September 2025 maturity from 1 July 2024 to 30 June 2025
#'
#' \itemize{
#'   \item Day: time variable in days since 1 July 2024
#'   \item Open: Price at market open in USD/ton
#'   \item High: Daily high price in USD/ton
#'   \item Low: Daily low price in USD/ton
#'   \item Close: Price at market close in USD/ton
#'   \item Adj Close: Same as Close for futures
#'   \item Volume: Daily number of trades
#'   \item Year: time variable as decimal year
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Finance_KansasCity_WheatFutures
#' @format csv file with 250 rows and 8 columns
#' @source https://finance.yahoo.com/quote (search for KE=F)
NULL

#' US Dollars per Australian Dollar
#'
#' Daily exchange rates from 1 July 2024 to 30 June 2025,
#'  with AUD as the base currency and USD as the quote currency
#'
#' \itemize{
#'   \item Day: time variable in days since 1 July 2024
#'   \item Open: Rate at session open in USD/AUD
#'   \item High: Daily high rate in USD/AUD
#'   \item Low: Daily low rate in USD/AUD
#'   \item Close: Rate at session close in USD/AUD
#'   \item Adj Close: Same as Close for exchange rates
#'   \item Year: time variable as decimal year
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Finance_USDAUD
#' @format csv file with 258 rows and 7 columns
#' @source https://finance.yahoo.com/quote (search for AUDUSD=X)
NULL
