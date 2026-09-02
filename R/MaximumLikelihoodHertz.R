library(R6)
library(plotly)
library(stringr)
library(clipr)

# roxygen ----
#' @title R6 class implementing Maximum Likelihood estimation of the Ornstein-Uhlenbeck Process.
#'
#' @description
#' Maximum Likelihood estimation using a modified Nelder-Mead algorithm with testing
#'  for simple hypotheses.
#'
#' @details # Formulas and Methods:
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
#'       alphau: identifies distribution of lnLu 0<=alphau<=1
#'       lnLr:   restricted log likelihood -inf<lnLr<=lnLu
#'       alphar: identifies distribution of lnLr 0<=alphar<=1
#'       m1:     number of observations 0<m1=m-1
#'
#' @details # Usage:
#' The MaximumLikelihood object must first be instantiated before its methods
#'  are called.  There are two ways.  The first way instantiates the OUProcess
#'  object and calls a function to get a pointer:
#'
#'       OUP <- OUProcess$new()
#'       A <- OUP$get_Analytical()
#'       FD <- OUP$get_FiniteDifference()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'
#' The MaximumLikelihood object will coordinate arguments to functions with the
#'  Analytical, FiniteDifference and MonteCarlo objects.  The second way
#'  instantiates the MaximumLikelihood object by itself with no coordination:
#'
#'       ML <- MaximumLikelihood$new()
#'
#' Once the object is instantiated, its methods can be called, to estimate the
#'  parameters of the Ornstein-Uhlenbeck Process, for example:
#'
#'       ML$Estimate()
#'
#' The plot methods can be used to customize the plots, with a title, for example:
#'
#'       ML$PlotEstimates(title="My Estimates")
#'
#' Other functions and methods are called in the same way.  To see all the
#'  possibilities and, in particular, how to read in data, check out the demos below.
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
#' The Nelder-Mead algorithm is used to maximize a Log Likelihood derived from
#'  the Transition Density of the Ornstein-Uhlenbeck Process. The original
#'  algorithm is designed to minimize and has been called the amoeba algorithm.
#'  If spread across a bumpy surface, the amoeba pulls itself over bumps and out
#'  of hollows and flows down to the lowest level.  There it contracts around
#'  its center.  If spread across a flat surface, the amoeba shrinks without
#'  flowing.  Once the amoeba contracts or shrinks, it is spread across the
#'  surface opposite from where it came.  Again it flows, contracts and shrinks.
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
#'  values are calculated from the data and can be updated automatically. The
#'  arguments for starting values should seldom be needed.
#'
#' Data are entered as a data frame with columns for tau, a vector of times, and
#'  z, a vector of states. Arguments taucol and zcol are columns to extract from
#'  a data frame containing many columns.
#'
#' To read a csv file into a data frame, create a MaximumLikelihood object,
#'  estimate parameters and test the goodness-of-fit:
#'
#'       df <- read.csv("data/MyData.csv")
#'       ML <- MaximumLikelihood$new()
#'       ML$Estimates(df)
#'       ML$GoodnessOfFit()
#'
#' In this example, the 'data' directory is in the working directory. The taucol
#'  and zcol arguments are optional and will default to 1 and 2. Continue on by
#'  restricting a parameter and performing a likelihood ratio test:
#'
#'       ML$Estimates(rhor=0.5)
#'       ML$LikelihoodRatioTest()
#'
#' You can also list and read from available data sets using:
#'
#'      OUPDataList()
#'      df <- OUPDataRead(27)
#'
#' Number 27 in the list is the data set "OUP_Convergence".  The available data
#'  sets are documented. Type one of the following categories:
#'
#'       ?Agric_
#'       ?Climate_
#'       ?Ecosys_
#'       ?Finance_
#'       ?OUP_
#'
#'  and then select from the drop-down list.  Or use the function:
#'
#'       OUPDataHelp(27)

# class ----
MaximumLikelihood <- R6::R6Class("MaximumLikelihood",
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
    #' Create a MaximumLikelihood object
    #' @param OUP pointer set by the OUProcess object
    #' @return A new MaximumLikelihood object
    initialize = function(OUP=NULL)
    {
      # pointer to container object ----
      if(!is.null(OUP) && class(OUP)[[1]] == "OUProcess") { private$OUP <- OUP }
      # arguments ----
      private$oup_params <- list(rho=0.459755421,mu=-22.8712176,sigma=13.7350886)
      private$oup_params_restr <- list(rhor=0.5,mur=-15,sigmar=15)
      private$oup_params_start <- list(rhos=0.5,mus=-15,sigmas=15)
      # time series ----
      tau <- c(20,20.05,20.1,20.15,20.2,20.25,20.3,20.35,20.55,20.6,20.65,20.7,20.75,20.8,20.95,21,21.05,21.1,21.15,21.2,21.25,21.3,21.35,21.4,21.45,21.5,21.55,21.6,21.65,21.7,21.75,21.8,21.85,21.9,21.95,22,22.05,22.1,22.15,22.2,22.25,22.3,22.35,22.4,22.45,22.5,22.55,22.6,22.65,22.7,22.75,22.8,22.85,22.9,22.95,23,23.05,23.45,23.5,23.55,23.6,23.65,23.7,23.75,23.8,23.85,23.9,23.95,24.1,24.15,24.2,24.25,24.3,24.35,24.4,24.45,24.5,24.6,24.65,24.75,24.8,24.95,25,25.05,25.1,25.15,25.2,25.25,25.3,25.35,25.4,25.45,25.5,25.55,25.6,25.95,26,26.05,26.1,26.15,26.2,26.25,26.3,26.35,26.4,26.45,26.5,26.55,26.6,26.65,26.7,26.75,26.85,26.9,26.95,27,27.05,27.1,27.15,27.2,27.25,27.3,27.35,27.4,27.45,27.5,27.55,27.6,27.65,27.7,27.75,27.8,27.85,27.9,27.95,28,28.05,28.1,28.15,28.2,28.25,28.3,28.35,28.4,28.45,28.5,28.55,28.6,28.65,28.7,28.75,28.8,28.85,28.9,28.95,29,29.05,29.1,29.15,29.2,29.25,29.3,29.35,29.4,29.45,29.5,29.55,29.6,29.65,29.7,29.75,29.8,29.85,29.9,29.95,30)
      z <- c(30,27.49,24.08,23.45,24.72,21.5,22.15,21.05,27.9,27.71,33.01,30.55,32.42,26.43,24.37,23.22,17.02,17.89,20.04,17.9,19.52,16.78,16.13,16.91,13.92,15.02,13.04,8.19,10.62,12.04,10.39,11.94,14.62,6.61,1.7,1.71,-0.22,-2.69,1.69,-2.98,0.19,3.36,-0.3,-1.24,-4.35,-5.14,-14.63,-14.25,-19.83,-20.51,-20.98,-17.14,-18.94,-21.53,-25.73,-23.74,-24.57,-20.15,-22.54,-20.45,-17.59,-19.29,-21.85,-23.17,-23.8,-27.92,-26.87,-24.97,-17.32,-13.92,-13.49,-11.5,-9.23,-10.3,-9.54,-7.25,-4.04,-9.64,-9.42,-7.87,-6.86,-6.68,-9.03,-20.15,-19.21,-20.75,-21.07,-17.12,-9.99,-12.91,-12.19,-11.42,-8.96,-13.32,-17.8,-21.72,-26.66,-29.16,-28.77,-28.06,-22.83,-24.47,-29.79,-31.28,-34.05,-36.76,-38.97,-39.09,-32.89,-32.89,-31.25,-32.98,-36.19,-36.12,-40.74,-40.57,-42.75,-40.7,-38.5,-38.91,-40.66,-45.96,-48.73,-49.76,-48.19,-47.9,-47.73,-43.19,-41.62,-41.96,-42.98,-44.98,-41.98,-40.65,-36.21,-33.61,-33.73,-29.53,-34.07,-37.05,-33.84,-29.71,-29.24,-30.83,-26.19,-28.73,-30.27,-32.39,-30.3,-25.83,-29.39,-27.64,-25.63,-23.81,-23.34,-18.32,-16.53,-13.82,-13.74,-13.28,-16.98,-13,-13.03,-14.26,-16.94,-15.62,-13.54,-19.63,-14.06,-17.23,-17.6,-19.67,-19.87,-20.65,-14.42,-10.23)
      private$timeseries <- data.frame(tau=tau,z=z)
      private$timeseries_info <- list(tbeg=20,tend=30,dataname="Default",timename="time",statename="state",estimation="Unrestricted")
      # plot info ----
      plotfont <- list(family="Cambria",size=14)
      plotfile <- list(format="png",width=640,height=480)
      plottheme <- list(name="light",opaque=0.0)
      if(Sys.getenv("RSTUDIO") == "1")
      {
        if(rstudioapi::getThemeInfo()$dark) { plottheme <- list(name="dark",opaque=1.0) }
      }
      private$plot_info <- list(plotfont=plotfont,plotfile=plotfile,plottheme=plottheme,plotlabels=TRUE)
      # flags ----
      private$flags <- list(plotit=FALSE,copyit=FALSE)
      # plot globals ----
      private$plot_colors <- private$rainbow(plottheme$name,plottheme$opaque)
      private$modebar_2D <- list(list("zoom2d","pan2d","resetScale2d"),list("toImage"))
      private$modebar_3D <- list(list("zoom3d","pan3d","orbitRotation","tableRotation"),list("resetCameraDefault3d","resetCameraLastSave3d"),list("hoverClosest3d"),list("toImage"))
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
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_oup_params(rho,mu,sigma,"ML") }
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
            private$theta <- NULL
            private$theta_l <- NULL
            private$goodness <- NULL
            private$likelyratio <- NULL
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
            private$theta <- NULL
            private$theta_l <- NULL
            private$goodness <- NULL
            private$likelyratio <- NULL
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
            private$theta <- NULL
            private$theta_l <- NULL
            private$goodness <- NULL
            private$likelyratio <- NULL
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
      if(!is.null(rhor) && is.null(rhosca)) { message("rhor not set.")}
      musca <- private$extract_scalar(mur)
      if(!is.null(mur) && is.null(musca)) { message("mur not set.")}
      sigmasca <- private$extract_scalar(sigmar)
      if(!is.null(sigmar) && is.null(sigmasca)) { message("sigmar not set.")}
      rhoparam <- private$oup_params_restr[[1]]
      muparam <- private$oup_params_restr[[2]]
      sigmaparam <- private$oup_params_restr[[3]]
      hit <- FALSE
      if(!is.null(rhosca) && is.null(rhoparam) || is.null(rhosca) && !is.null(rhoparam)) { hit <- TRUE }
      else if(!is.null(musca) && is.null(muparam) || is.null(musca) && !is.null(muparam)) { hit <- TRUE }
      else if(!is.null(sigmasca) && is.null(sigmaparam) || is.null(sigmasca) && !is.null(sigmaparam)) { hit <- TRUE }
      else if(!is.null(rhosca) && (rhosca <= 0 && rhoparam > 0 || rhosca > 0 && rhosca != rhoparam)) { hit <- TRUE }
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
              if(!is.null(tc) && !is.null(zc))
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
              else if(is.null(tc) && is.null(zc))
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
                dataname <- df$dfname[1]
                df <- data.frame(df[,c(tc,zc)])
                df <- cleandata(df)
                nrows <- nrow(df)
                if(nrows > 2)
                {
                  i <- 1
                  while(i < nrows && OK == TRUE)
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
                    private$timeseries_info <- list(tbeg=df[1,1],tend=df[nrows,1],dataname=dataname,timename=framenames[1],statename=framenames[2],estimation=NULL)
                    private$oup_params_restr <- list(rhor=NULL,mur=NULL,sigmar=NULL)
                    private$oup_params_start <- list(rhos=NULL,mus=NULL,sigmas=NULL)
                    private$theta <- NULL
                    private$theta_l <- NULL
                    private$theta_u <- NULL
                    private$theta_r <- NULL
                    private$goodness <- NULL
                    private$likelyratio <- NULL
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
        if(!is.null(taucol) || !is.null(zcol))
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
        if(tbeg == -Inf) { private$timeseries_info$tbeg <- tau[1] }
        else if(tbeg == Inf) { private$timeseries_info$tbeg <- tau[m] }
        else
        {
          sca <- private$extract_scalar(tbeg)
          if(!is.null(sca))
          {
            if(sca < tau[1]) { private$timeseries_info$tbeg <- tau[1] }
            else if(sca > tau[m]) { private$timeseries_info$tbeg <- tau[m] }
            else { private$timeseries_info$tbeg <- sca }
          }
          else { message("beg not set.") }
        }
      }
      if(!is.null(tend))
      {
        if(tend == Inf) { private$timeseries_info$tend <- tau[m] }
        else if(tend == -Inf) { private$timeseries_info$tend <- private$timeseries_info$tbeg }
        else
        {
          sca <- private$extract_scalar(tend)
          if(!is.null(sca))
          {
            if(sca > tau[m]) { private$timeseries_info$tend <- tau[m] }
            else if(sca < private$timeseries_info$tbeg) {  private$timeseries_info$tend <- private$timeseries_info$tbeg }
            else { private$timeseries_info$tend <- sca }
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
    #' @param who        object id of sender
    #' @return list(font,file,theme,labels)
    set_plot_info = function(fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,labels=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,NULL,NULL,labels,"ML") }
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
      if(!is.null(labels))
      {
        bool <- private$extract_boolean(labels)
        if(!is.null(bool)) { private$plot_info$plotlabels <- bool  }
        else { message("labels not set.") }
      }
      return(private$plot_info)
    },
    #' @description
    #' Set flags for plotting and copying
    #' @param plotit automatic plot after calculation TRUE or FALSE
    #' @param copyit copy to clipboard TRUE or FALSE
    #' @param who    object id of sender
    #' @return list(plotit,copyit)
    set_flags = function(plotit=NULL,copyit=NULL,who=NULL)
    {
      if(is.null(who) && !is.null(private$OUP)) { private$OUP$send_flags(plotit,copyit,"ML") }
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
    #' Get all arguments, time series and information
    #' @return list(oup_params,oup_params_restr,oup_params_start,timeseries,timeseries_info,plot_info)
    get_all = function()
    {
      all <- list(oup_params = private$oup_params,
        oup_params_restr = private$oup_params_restr,
        oup_params_start = private$oup_params_start,
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
    #' Get OUP parameter restrictions
    #' @return list(rhor,mur,sigmar)
    get_oup_params_restr = function() { return(private$oup_params_restr) },
    #' @description
    #' Get OUP parameter starting values
    #' @return list(rhos,mus,sigmas)
    get_oup_params_start = function() { return(private$oup_params_start) },
    #' @description
    #' Get time series data for time tau and state z
    #' @return dataframe(tau,z)
    get_timeseries = function() { return(private$timeseries) },
    #' @description
    #' Get information for times series
    #' @return list(tbeg,tend,datename,timename,statename,estimation)
    get_timeseries_info = function() { return(private$timeseries_info) },
    #' @description
    #' Get information for plotting
    #' @return list(font,file,theme,labels)
    get_plot_info = function() { return(private$plot_info) },
    #' @description
    #' Get colors for plotting
    #' @return list(red,ylw,grn,cyn,blu,mgn,gry,background,font)
    get_plot_colors = function() { return(private$plot_colors) },
    #' @description
    #' Get flags for plotting and copying
    #' @return list(plotit,copyit)
    get_flags = function() { return(private$flags) },
    # public calculate methods ----
    #' @description
    #' Calculate the log likelihood of estimates
    #' @param rho    rate parameter 0<=rho<inf
    #' @param mu     location parameter -inf<mu<inf
    #' @param sigma  scale parameter -inf<sigma<inf
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index of a column containing z
    #' @return list(rho,mu,sigma,lnL,k,alpha,m1)
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
      copyit <- private$flags[[2]]
      # calculate ----
      self$set_timeseries_info(NULL,NULL,NULL,NULL,NULL,"LogLikely")
      theta <- private$theta_l
      if(is.null(theta))
      {
        logL <- RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
        theta <- list(rho=rho,mu=mu,sigma=sigma,lnL=logL[1],k=logL[2],alpha=logL[3],m1=logL[4])
        private$theta_l <- theta
      }
      private$theta <- theta
      # copy ----
      if(copyit == TRUE)
      {
        dataname <- private$timeseries_info[[3]]
        timename <- private$timeseries_info[[4]]
        statename <- private$timeseries_info[[5]]
        estimation <- private$timeseries_info[[6]]
        clip <- rbind(c("Maximum Likelihood",""),c("Log Likelihood",""),c("File:",dataname),c("Time:",timename),c("State:",statename),c("rho",theta[1]),c("mu",theta[2]),c("sigma",theta[3]),c("lnL",theta[4]),c("k",theta[5]),c("alpha",theta[6]),c("m1",theta[7]))
        private$CopyToClipboard(clip)
      }
      return(theta)
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
    #' @return list(rhohat,muhat,sigmahat,lnLu,ku,alphau,m1) or list(rhor,mur,sigmar,lnLr,kr,alphar,m1)
    Estimates = function(df=NULL,taucol=NULL,zcol=NULL,rhor=NULL,mur=NULL,sigmar=NULL,rhos=NULL,mus=NULL,sigmas=NULL)
    {
      # set / get ----
      self$set_timeseries(df,taucol,zcol)
      self$set_oup_params_start(rhos,mus,sigmas)
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      rhos <- private$oup_params_start[[1]]
      mus <- private$oup_params_start[[2]]
      sigmas <- private$oup_params_start[[3]]
      if(is.null(rhor) && is.null(mur) && is.null(sigmar))
      {
        self$set_timeseries_info(NULL,NULL,NULL,NULL,NULL,"Unrestricted")
        theta <- private$theta_u
      }
      else
      {
        self$set_timeseries_info(NULL,NULL,NULL,NULL,NULL,"Restricted")
        self$set_oup_params_restr(rhor,mur,sigmar)
        rhor <- private$oup_params_restr[[1]]
        mur <- private$oup_params_restr[[2]]
        sigmar <- private$oup_params_restr[[3]]
        theta <- private$theta_r
      }
      plotit <- private$flags[[1]]
      copyit <- private$flags[[2]]
      # calculate ----
      if(is.null(theta))
      {
        nmstart <- RcppOUPMLNMStart(tau,z,rhor,mur,sigmar,rhos,mus,sigmas)
        nk <- nmstart[2,4]
        nm <- RcppOUPMLNelderMead(tau,z,nmstart[1,1:3],nmstart[2,1:3])
        rho <- nm[1]
        mu <- nm[2]
        sigma <- nm[3]
        logL <- RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
        lnL <- logL[1]
        alpha <- logL[3]
        m1 <-logL[4]
      # reset and set ----
        if(nmstart[2,1] != 0 && nmstart[2,2] != 0 && nmstart[2,3] != 0)
        {
          private$oup_params_start$rhos <- rho
          private$oup_params_start$mus <- mu
          private$oup_params_start$sigmas <- sigma
          theta <- list(rhohat=rho,muhat=mu,sigmahat=sigma,lnLu=lnL,ku=nk,alphau=alpha,m1=m1)
          private$theta_u <- theta
        }
        else
        {
          if(nmstart[2,1] != 0) { theta <- list(rhobar=rho) }
          else { theta <- list(rhor=rho) }
          if(nmstart[2,2] != 0) { theta <- append(theta,list(mubar=mu)) }
          else { theta <- append(theta,list(mur=mu)) }
          if(nmstart[2,3] != 0) { theta <- append(theta,list(sigmabar=sigma)) }
          else { theta <- append(theta,list(sigmar=sigma)) }
          theta <- append(theta,list(lnLr=lnL,kr=nk,alphar=alpha,m1=m1))
          private$theta_r <- theta
        }
      }
      self$set_oup_params(rho=theta[[1]],mu=theta[[2]],sigma=theta[[3]])
      private$theta <- theta
      # plot or copy ----
      if(plotit == TRUE) { print(self$PlotEstimates()) }
      else if(copyit == TRUE)
      {
        names <- names(theta)
        theta <- unname(theta)
        dataname <- private$timeseries_info[[3]]
        timename <- private$timeseries_info[[4]]
        statename <- private$timeseries_info[[5]]
        estimation <- private$timeseries_info[[6]]
        clip <- rbind(c("Maximum Likelihood",""),c("Estimates",""),c("File:",dataname),c("Time:",timename),c("State:",statename),c(estimation,""),c(names[1],theta[1]),c(names[2],theta[2]),c(names[3],theta[3]),c(names[4],theta[4]),c(names[5],theta[5]),c(names[6],theta[6]),c(names[7],theta[7]))
        private$CopyToClipboard(clip)
      }
      return(theta)
    },
    #' @description
    #' Compare maximum likeliood estimates with invariant and scaled brownian motion estimates
    #' @param rho    rate parameter 0<=rho<inf
    #' @param mu     location parameter -inf<mu<inf
    #' @param sigma  scale parameter -inf<sigma<inf
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index of a column containing z
    #' @return list(theta,theta_i,theta_s,Inv,SBM) with Inv=list(R2,Pval) and SBM=list(R2,Pval)
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
      copyit <- private$flags[[2]]
      # calculate ----
      goodness <- private$goodness
      if(is.null(goodness))
      {
        theta <- private$theta #last estimate (user, unrestricted or restricted)
        if(is.null(theta))
        {
          logL <- RcppOUPMLLogLikelihood(tau,z,rho,mu,sigma)
          theta <- list(rho=rho,mu=mu,sigma=sigma,lnL=logL[1],k=logL[2],alpha=logL[3],m1=logL[4])
          private$theta_l <- theta
          private$theta <- theta
        }
        gof <- RcppOUPMLGoodnessOfFit(tau,z,theta[[4]],theta[[6]])
        theta_i <- list(rhor=gof[1,1],mu=gof[1,2],sigma=gof[1,3],lnLr=gof[1,4],k=gof[1,5],alphar=gof[1,6],m1=gof[1,7])
        Inv <- list(R2=gof[1,8],PVal=gof[1,9])
        theta_s <- list(rhor=gof[2,1],mur=gof[2,2],sigma=gof[2,3],lnLr=gof[2,4],k=gof[2,5],alphar=gof[2,6],m1=gof[2,7])
        SBM <- list(R2=gof[2,8],Pval=gof[2,9])
        goodness <- list(theta=theta,theta_i=theta_i,theta_s=theta_s,Inv=Inv,SBM=SBM)
        private$goodness <- goodness
      }
      # copy ----
      if(copyit == TRUE)
      {
        theta <- goodness[[1]]
        theta_i <- goodness[[2]]
        theta_s <- goodness[[3]]
        Inv <- goodness[[4]]
        SBM <- goodness[[5]]
        dataname <- private$timeseries_info[[3]]
        timename <- private$timeseries_info[[4]]
        statename <- private$timeseries_info[[5]]
        estimation <- private$timeseries_info[[6]]
        clip <- rbind(c("Maximum Likelihood",rep("",3)),c("Goodness of Fit",rep("",3)),c("File:",dataname,rep("",2)),c("Time:",timename,rep("",2)),c("State:",statename,rep("",2)),c("Inv R2",Inv[[1]],rep("",2)),c("Inv 1-P",Inv[[2]],rep("",2)),c("SBM R2",SBM[[1]],rep("",2)),c("SBM 1-P",SBM[[2]],rep("",2)),c("",estimation,"Invariant","Scaled BM"),c("rho",theta[1],theta_i[1],theta_s[1]),c("mu",theta[2],theta_i[2],theta_s[2]),c("sigma",theta[3],theta_i[3],theta_s[3]),c("lnL",theta[4],theta_i[4],theta_s[4]),c("k",theta[5],theta_i[5],theta_s[5]),c("alpha",theta[6],theta_i[6],theta_s[6]),c("m1",theta[7],theta_i[7],theta_s[7]))
        private$CopyToClipboard(clip)
      }
      return(goodness)
    },
    #' @description
    #' Test for significant difference between restricted and unrestricted likelihoods
    #' @param rho    rate parameter 0<=rho<inf
    #' @param mu     location parameter -inf<mu<inf
    #' @param sigma  scale parameter -inf<sigma<inf
    #' @param df     data frame containing columns for tau and z
    #' @param taucol index of a column containing tau
    #' @param zcol   index of a column containing z
    #' @return list(theta_u,theta,R2,Pval)
    LikelihoodRatioTest = function(rho=NULL,mu=NULL,sigma=NULL,df=NULL,taucol=NULL,zcol=NULL)
    {
      # set / get ----
      oup_params <- self$set_oup_params(rho,mu,sigma)
      self$set_timeseries(df,taucol,zcol)
      copyit <- private$flags[[2]]
      # calculate ----
      likelyratio <- private$likelyratio
      estimation <- private$timeseries_info[[6]]
      if(is.null(likelyratio))
      {
        theta <- private$theta #last estimate (loglikely, unrestricted or restricted)
        if(is.null(theta))
        {
          theta <- self$LogLikelihood()
          estimation <- private$timeseries_info[[6]]
        }
        theta_u <- private$theta_u
        if(is.null(theta_u))
        {
          theta_u <- self$Estimates()
          private$theta <- theta
          private$timeseries_info[[6]] <- estimation
        }
        lnLu <- theta_u[[4]]
        alpha <- theta_u[[6]]
        m <- theta_u[[7]]+1
        lnLr <- theta[[4]]
        lrt <- RcppOUPMLLikelihoodRatioTest(lnLu,alpha,m,lnLr)
        likelyratio <- list(theta_u=theta_u,theta=theta,R2=lrt[1],Pval=lrt[2])
        private$likelyratio <- likelyratio
      }
      # copy ----
      if(copyit == TRUE)
      {
        R2 <- likelyratio[[3]]
        Pval <- likelyratio[[4]]
        dataname <- private$timeseries_info[[3]]
        timename <- private$timeseries_info[[4]]
        statename <- private$timeseries_info[[5]]
        estimation <- private$timeseries_info[[6]]
        clip <- rbind(c("Maximum Likelihood",rep("",2)),c("Likelihood Ratio Test",rep("",2)),c("File:",dataname,""),c("Time:",timename,""),c("State:",statename,""),c("R2",R2,""),c("1-P",Pval,""),c("","Unrestricted",estimation),c("rho",theta_u[1],theta[1]),c("mu",theta_u[2],theta[2]),c("sigma",theta_u[3],theta[3]),c("lnL",theta_u[4],theta[4]),c("k",theta_u[5],theta[5]),c("alpha",theta_u[6],theta[6]),c("m1",theta_u[7],theta[7]))
        private$CopyToClipboard(clip)
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
      tbeg <- private$timeseries_info[[1]]
      tend <- private$timeseries_info[[2]]
      dataname <- private$timeseries_info[[3]]
      timename <- private$timeseries_info[[4]]
      statename <- private$timeseries_info[[5]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      m <- length(tau)
      Inx <- index(tau,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        tau <- tau[Ixbeg:Ixend]
        z <- z[Ixbeg:Ixend]
      }
      # copy ----
      if(copyit == TRUE)
      {
        clip <- rbind(c("Maximum Likelihood",""),c("Time Series",""),c("File:",dataname),c(timename,statename),cbind(tau,z))
        private$CopyToClipboard(clip)
      }
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
        add_trace(.,type="scattergl",x=tau,y=z,mode="lines+markers",line=zline,marker=zmarker) %>%
        config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
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
      tau <- private$timeseries[[1]]
      z <- private$timeseries[[2]]
      rho <- private$oup_params[[1]]
      mu <- private$oup_params[[2]]
      sigma <- private$oup_params[[3]]
      tbeg <- private$timeseries_info[[1]]
      tend <- private$timeseries_info[[2]]
      dataname <- private$timeseries_info[[3]]
      timename <- private$timeseries_info[[4]]
      statename <- private$timeseries_info[[5]]
      estimation <- private$timeseries_info[[6]]
      font <- list(family=private$plot_info$plotfont$family,size=private$plot_info$plotfont$size,color=private$plot_colors$font)
      file <- private$plot_info$plotfile
      labels <- private$plot_info$plotlabels
      red <- private$plot_colors$red
      grn <- private$plot_colors$grn
      cyn <- private$plot_colors$cyn
      mgn <- private$plot_colors$mgn
      gry <- private$plot_colors$gry
      background <- private$plot_colors$background
      copyit <- private$flags[[2]]
      m <- length(tau)
      Inx <- index(tau,tbeg,tend)
      Ixbeg <- Inx[[1]]
      Ixend <- Inx[[2]]
      if(Ixbeg > 1 || Ixend < m)
      {
        tau <- tau[Ixbeg:Ixend]
        z <- z[Ixbeg:Ixend]
        m <- length(tau)
      }
      taumv <- tau[1:(m-1)]
      # calculate ----
      mean <- vector("double",m-1)
      variance <- vector("double",m-1)
      resid <- vector("double",m-1)
      if(rho < 0.0000000001)
      {
        for(i in 1:(m-1))
        {
          mean[i] <- z[i]
          variance[i] <- sigma^2*(tau[i+1]-tau[i])
          if(abs(sigma) < 0.0000000001) { resid[i] <- 0 }
          else { resid[i] <- (z[i+1]-mean[i])/variance[i]^0.5 }
        }
      }
      else
      {
        for(i in 1:(m-1))
        {
          mean[i] <- mu+(z[i]-mu)*exp(-rho*(tau[i+1]-tau[i]))
          variance[i] <- sigma^2/(2*rho)*(1-exp(-2*rho*(tau[i+1]-tau[i])))
          if(abs(sigma) < 0.0000000001) { resid[i] <- 0 }
          else { resid[i] <- (z[i+1]-mean[i])/variance[i]^0.5 }
        }
      }
      # copy ----
      if(copyit == TRUE)
      {
        theta <- private$theta
        if(is.null(theta)) { theta <- self$LogLikelihood() }
        names <- names(theta)
        theta <- unname(theta)
        clip <- rbind(c("Maximum Likelihood",""),c("Estimates",""),c("File:",dataname),c("Time:",timename),c("State:",statename),c(estimation,""),c(names[1],theta[1]),c(names[2],theta[2]),c(names[3],theta[3]),c(names[4],theta[4]),c(names[5],theta[5]),c(names[6],theta[6]),c(names[7],theta[7]))
        private$CopyToClipboard(clip)
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
        if(is.null(title)) { title <- paste(sep=":  ",dataname,estimation) }
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
      legendpos <- list(orientation="h",x=0.5,y=1.0,xanchor="center")
      fig <- plot_ly()  %>%
        add_trace(.,type="scattergl",x=tau,y=z,name="<i>z</i>",mode="markers",marker=zmarker) %>%
        add_trace(.,type="scattergl",x=taumv,y=mean,name="<i>G</i>",mode="lines+markers",line=meanline,marker=meanmarker) %>%
        add_trace(.,type="scattergl",x=taumv,y=variance,name="<i>H</i><sup>2</sup>",mode="lines+markers",line=varianceline,marker=variancemarker,visible="legendonly") %>%
        add_trace(.,type="scattergl",x=taumv,y=resid,name="<i>v</i>",mode="lines+markers",line=residline,marker=residmarker,visible="legendonly") %>%
        config(.,toImageButtonOptions=imageoptions,modeBarButtons=private$modebar_2D,displaylogo=FALSE) %>%
        layout(.,title=lookup,legend=legendpos,font=font,paper_bgcolor=background,plot_bgcolor=background,xaxis=horz,yaxis=vert,margin=list(t=50,r=40,b=100,l=40))

      return(fig)
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    OUP = NULL,
    # private attributes ----
    oup_params = NULL,
    oup_params_restr = NULL,
    oup_params_start = NULL,
    timeseries = NULL,
    timeseries_info = NULL,
    plot_info = NULL,
    flags = NULL,
    # private output fields ----
    theta = NULL,
    theta_l = NULL,
    theta_u = NULL,
    theta_r = NULL,
    goodness = NULL,
    likelyratio = NULL,
    # private globals ----
    modebar_2D = NULL,
    modebar_3D = NULL,
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
                if(is.finite(df[i,1]) && !is.na(df[i,1]) && is.finite(df[i,2]) && !is.na(df[i,2]))
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
    # private clipboard methods ----
    CopyToClipboard = function(clip)
    {
      if(!is.null(private$OUP)) { OUP$CopyToClipboard(clip) }
      else if(interactive() && clipr_available()) { write_clip(clip,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    }
  )
)
