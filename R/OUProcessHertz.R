library(R6)
library(shiny)
library(clipr)

# roxygen ----
#' @title R6 class implementing a container to synchronize OUP objects.
#'
#' @description
#' The Analytical, FiniteDifference, MaximumLikelihood and MonteCarlo objects
#'  are a complete set of functions for maximum likelihood estimation and
#'  for the calculation of probabilities, option prices, visiting times, first
#'  passage times, decision thresholds and more--everything for a Real Options
#'  Analysis.  Each object can be used by itself.  This object, OUProcess, will
#'  instantiate and synchronize the other objects.
#'
#' @details # Usage:
#' This object, OUProcess, and the Analytical, FiniteDifference,
#'  MaximumLikelihood and MonteCarlo objects can be instantiated together
#'  by:
#'
#'       OUP <- OUProcess$new()
#'
#' Then pointers to individual objects can be accessed by:
#'
#'       A <- OUP$get_Analytical()
#'       FD <- OUP$get_FiniteDifference()
#'       ML <- OUP$get_MaximumLikelihood()
#'       MC <- OUP$get_MonteCarlo()
#'
#' OUP has other public methods, but these are called by individual objects
#'  and won't do anything if called by a user.
#'
#' @details # Discussion
#'
#' OUP will have pointers to A, FD, ML and MC.  A, FD, ML and MC will each
#'  have a pointer to OUP. Thus A, FD, ML and MC can follow the pointers to
#'  synchronize with each other. For example, rho, mu and sigma estimated in ML
#'  will be propagated to the other objects:
#'
#'       ML$Estimates()
#'       A$Option()
#'
#' Alternately, each object can be instantiated individually:
#'
#'       A <- Analytical$new()
#'       FD <- FiniteDifference$new()
#'       ML <- MaximumLikelihood$new()
#'       MC <- MonteCarlo$new()
#'
#' In this case, A, FD, ML and MC will not have pointers and will not
#'  synchronize. To estimate rho, mu and sigma and calculate option prices
#'  would require:
#'
#'       estimates <- ML$Estimates()
#'       rho <- estimates$rho
#'       mu <- estimates$mu
#'       sigma <- estimates$sigma
#'       A$Option(rho=rho,mu=mu,sigma=sigma)
#'
#' Upon initialization, OUP will ask if it is being instantiated by an RShiny
#'  app.  If so, it will store the session pointer and use it to coordinate
#'  clipboard access between RShiny and the A, FD, ML and MC objects.

# class ----
#' @importFrom clipr clipr_available write_clip
#' @export
OUProcess <- R6::R6Class("OUProcess",
  portable = FALSE,
  cloneable = FALSE,
  #public members ----
  public = list(
    # constructor ----
    #' @description
    #' Create an OUProcess object as a container for other objects
    #' @return A new OUProcess object
    #' @param session session pointer to RShiny app
    initialize = function(session=NULL)
    {
      # instantiate objects
      private$A <- Analytical$new(self)
      private$FD <- FiniteDifference$new(self)
      private$ML <- MaximumLikelihood$new(self)
      private$MC <- MonteCarlo$new(self)
      if(!is.null(session) && requireNamespace("shiny",quietly=TRUE) && shiny::isRunning()) { private$session <- session }
    },
    # public pointers ----
    #' @description
    #' Pointer to Analytical object
    #' @return pointer
    get_Analytical = function() { return(private$A) },
    #' @description
    #' Pointer to FiniteDifference object
    #' @return pointer
    get_FiniteDifference = function() { return(private$FD) },
    #' @description
    #' Pointer to MaximumLikelihood object
    #' @return pointer
    get_MaximumLikelihood = function() { return(private$ML) },
    #' @description
    #' Pointer to MonteCarlo object
    #' @return pointer
    get_MonteCarlo = function() { return(private$MC) },
    # public send ----
    #' @description
    #' Function called by modules to coordinate parameters
    #' @param rho   rate parameter 0<=rho<inf
    #' @param mu    location parameter -inf<mu<inf
    #' @param sigma scale parameter -inf<sigma<inf
    #' @param who   object id of sender
    send_oup_params = function(rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_oup_params(rho,mu,sigma,who)
          private$ML$set_oup_params(rho,mu,sigma,who)
          private$MC$set_oup_params(rho,mu,sigma,who)
        }
        else if(who == "FD")
        {
          private$A$set_oup_params(rho,mu,sigma,who)
          private$ML$set_oup_params(rho,mu,sigma,who)
          private$MC$set_oup_params(rho,mu,sigma,who)
        }
        else if(who == "ML")
        {
          private$A$set_oup_params(rho,mu,sigma,who)
          private$FD$set_oup_params(rho,mu,sigma,who)
          private$MC$set_oup_params(rho,mu,sigma,who)
        }
        else if(who == "MC")
        {
          private$A$set_oup_params(rho,mu,sigma,who)
          private$FD$set_oup_params(rho,mu,sigma,who)
          private$ML$set_oup_params(rho,mu,sigma,who)
        }
      }
    },
    #' @description
    #' Function called by modules to coordinate arguments for y as a stochastic state
    #' @param t   vector of m times -inf<t<inf
    #' @param y   vector of n states -inf<y<inf
    #' @param x   initial state -inf<x<inf
    #' @param psi <=0 for integral -inf to y, >0 for integral y to inf
    #' @param who object id of sender
    send_y_stoch_args = function(t=NULL,y=NULL,x=NULL,psi=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A") { private$MC$set_y_stoch_args(t,y,x,psi,who) }
        else if(who == "MC") { private$A$set_y_stoch_args(t,y,NULL,x,psi,NULL,who) }
      }
    },
    #' @description
    #' Function called by modules to coordinate x as a stochastic state and its arguments
    #' @param s   vector of m times -inf<s<t
    #' @param x   vector of n states -inf<x<inf
    #' @param y   terminal state -inf<y<inf
    #' @param r   discount rate -inf<r<inf
    #' @param phi <=0 for exit option, >0 for entry option
    #' @param who object id of sender
    send_x_stoch_args = function(s=NULL,x=NULL,y=NULL,r=NULL,phi=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_x_stoch_args(s,x,NULL,r,phi,NULL,NULL,who)
          private$MC$set_x_stoch_args(s,x,y,r,phi,who)
        }
        else if(who == "FD")
        {
          private$A$set_x_stoch_args(s,x,NULL,y,r,phi,NULL,NULL,who)
          private$MC$set_x_stoch_args(s,x,y,r,phi,who)
        }
        else if(who == "MC")
        {
          private$A$set_x_stoch_args(s,x,NULL,y,r,phi,NULL,NULL,who)
          private$FD$set_x_stoch_args(s,x,NULL,r,phi,NULL,NULL,who)
        }
      }
    },
    #' @description
    #' Function called by modules to coordinate t stochastic arguments
    #' @param t     vector of m times -inf<t<inf
    #' @param k     threshold -inf<k<inf
    #' @param x     initial state -inf<x<inf
    #' @param omega degree of irreversibility 0<=omega<=1
    #' @param Ppct  passage time probability for a percentile  0.01<=Ppct<=0.99
    #' @param who   object id of sender
    send_t_stoch_args = function(t=NULL,k=NULL,x=NULL,omega=NULL,Ppct=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A") { private$MC$set_t_stoch_args(t,k,x,omega,Ppct,who) }
        else if(who == "MC") { private$A$set_t_stoch_args(t,k,NULL,x,NULL,omega,Ppct,who) }
      }
    },
    #' @description
    #' Function called by modules to coordinate arguments for plotting
    #' @param pmax  maximum transition density
    #' @param ptmax maximum visiting time and first passage time densities
    #' @param who   object id of sender
    send_plot_args = function(pmax=NULL,ptmax=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A") { private$MC$set_plot_args(pmax,ptmax,NULL,NULL,NULL,NULL,who) }
        else if(who == "MC") { private$A$set_plot_args(pmax,ptmax,who) }
      }
    },
    #' @description
    #' Function called by modules to coordinate information for plotting
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
    send_plot_info = function(fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
          private$MC$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
        }
        else if(who == "FD")
        {
          private$A$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
          private$MC$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
        }
        else if(who == "ML")
        {
          private$A$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$FD$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$MC$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
        }
        else if(who == "MC")
        {
          private$A$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$FD$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
        }
      }
    },
    #' @description
    #' Function called by modules to coordinate flags for plotting and copying
    #' @param plotit automatic plot after calculation TRUE or FALSE
    #' @param copyit copy to clipboard TRUE or FALSE
    #' @param who    object id of sender
    send_flags = function(plotit=NULL,copyit=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_flags(plotit,copyit,who)
          private$ML$set_flags(plotit,copyit,who)
          private$MC$set_flags(plotit,copyit,who)
        }
        else if(who == "FD")
        {
          private$A$set_flags(plotit,copyit,who)
          private$ML$set_flags(plotit,copyit,who)
          private$MC$set_flags(plotit,copyit,who)
        }
        else if(who == "ML")
        {
          private$A$set_flags(plotit,copyit,who)
          private$FD$set_flags(plotit,copyit,who)
          private$MC$set_flags(plotit,copyit,who)
        }
        else if(who == "MC")
        {
          private$A$set_flags(plotit,copyit,who)
          private$FD$set_flags(plotit,copyit,who)
          private$ML$set_flags(plotit,copyit,who)
        }
      }
    },
    # public clipboard methods ----
    #' @description
    #' Write to clipboard
    #' @param clip data frame, matrix or vector
    CopyToClipboard = function(clip)
    {
      if(!is.null(private$session))
      {
        if(is.data.frame(clip) || is.matrix(clip))
        {
          connect <- textConnection("str","w",local=TRUE)
          write.table(clip,file=connect,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
          close(connect)
        }
        else { str <- as.character(clip) }
        txt <- paste(str,collapse="\n")
        private$session$sendCustomMessage("sendToShiny",txt)
      }
      else if(interactive() && clipr_available()) { write_clip(clip,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    A = NULL,
    FD = NULL,
    ML = NULL,
    MC = NULL,
    session = NULL
  )
)
