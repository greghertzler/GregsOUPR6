library(R6)

# roxygen ----
#' R6 Class implementing OUProcess, a container to synchronize OUP objects.
#'
#' @description
#' The Analytical, FiniteDifference, MaximumLikelihood and MonteCarlo objects
#'  are a complete set of functions for maximum likelihood estimation and
#'  for the calculation of probabilities, option prices, visiting times, first
#'  passage times, decision thresholds and more--everything for a Real Options
#'  Analysis.
#' Each object can be used by itself.  This object, OUProcess, will instantiate
#'  and synchronize the other objects.
#'
#' @details # Functions are categorized by:
#'     Solution methods,
#'     Stochastic variables and
#'     Random variables.
#'
#' @details # Solution methods include:
#'     Analytical formulas,
#'     Finite Difference method,
#'     Maximum Likelihood estimation, and
#'     Monte Carlo simulation.
#'
#' @details # Stochastic variables can be:
#'     z, generic states of the system,
#'     y, future or forward states of the system,
#'     x, current or backward states of the system, and
#'     t, future or forward times.
#'
#' @details # Random variables can be::
#'     rho, mu and sigma, the OUP parameters.
#'
#' @details # Analytical formulas:
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
#'     t stochastic
#'       PassageTimeMode
#'       PassageTimeMedian
#'       PassageTimeMean
#'       PassageTimeVariance
#'       PassageTimePercentile
#'       PassageTimeDensity
#'       PassageTimeProbability
#'
#' @details # Finite Difference method:
#'     x stochastic
#'       Drift
#'       Diffusion
#'       TerminalValue_Linear
#'       TerminalValue_Kinked
#'       TerminalValue_Stepped
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
#' @details # Maximum Likelihood estimation:
#'     rho, mu and sigma random
#'       LogLikelihood
#'       Estimates
#'       GoodnessOfFit
#'       LikelihoodRatioTest
#'
#' @details # Monte Carlo simulation:
#'     y stochastic
#'       SamplePaths
#'       BoundedSamplePaths
#'       Probabilities
#'     x stochastic
#'       Option
#'     t stochastic
#'       PassageTimeProbabilities
#'
#' @details # Using the formulas and methods:
#' This object, OUProcess, and the Analytical, FiniteDifference,
#'  MaximumLikelihood and Monte-Carlo objects can be instantiated together
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
#'  Then methods of individual objects can be called by, for example:
#'
#'       A$Drift()
#'       FD$Drift()
#'       ML$LogLikelihood()
#'       MC$ForwardPaths()
#'
#' OUP will have pointers to A, FD, ML and MC.  A, FD, ML and MC will each
#'  have a pointer to OUP. Thus A, FD, ML and MC can follow the pointers to
#'  synchronize with each other. For example,rho, mu and sigma estimated in ML
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
#' The functions all return named lists and estimates is a list of rho, mu and
#'  sigma.  The functions will also accept named lists as arguments:
#'
#'       estimates <- ML$Estimates()
#'       rho <- estimates[1]
#'       mu <- estimates[2]
#'       sigma <- estimates[3]
#'       A$Option(rho=rho,mu=mu,sigma=sigma)
#'
#' The single brackets assign named lists to the left-hand variables. If a
#'  function is called with a named list as an argument, it strips the names
#'  before using the numbers in its calculations.  As a general rule, named lists
#'  as arguments do not work.  If you wish to use numbers for the arguments you
#'  should use double brackets:
#'
#'       estimates <- ML$Estimates()
#'       rho <- estimates$[[1]]
#'       mu <- estimates$[[2]]
#'       sigma <- estimates$[[3]]
#'       A$Option(rho,mu,sigma)
#'
#' This will strip the name and the variables on the left-hand side are just
#'  numbers.  The arguments do not have to be accompanied by their names, as in
#'  rho=rho, where the first rho is the name and the second rho is the variable.
#'  They can also be entered by their position, with rho first, mu second and
#'  sigma third.
#'
#' These examples will both print a table of option prices and draw a plot.
#'  The plot is optional:
#'
#'       A$Option(plotit=FALSE)
#'
#' The table of numbers is also optional:
#'
#'       options <- A$Option()
#'
#' or:
#'
#'       A$set_oup_params(rho=rho,mu=mu,sigma=sigma)
#'       A$PlotOption()
#'
#' Most methods have a corresponding plot routine. By design, the plot routines
#'  do not take parameters as arguments, which must be set separately.
#'
#' If rho, mu or sigma is changed in one method, it is changed for all methods.
#'  If a data frame is read, it is available to all methods. For example, a
#'  time series may be entered to estimate parameters and then used to calculate
#'  the goodness-of-fit:
#'
#'       df <- read.csv("data/OUP_SimulatedData.csv")
#'       ML$Estimates(df=df,taucol=1,zcol=3)
#'       ML$GoodnessOfFit()
#'
#' If the state variable is changed in the Analytical object, it may also
#'  change in the FiniteDifference object. Plot information is also shared across
#'  objects.  This complicated communication allows code to be split into modules.
#'
#' The communication within each module is also complicated.  The modules are
#'  reactive.  In other words, don't calculate anything twice.  There is a cache
#'  of inputs and a cache of outputs, managed by set/get functions.  In the
#'  diagram to follow, A are arguments, I are inputs, O are outputs, S is a set
#'  function and G is a get function.
#'
#'                     A
#'            cache    |    cache
#'       G <--  I  <-> S -->  O
#'       |             |
#'       I             I
#'
#' Arguments enter the set function which validates the arguments and compares
#'  the arguments to the input cache.  If there is a change, the input cache is
#'  set and all outputs in the output cache which depend on that argument are
#'  set to NULL.  Finally, the set function returns the new input cache.  The
#'  get function will return the input cache without checking or changes.
#'
#' Once the arguments are set, other functions can be called.  In the diagram
#'  below, C is a calculation function called without arguments.
#'
#'
#'             |    cache
#'       G --> C <->  O
#'             |
#'             O
#'
#' The calculation function gets the input cache and checks the output cache.
#'  If the output cache is NULL, outputs are calculated, the cache is set and
#'  then returned.  If the output cache is not NULL, it is returned.  Thus a
#'  calculation function calculates once and then acts like a get function.
#'
#' For convenience, arguments can be entered directly into the calculation
#'  functions.
#'
#'             A
#'             |    cache
#'       S <-> C <->  O
#'             |
#'             O
#'
#' A calculation function sends the arguments to the set function which does
#'  its work and returns the input cache.  Then the calculation function checks
#'  the output cache, does its work and returns the output cache.
#'
#' Most calculation functions have a corresponding plot function.  There is a
#'  complicated relationship between them.  In the R console, you probably
#'  want calculations.  In RShiny, you get the plots.  In the R console, the
#'  calculation functions are in charge and in RShiny, the plot functions
#'  are in charge.
#'
#' In the diagrams below, P is a plot function and T is a flag set to true and
#'  F is the same flag set to false.  By default, the flag is set to T.  There
#'  are three possibilities if the calculation functions are in charge.
#'
#'       F        T
#'       |        |              |
#'       C        C --> P        C --> P
#'       |        |     |        |     |
#'       O        O    plot      O    plot
#'
#' The flag is an argument to the calculation function.  If the calculation
#'  function is called with a false flag, it returns the output cache.  If
#'  the calculation function is called with a true flag, it returns the output
#'  cache and calls the plot function which returns a plot.
#'
#' To prevent a circular reference, the flag is always false if a plot function
#'  is in charge.
#'
#'               |              |
#'       C <-F-> P      C <-F-> P <-F-> C
#'               |              |
#'              plot           plots
#'
#' The plot function calls the calculation function with a false flag, gets the
#'  output cache and returns a plot.  Some plot functions manage more than one
#'  calculation function and return several different types of plots.
#'
#' In summary, methods include set, get, calculation and plot functions.  Help
#'  for a method is in the help for its object:
#'
#'       ?Analytical
#'       ?FiniteDifference
#'       ?MaximumLikelihood
#'       ?MonteCarlo
#'
#' These and a more complete html help can be accessed by:
#'
#'       OUPHelp()
#'
#' Scripts to demonstrate the methods are in files in the 'demo' directory. To
#'  see a list of the demos, type:
#'
#'       demo()
#'
#' To run the 'A_Density' demo, type:
#'
#'       demo(A_Density)
#'
#' or type:
#'
#'       source('demo/A_Density.R',echo=TRUE)
#'
#' or, open the A_Density.R file in RStudio and Source it from the menu.
#'
#' Available data sets for Maximum Likelihood estimation have their own help.
#'  Type one of the categories:
#'
#'       ?Agric_
#'       ?Climate_
#'       ?Ecosys_
#'       ?Finance_
#'       ?OUP_
#'
#'  and then select from the drop-down list.
#'
#' Examples in R6 don't work the same as in S3 and S4 modules.  There is only
#'  one example for an R6 object, not one for each method in the object.
#'  To run examples, the devtools::run_examples() works, but the R command
#'  example("OUProcess") doesn't.  You can copy commands to the clipboard,
#'  paste into the console and press Enter. Examples in this help and a simple
#'  example at the bottom can be run in this way.
#'
#'  But really, the examples are included to pass CRAN tests.  Start with
#'   RShiny.  Then if you decide to work in the console, run the demos.

# class ----
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
    #' @examples
    #'   OUP <- OUProcess$new()
    #'   A <- OUP$get_Analytical()
    #'   FD <- OUP$get_FiniteDifference()
    #'   ML <- OUP$get_MaximumLikelihood()
    #'   MC <- OUP$get_MonteCarlo()
    #'   ML$Estimates()
    #'   MC$get_oup_params()
    #'   A$Option()
    #'   FD$TerminalValue_Mitscherlich(plotit=FALSE)
    #'   FD$Option()
    initialize = function()
    {
      # instantiate objects
      private$A <- Analytical$new(self)
      private$FD <- FiniteDifference$new(self)
      private$ML <- MaximumLikelihood$new(self)
      private$MC <- MonteCarlo$new(self)
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
    #' Send OUP parameters
    #' @param rho   rate parameter 0<=rho<inf
    #' @param mu    location parameter -inf<mu<inf
    #' @param sigma scale parameter -inf<sigma<inf
    #' @param who   identifier for object sending the parameters
    send_oup_params = function(rho=NULL,mu=NULL,sigma=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$ML$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$MC$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
        }
        else if(who == "FD")
        {
          private$A$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$ML$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$MC$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
        }
        else if(who == "ML")
        {
          private$A$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$FD$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$MC$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
        }
        else if(who == "MC")
        {
          private$A$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$FD$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
          private$ML$set_oup_params(rho=rho,mu=mu,sigma=sigma,who=who)
        }
      }
    },
    #' @description
    #' send x as a stochastic state and its arguments
    #' @param s   vector of m times -inf<s<t
    #' @param x   vector of n states -inf<x<inf
    #' @param r   discount rate -inf<r<inf
    #' @param phi <=0 for exit option, >0 for entry option
    #' @param who identifier for object sending the parameters
    send_x_stoch_args = function(s=NULL,x=NULL,r=NULL,phi=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A") { private$FD$set_x_stoch_args(s=s,x=x,r=r,phi=phi,who=who) }
        else if(who == "FD") { private$A$set_x_stoch_args(s=s,x=x,r=r,phi=phi,who=who) }
      }
    },
    #' @description
    #' Send information for plotting
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
    send_plot_info = function(type=NULL,pmax=NULL,ptmax=NULL,fontfamily=NULL,fontsize=NULL,fileformat=NULL,filewidth=NULL,fileheight=NULL,theme=NULL,opaque=NULL,walls=NULL,floor=NULL,labels=NULL,who=NULL)
    {
      if(is.character(who))
      {
        if(who == "A")
        {
          private$FD$set_plot_info(type,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
          private$MC$set_plot_info(type,pmax,ptmax,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
        }
        else if(who == "FD")
        {
          private$A$set_plot_info(type,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
          private$MC$set_plot_info(type,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
        }
        else if(who == "ML")
        {
          private$A$set_plot_info(NULL,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,NULL,NULL,labels,who)
          private$FD$set_plot_info(NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,NULL,NULL,labels,who)
          private$MC$set_plot_info(NULL,NULL,NULL,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,NULL,NULL,labels,who)
        }
        else if(who == "MC")
        {
          private$A$set_plot_info(type,pmax,ptmax,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$FD$set_plot_info(type,fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,walls,floor,labels,who)
          private$ML$set_plot_info(fontfamily,fontsize,fileformat,filewidth,fileheight,theme,opaque,labels,who)
        }
      }
    }
  ),
  # private members ----
  private = list(
    # private pointers ----
    A = NULL,
    FD = NULL,
    ML = NULL,
    MC = NULL
  )
  # class end ----
)
# help launch ----
#' Function to launch roxygen help for GregsOUPR6
#'
#' @description
#' Launch ?OUProcess, ?Analytical, ?MaximumLikelihood and ?MonteCarlo in Help Viewer
#' @return Help instance
#' @export
#' @examples
#'   OUPHelp()
OUPHelp = function()
{
  utils::help("MonteCarlo","GregsOUPR6")
  utils::help("MaximumLikelihood","GregsOUPR6")
  utils::help("FiniteDifference","GregsOUPR6")
  utils::help("Analytical","GregsOUPR6")
  utils::help("OUProcess","GregsOUPR6")
}
#' Function to launch Ribbon Help for GregsOUPR6
#'
#' @description
#' Launch OUP_Help in the default browser.
#' @return Browser instance
#' @export
#' @examples
#'   OUPRibbonHelp()
OUPRibbonHelp = function() { utils::browseURL(system.file("ribbonhelp/OUP_Help.html",package="GregsOUPR6")) }

