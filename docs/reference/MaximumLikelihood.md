# R6 class implementing Maximum Likelihood estimation of the Ornstein-Uhlenbeck Process.

Maximum Likelihood estimation using a modified Nelder-Mead algorithm
with testing for simple hypotheses.

## Formulas and Methods:

    rho, mu and sigma random
      LogLikelihood
      Estimates
      GoodnessOfFit
      LikelihoodRatioTest

## Plots:

      PlotTimeSeries
      PlotEstimates

## Arguments of functions:

      All arguments are optional in all functions.
      rho:    rate parameter 0<=rho<inf
      mu:     location parameter -inf<mu<inf
      sigma:  scale parameter -inf<sigma<inf
      tau:    vector of observed times -inf<tau<inf
      z:      vector of observed states -inf<z<inf
      df:     data frame containing columns for tau and z
      taucol: index of a column containing tau
      zcol:   index of a column containing z
      rhor:   constant to fix the rate parameter 0<=rhor<inf
      mur:    constant to fix the location parameter -inf<mur<inf
      sigmar: constant to fix scale parameter -inf<sigmar<inf
      rhos:   starting value for the rate parameter 0<=rhos<inf
      mus:    starting value for the location parameter -inf<mus<inf
      sigmas: starting value the scale parameter -inf<sigmas<inf
      lnLu:   unrestricted log likelihood -inf<lnLu<=0
      alphau: identifies distribution of lnLu 0<=alphau<=1
      lnLr:   restricted log likelihood -inf<lnLr<=lnLu
      alphar: identifies distribution of lnLr 0<=alphar<=1
      m1:     number of observations 0<m1=m-1

## Usage:

The MaximumLikelihood object must first be instantiated before its
methods are called. There are two ways. The first way instantiates the
OUProcess object and calls a function to get a pointer:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      FD <- OUP$get_FiniteDifference()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()

The MaximumLikelihood object will coordinate arguments to functions with
the Analytical, FiniteDifference and MonteCarlo objects. The second way
instantiates the MaximumLikelihood object by itself with no
coordination:

      ML <- MaximumLikelihood$new()

Once the object is instantiated, its methods can be called, to estimate
the parameters of the Ornstein-Uhlenbeck Process, for example:

      ML$Estimate()

The plot methods can be used to customize the plots, with a title, for
example:

      ML$PlotEstimates(title="My Estimates")

Other functions and methods are called in the same way. To see all the
possibilities and, in particular, how to read in data, check out the
demos below.

## Demos:

Demonstration scripts are in files in the 'demo' directory. The most
convenient way to list and run demos are the commands:

      OUPDemoList()
      OUPDemoRun(<number of demo in list>)

Entering the demos by number in the list saves typing.

## Discussion:

The Nelder-Mead algorithm is used to maximize a Log Likelihood derived
from the Transition Density of the Ornstein-Uhlenbeck Process. The
original algorithm is designed to minimize and has been called the
amoeba algorithm. If spread across a bumpy surface, the amoeba pulls
itself over bumps and out of hollows and flows down to the lowest level.
There it contracts around its center. If spread across a flat surface,
the amoeba shrinks without flowing. Once the amoeba contracts or
shrinks, it is spread across the surface opposite from where it came.
Again it flows, contracts and shrinks. Spreading, flowing, contracting
and shrinking continue until the amoeba contracts around the same point
or shrinks at the same level at least twice.

In this implementation, the Nelder-Mead algorithm has been modified in
five ways.

1.  It is modified to maximize;

2.  It can set parameters rho, mu and/or sigma to constants for simple
    hypothesis tests;

3.  It has a tie-breaking condition to prevent cycling or freezing;

4.  It can accelerate over long distances for searching flat
    likelihoods;

5.  It checks for log likelihoods greater than zero to prevent errors in
    case the data is faked or artificially smoothed.

Once the unrestricted estimates are found, the arguments for restricting
the parameters can be set for simple hypothesis tests. By default,
starting values are calculated from the data and can be updated
automatically. The arguments for starting values should seldom be
needed.

Data are entered as a data frame with columns for tau, a vector of
times, and z, a vector of states. Arguments taucol and zcol are columns
to extract from a data frame containing many columns.

To read a csv file into a data frame, create a MaximumLikelihood object,
estimate parameters and test the goodness-of-fit:

      df <- read.csv("data/MyData.csv")
      ML <- MaximumLikelihood$new()
      ML$Estimates(df)
      ML$GoodnessOfFit()

In this example, the 'data' directory is in the working directory. The
taucol and zcol arguments are optional and will default to 1 and 2.
Continue on by restricting a parameter and performing a likelihood ratio
test:

      ML$Estimates(rhor=0.5)
      ML$LikelihoodRatioTest()

You can also list and read from available data sets using:

     OUPDataList()
     df <- OUPDataRead(27)

Number 27 in the list is the data set "OUP_Convergence". The available
data sets are documented. Type one of the following categories:

      ?Agric_
      ?Climate_
      ?Ecosys_
      ?Finance_
      ?OUP_

and then select from the drop-down list. Or use the function:

      OUPDataHelp(27)

## Methods

### Public methods

- [`MaximumLikelihood$new()`](#method-MaximumLikelihood-initialize)

- [`MaximumLikelihood$set_oup_params()`](#method-MaximumLikelihood-set_oup_params)

- [`MaximumLikelihood$set_oup_params_restr()`](#method-MaximumLikelihood-set_oup_params_restr)

- [`MaximumLikelihood$set_oup_params_start()`](#method-MaximumLikelihood-set_oup_params_start)

- [`MaximumLikelihood$set_timeseries()`](#method-MaximumLikelihood-set_timeseries)

- [`MaximumLikelihood$set_timeseries_info()`](#method-MaximumLikelihood-set_timeseries_info)

- [`MaximumLikelihood$set_plot_info()`](#method-MaximumLikelihood-set_plot_info)

- [`MaximumLikelihood$set_flags()`](#method-MaximumLikelihood-set_flags)

- [`MaximumLikelihood$get_all()`](#method-MaximumLikelihood-get_all)

- [`MaximumLikelihood$get_oup_params()`](#method-MaximumLikelihood-get_oup_params)

- [`MaximumLikelihood$get_oup_params_restr()`](#method-MaximumLikelihood-get_oup_params_restr)

- [`MaximumLikelihood$get_oup_params_start()`](#method-MaximumLikelihood-get_oup_params_start)

- [`MaximumLikelihood$get_timeseries()`](#method-MaximumLikelihood-get_timeseries)

- [`MaximumLikelihood$get_timeseries_info()`](#method-MaximumLikelihood-get_timeseries_info)

- [`MaximumLikelihood$get_plot_info()`](#method-MaximumLikelihood-get_plot_info)

- [`MaximumLikelihood$get_plot_colors()`](#method-MaximumLikelihood-get_plot_colors)

- [`MaximumLikelihood$get_flags()`](#method-MaximumLikelihood-get_flags)

- [`MaximumLikelihood$LogLikelihood()`](#method-MaximumLikelihood-LogLikelihood)

- [`MaximumLikelihood$Estimates()`](#method-MaximumLikelihood-Estimates)

- [`MaximumLikelihood$GoodnessOfFit()`](#method-MaximumLikelihood-GoodnessOfFit)

- [`MaximumLikelihood$LikelihoodRatioTest()`](#method-MaximumLikelihood-LikelihoodRatioTest)

- [`MaximumLikelihood$PlotTimeSeries()`](#method-MaximumLikelihood-PlotTimeSeries)

- [`MaximumLikelihood$PlotEstimates()`](#method-MaximumLikelihood-PlotEstimates)

------------------------------------------------------------------------

### `MaximumLikelihood$new()`

Create a MaximumLikelihood object

#### Usage

    MaximumLikelihood$new(OUP = NULL)

#### Arguments

- `OUP`:

  pointer set by the OUProcess object

#### Returns

A new MaximumLikelihood object

------------------------------------------------------------------------

### `MaximumLikelihood$set_oup_params()`

Set OUP parameters

#### Usage

    MaximumLikelihood$set_oup_params(
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of sender

#### Returns

list(rho,mu,sigma) Set OUP parameter restrictions

------------------------------------------------------------------------

### `MaximumLikelihood$set_oup_params_restr()`

#### Usage

    MaximumLikelihood$set_oup_params_restr(rhor = NULL, mur = NULL, sigmar = NULL)

#### Arguments

- `rhor`:

  rate parameter 0\<=rhor\<inf

- `mur`:

  location parameter -inf\<mur\<inf

- `sigmar`:

  scale parameter -inf\<sigmar\<inf

#### Returns

list(rhor,mur,sigmar) Set OUP parameter starting values

------------------------------------------------------------------------

### `MaximumLikelihood$set_oup_params_start()`

#### Usage

    MaximumLikelihood$set_oup_params_start(rhos = NULL, mus = NULL, sigmas = NULL)

#### Arguments

- `rhos`:

  rate parameter 0\<=rhos\<inf

- `mus`:

  location parameter -inf\<mus\<inf

- `sigmas`:

  scale parameter -inf\<sigmas\<inf

#### Returns

list(rhos,mus,sigmas)

------------------------------------------------------------------------

### `MaximumLikelihood$set_timeseries()`

Set time series data for time tau and state z

#### Usage

    MaximumLikelihood$set_timeseries(df = NULL, taucol = NULL, zcol = NULL)

#### Arguments

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index a column containing z

#### Returns

dataframe(tau,z)

------------------------------------------------------------------------

### `MaximumLikelihood$set_timeseries_info()`

Set information for plotting times series and estimates

#### Usage

    MaximumLikelihood$set_timeseries_info(
      tbeg = NULL,
      tend = NULL,
      dataname = NULL,
      timename = NULL,
      statename = NULL,
      estimation = NULL
    )

#### Arguments

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `dataname`:

  name for the data

- `timename`:

  name for time

- `statename`:

  name for state

- `estimation`:

  description of estimation

#### Returns

list(tbeg,Ixbeg,tend,Ixend,datename,timename,statename,estimation)

------------------------------------------------------------------------

### `MaximumLikelihood$set_plot_info()`

Set information for plotting

#### Usage

    MaximumLikelihood$set_plot_info(
      fontfamily = NULL,
      fontsize = NULL,
      fileformat = NULL,
      filewidth = NULL,
      fileheight = NULL,
      theme = NULL,
      opaque = NULL,
      labels = NULL,
      who = NULL
    )

#### Arguments

- `fontfamily`:

  font family for plot labels

- `fontsize`:

  font size for plot labels

- `fileformat`:

  'png' or 'svg'

- `filewidth`:

  pixel width of 2D plot, pixel width and height of 3D plot

- `fileheight`:

  pixel height of 2D plot

- `theme`:

  'light' or 'dark'

- `opaque`:

  transparent to opaque background 0.0\<=opaque\<=1.0

- `labels`:

  title and parameters TRUE or FALSE

- `who`:

  object id of sender

#### Returns

list(font,file,theme,labels)

------------------------------------------------------------------------

### `MaximumLikelihood$set_flags()`

Set flags for plotting and copying

#### Usage

    MaximumLikelihood$set_flags(plotit = NULL, copyit = NULL, who = NULL)

#### Arguments

- `plotit`:

  automatic plot after calculation TRUE or FALSE

- `copyit`:

  copy to clipboard TRUE or FALSE

- `who`:

  object id of sender

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `MaximumLikelihood$get_all()`

Get all arguments, time series and information

#### Usage

    MaximumLikelihood$get_all()

#### Returns

list(oup_params,oup_params_restr,oup_params_start,timeseries,timeseries_info,plot_info)

------------------------------------------------------------------------

### `MaximumLikelihood$get_oup_params()`

Get OUP parameters

#### Usage

    MaximumLikelihood$get_oup_params()

#### Returns

list(rho,mu,sigma)

------------------------------------------------------------------------

### `MaximumLikelihood$get_oup_params_restr()`

Get OUP parameter restrictions

#### Usage

    MaximumLikelihood$get_oup_params_restr()

#### Returns

list(rhor,mur,sigmar)

------------------------------------------------------------------------

### `MaximumLikelihood$get_oup_params_start()`

Get OUP parameter starting values

#### Usage

    MaximumLikelihood$get_oup_params_start()

#### Returns

list(rhos,mus,sigmas)

------------------------------------------------------------------------

### `MaximumLikelihood$get_timeseries()`

Get time series data for time tau and state z

#### Usage

    MaximumLikelihood$get_timeseries()

#### Returns

dataframe(tau,z)

------------------------------------------------------------------------

### `MaximumLikelihood$get_timeseries_info()`

Get information for times series

#### Usage

    MaximumLikelihood$get_timeseries_info()

#### Returns

list(tbeg,tend,datename,timename,statename,estimation)

------------------------------------------------------------------------

### `MaximumLikelihood$get_plot_info()`

Get information for plotting

#### Usage

    MaximumLikelihood$get_plot_info()

#### Returns

list(font,file,theme,labels)

------------------------------------------------------------------------

### `MaximumLikelihood$get_plot_colors()`

Get colors for plotting

#### Usage

    MaximumLikelihood$get_plot_colors()

#### Returns

list(red,ylw,grn,cyn,blu,mgn,gry,background,font)

------------------------------------------------------------------------

### `MaximumLikelihood$get_flags()`

Get flags for plotting and copying

#### Usage

    MaximumLikelihood$get_flags()

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `MaximumLikelihood$LogLikelihood()`

Calculate the log likelihood of estimates

#### Usage

    MaximumLikelihood$LogLikelihood(
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      df = NULL,
      taucol = NULL,
      zcol = NULL
    )

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index of a column containing z

#### Returns

list(rho,mu,sigma,lnL,k,alpha,m1)

------------------------------------------------------------------------

### `MaximumLikelihood$Estimates()`

Calculate unrestricted or restricted maximum likelihood estimates

#### Usage

    MaximumLikelihood$Estimates(
      df = NULL,
      taucol = NULL,
      zcol = NULL,
      rhor = NULL,
      mur = NULL,
      sigmar = NULL,
      rhos = NULL,
      mus = NULL,
      sigmas = NULL
    )

#### Arguments

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index of a column containing z

- `rhor`:

  constant to fix the rate parameter 0\<=rhor\<inf

- `mur`:

  constant to fix the location parameter -inf\<mur\<inf

- `sigmar`:

  constant to fix scale parameter -inf\<sigmar\<inf

- `rhos`:

  starting value for the rate parameter 0\<=rhos\<inf

- `mus`:

  starting value for the location parameter -inf\<mus\<inf

- `sigmas`:

  starting value the scale parameter -inf\<sigmas\<inf

#### Returns

list(rhohat,muhat,sigmahat,lnLu,ku,alphau,m1) or
list(rhor,mur,sigmar,lnLr,kr,alphar,m1)

------------------------------------------------------------------------

### `MaximumLikelihood$GoodnessOfFit()`

Compare maximum likeliood estimates with invariant and scaled brownian
motion estimates

#### Usage

    MaximumLikelihood$GoodnessOfFit(
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      df = NULL,
      taucol = NULL,
      zcol = NULL
    )

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index of a column containing z

#### Returns

list(theta,theta_i,theta_s,Inv,SBM) with Inv=list(R2,Pval) and
SBM=list(R2,Pval)

------------------------------------------------------------------------

### `MaximumLikelihood$LikelihoodRatioTest()`

Test for significant difference between restricted and unrestricted
likelihoods

#### Usage

    MaximumLikelihood$LikelihoodRatioTest(
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      df = NULL,
      taucol = NULL,
      zcol = NULL
    )

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index of a column containing z

#### Returns

list(theta_u,theta,R2,Pval)

------------------------------------------------------------------------

### `MaximumLikelihood$PlotTimeSeries()`

Plot time series

#### Usage

    MaximumLikelihood$PlotTimeSeries(
      df = NULL,
      taucol = NULL,
      zcol = NULL,
      tbeg = NULL,
      tend = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL
    )

#### Arguments

- `df`:

  data frame containing columns for tau and z

- `taucol`:

  index of a column containing tau

- `zcol`:

  index a column containing z

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

#### Returns

plot

------------------------------------------------------------------------

### `MaximumLikelihood$PlotEstimates()`

Plot time series with estimated means and variances

#### Usage

    MaximumLikelihood$PlotEstimates(
      tbeg = NULL,
      tend = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL
    )

#### Arguments

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

#### Returns

plot
