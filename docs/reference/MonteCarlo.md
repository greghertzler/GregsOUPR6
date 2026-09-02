# R6 class implementing Monte Carlo simulation of the Ornstein-Uhlenbeck Process

Monte Carlo simulations of Forward Paths, Bounded Paths, Backward Paths,
Probabilities, Options and Passage Times. Forward, Backward and Bounded
Paths are simulated by a 4th order Runge-Kutta method and by using the
stochastic integral equation. Probabilities, Options and Passage Times
are various ways of binning and counting the Paths.

## Methods:

    y stochastic
      ForwardPaths
      Mean
      Variance
      Density
      Probability
      DoubleIntegral
    x stochastic
      BackwardPaths
      Option
    t stochastic
      ForwardPaths
      BoundedPaths
      VisitingTimeModeMedianMean
      VisitingTimePercentiles
      VisitingTimeDensity
      VisitingTimeProbability
      FirstPassageTimeModeMedianMean
      FirstPassageTimePercentiles
      FirstPassageTimeDensity
      FirstPassageTimeProbability

## Plots:

      PlotForwardPaths
      PlotBackwardPaths
      PlotBoundedPaths
      PlotMean
      PlotVariance
      PlotDensity
      PlotProbability
      PlotDoubleIntegral
      PlotOption
      PlotVisitingTimeModeMedianMean
      PlotVisitingTimePercentiles
      PlotVisitingTimeDensity
      PlotVisitingTimeProbability
      PlotFirstPassageTimeModeMedianMean
      PlotFirstPassageTimePercentiles
      PlotFirstPassageTimeDensity
      PlotFirstPassageTimeProbability

## Arguments of functions:

      All arguments are optional in all functions.
    OUP parameters
      rho:    rate parameter 0<=rho<inf
      mu:     location parameter -inf<mu<inf
      sigma:  scale parameter -inf<sigma<inf
    y stochastic
      t:      vector of times s<=t<inf
      x:      initial state -inf<x<inf
      psi:    <=0 for integral -inf to y,
               >0 for integral y to inf
    x stochastic
      s:      vector of times -inf<s<t
      y:      terminal state -inf<y<inf
      r:      discount rate -inf<r<inf
      phi:    <=0 for exit option,
               >0 for entry option
    t stochastic
      t:      vector of times s<=t<inf
      k:      decision threshold -inf<k<int
      x:      initial state -inf<x<inf
      Ppct:   probability for a percentile 0.01<=Ppct<=0.99
    path
      paths:  number of paths 1<paths<1,000,000
      skip:   subdivide time interval but report at times t 1<=skip<=50
      seed:   seed for random number generators -inf<seed<inf
      method: 4 for 4th order Runge-Kutta, otherwise integral equation

## Usage:

The MonteCarlo object must first be instantiated before its methods are
called. There are two ways. The first way instantiates the OUProcess
object and calls a function to get a pointer:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      FD <- OUP$get_FiniteDifference()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()

The MonteCarlo object will coordinate arguments to functions with the
Analytial, FiniteDifference and MaximumLikelihood objects. The second
way instantiates the MonteCarlo object by itself with no coordination:

      MC <- MonteCarlol$new()

Once the object is instantiated, its methods can be called, to simulate
100,000 Forward Paths, for example:

      MC$ForwardPaths(paths=100000)

The plot methods can be used to customize the plots, with a title, for
example:

      A$PlotForwardPaths(title="My Paths Forward")

An attempt to plot 100,000 paths would choke the computer, so there are
tricks. One is to select a hundred or so paths for the plot. Another is
to plot heat maps, just like in a weather report.

Other functions and methods are called in the same way. To see all the
possibilities, check out the demos below.

## Demos:

Demonstration scripts are in files in the 'demo' directory. The most
convenient way to list and run demos are the commands:

      OUPDemoList()
      OUPDemoRun(<number of demo in list>)

Entering the demos by number in the list saves typing.

## Discussion:

Monte Carlo simulation of the Ornstein-Uhlenbeck Process can be done
with either of two methods: numerically integrating the stochastic
differential equation, or calculating the stochastic integral equation.
The stochastic differential equation is shocked by a Wiener Process,
simulated as sigma \* dt^0.5 \* epsilon, where sigma \* dt^0.5 is the
square root of the instantaneous variance and epsilon are draws from a
standard normal density. The stochastic integral equation is shocked by
the integral of the Wiener Process, or H \* epsilon, where H is the
square-root of the variance over a longer time interval.

Numerically integrating the stochastic differential equation uses the
Euler, Marayuma or Runge-Kutta schemes. The Euler and Marayuma schemes
are first and third-order schemes and less accurate than the
fourth-order Runge-Kutta scheme.

The fourth-order Runge-Kutta scheme is standard practice. In tests,
shocked by the same draws from a standard normal density, the paths from
the stochastic integral equation and the fourth-order Runge-Kutta scheme
were the same to within five or six significant digits. To compare the
fourth-order Runge- Kutta scheme with the integral equation:

      MC <- MonteCarlo$new()
      rk <- MC$ForwardPaths(paths=100000,skip=10,method=4)[[1]]
      ie <- MC$ForwardPaths(method=1)[[1]]
      dif <- rk-ie
      max(dif)
      min(dif)
      sum(dif)

A Forward Path starts from the backward state at the backward time and
goes forward. A single path, sampled from all possible paths, is a
Sample Path. Just like the flea trying to understand the elephant, a
sample Path is enough for Maximum Likelihood Estimation of the
Ornstein-Uhlenbeck Process. An ensemble of paths can be counted to
approximate Transition Densities and Probabilities and Visiting Time
Densities and Probabilities. The larger the ensemble, the better the
approximations.

Forward Paths are continuous in the state and time which makes them very
hard to count. Instead, Forward Paths are treated as discrete and put
into bins. Counting the number of Forward Paths in each bin and dividing
by the total number of paths approximates Transition Densities. Summing
the Transition Densities approximates Transition Probabilities. Summing
again approximates Double Integrals. Double Integrals are a curiosity.
If time runs backwards, they become Options.

A Forward Path begins from a known state and travels forward into an
uncertain future. A Backward Path ends with a known state and trudges
backward into an uncertain past. Turning around and travelling back to
the future resolves the uncertainty over time. An example is a Bayesian
analysis which begins with a Diffuse Prior and ends with certainty.
Another example is an Option. Simulating and counting Backward Paths
approximates Prior Densities, Prior Probabilities and Options. To
compare Monte Carlo and Analytical Options:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()
      ML$Estimates()
      ao <- A$Option()[[1]]
      mco <- MC$Option(paths=1000000)[[1]]
      dif <- ao-mco
      max(dif)
      min(dif)
      sum(dif)

If we count at right angles–in the time direction instead of the state
direction–Forward Paths become Visiting Time Densities and
Probabilities. Bounded Paths become First Passage Time Densities and
Probabilities. Monte Carlo and Analytical Visiting and First Passage
Times can also be compared:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      MC <- OUP$get_MonteCarlo()
      ad <- A$PassageTimeDensity(omega=0)[[1]]
      mcd <- MC$VisitingTimeDensity(paths=1000000)[[1]]
      dif <- ad-mcd
      max(dif)
      min(dif)
      sum(dif)

Of course, the question is, 'Why bother?' We have Analytical formulas to
do the counting. One reason is to explain the formulas. First Passage
Times make more sense if you plot Bounded Paths and count the number of
paths that have crossed the threshold. Even in prestigious journal
articles, the first and, possibly, only plot will be a Monte Carlo
simulation.

Another reason is to validate the formulas. Although an Analytical
formula may calculate a thousand times faster than a Monte Carlo
simulation, arriving at approximately the same answer both ways is
reassuring.

## Methods

### Public methods

- [`MonteCarlo$new()`](#method-MonteCarlo-initialize)

- [`MonteCarlo$set_oup_params()`](#method-MonteCarlo-set_oup_params)

- [`MonteCarlo$set_y_stoch_args()`](#method-MonteCarlo-set_y_stoch_args)

- [`MonteCarlo$set_x_stoch_args()`](#method-MonteCarlo-set_x_stoch_args)

- [`MonteCarlo$set_t_stoch_args()`](#method-MonteCarlo-set_t_stoch_args)

- [`MonteCarlo$set_path_args()`](#method-MonteCarlo-set_path_args)

- [`MonteCarlo$set_plot_args()`](#method-MonteCarlo-set_plot_args)

- [`MonteCarlo$set_plot_info()`](#method-MonteCarlo-set_plot_info)

- [`MonteCarlo$set_plot_type()`](#method-MonteCarlo-set_plot_type)

- [`MonteCarlo$set_flags()`](#method-MonteCarlo-set_flags)

- [`MonteCarlo$get_all()`](#method-MonteCarlo-get_all)

- [`MonteCarlo$get_oup_params()`](#method-MonteCarlo-get_oup_params)

- [`MonteCarlo$get_y_stoch_args()`](#method-MonteCarlo-get_y_stoch_args)

- [`MonteCarlo$get_x_stoch_args()`](#method-MonteCarlo-get_x_stoch_args)

- [`MonteCarlo$get_t_stoch_args()`](#method-MonteCarlo-get_t_stoch_args)

- [`MonteCarlo$get_path_args()`](#method-MonteCarlo-get_path_args)

- [`MonteCarlo$get_plot_args()`](#method-MonteCarlo-get_plot_args)

- [`MonteCarlo$get_plot_info()`](#method-MonteCarlo-get_plot_info)

- [`MonteCarlo$get_plot_colors()`](#method-MonteCarlo-get_plot_colors)

- [`MonteCarlo$get_plot_types()`](#method-MonteCarlo-get_plot_types)

- [`MonteCarlo$get_flags()`](#method-MonteCarlo-get_flags)

- [`MonteCarlo$axes_y_stoch()`](#method-MonteCarlo-axes_y_stoch)

- [`MonteCarlo$axes_x_stoch()`](#method-MonteCarlo-axes_x_stoch)

- [`MonteCarlo$axes_t_stoch()`](#method-MonteCarlo-axes_t_stoch)

- [`MonteCarlo$sync_yxt_stoch()`](#method-MonteCarlo-sync_yxt_stoch)

- [`MonteCarlo$undo_clear()`](#method-MonteCarlo-undo_clear)

- [`MonteCarlo$undo_save()`](#method-MonteCarlo-undo_save)

- [`MonteCarlo$undo_undo()`](#method-MonteCarlo-undo_undo)

- [`MonteCarlo$ForwardPaths()`](#method-MonteCarlo-ForwardPaths)

- [`MonteCarlo$BackwardPaths()`](#method-MonteCarlo-BackwardPaths)

- [`MonteCarlo$BoundedPaths()`](#method-MonteCarlo-BoundedPaths)

- [`MonteCarlo$Mean()`](#method-MonteCarlo-Mean)

- [`MonteCarlo$Variance()`](#method-MonteCarlo-Variance)

- [`MonteCarlo$Density()`](#method-MonteCarlo-Density)

- [`MonteCarlo$Probability()`](#method-MonteCarlo-Probability)

- [`MonteCarlo$DoubleIntegral()`](#method-MonteCarlo-DoubleIntegral)

- [`MonteCarlo$Option()`](#method-MonteCarlo-Option)

- [`MonteCarlo$VisitingTimeModeMedianMean()`](#method-MonteCarlo-VisitingTimeModeMedianMean)

- [`MonteCarlo$VisitingTimePercentiles()`](#method-MonteCarlo-VisitingTimePercentiles)

- [`MonteCarlo$VisitingTimeDensity()`](#method-MonteCarlo-VisitingTimeDensity)

- [`MonteCarlo$VisitingTimeProbability()`](#method-MonteCarlo-VisitingTimeProbability)

- [`MonteCarlo$FirstPassageTimeModeMedianMean()`](#method-MonteCarlo-FirstPassageTimeModeMedianMean)

- [`MonteCarlo$FirstPassageTimePercentiles()`](#method-MonteCarlo-FirstPassageTimePercentiles)

- [`MonteCarlo$FirstPassageTimeDensity()`](#method-MonteCarlo-FirstPassageTimeDensity)

- [`MonteCarlo$FirstPassageTimeProbability()`](#method-MonteCarlo-FirstPassageTimeProbability)

- [`MonteCarlo$PlotForwardPaths()`](#method-MonteCarlo-PlotForwardPaths)

- [`MonteCarlo$PlotBackwardPaths()`](#method-MonteCarlo-PlotBackwardPaths)

- [`MonteCarlo$PlotBoundedPaths()`](#method-MonteCarlo-PlotBoundedPaths)

- [`MonteCarlo$PlotMean()`](#method-MonteCarlo-PlotMean)

- [`MonteCarlo$PlotVariance()`](#method-MonteCarlo-PlotVariance)

- [`MonteCarlo$PlotDensity()`](#method-MonteCarlo-PlotDensity)

- [`MonteCarlo$PlotProbability()`](#method-MonteCarlo-PlotProbability)

- [`MonteCarlo$PlotDoubleIntegral()`](#method-MonteCarlo-PlotDoubleIntegral)

- [`MonteCarlo$PlotOption()`](#method-MonteCarlo-PlotOption)

- [`MonteCarlo$PlotVisitingTimeModeMedianMean()`](#method-MonteCarlo-PlotVisitingTimeModeMedianMean)

- [`MonteCarlo$PlotVisitingTimePercentiles()`](#method-MonteCarlo-PlotVisitingTimePercentiles)

- [`MonteCarlo$PlotVisitingTimeDensity()`](#method-MonteCarlo-PlotVisitingTimeDensity)

- [`MonteCarlo$PlotVisitingTimeProbability()`](#method-MonteCarlo-PlotVisitingTimeProbability)

- [`MonteCarlo$PlotFirstPassageTimeModeMedianMean()`](#method-MonteCarlo-PlotFirstPassageTimeModeMedianMean)

- [`MonteCarlo$PlotFirstPassageTimePercentiles()`](#method-MonteCarlo-PlotFirstPassageTimePercentiles)

- [`MonteCarlo$PlotFirstPassageTimeDensity()`](#method-MonteCarlo-PlotFirstPassageTimeDensity)

- [`MonteCarlo$PlotFirstPassageTimeProbability()`](#method-MonteCarlo-PlotFirstPassageTimeProbability)

------------------------------------------------------------------------

### `MonteCarlo$new()`

Create a MonteCarlo object

#### Usage

    MonteCarlo$new(OUP = NULL)

#### Arguments

- `OUP`:

  pointer set by the OUProcess object

#### Returns

A new MonteCarlo object

------------------------------------------------------------------------

### `MonteCarlo$set_oup_params()`

Set OUP parameters

#### Usage

    MonteCarlo$set_oup_params(rho = NULL, mu = NULL, sigma = NULL, who = NULL)

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

list(rho,mu,sigma)

------------------------------------------------------------------------

### `MonteCarlo$set_y_stoch_args()`

Set arguments for y as a stochastic state

#### Usage

    MonteCarlo$set_y_stoch_args(
      t = NULL,
      y = NULL,
      x = NULL,
      psi = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `x`:

  initial state -inf\<x\<inf

- `psi`:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- `who`:

  object id of sender

#### Returns

list(t,y,x,psi)

------------------------------------------------------------------------

### `MonteCarlo$set_x_stoch_args()`

Set arguments for x as a stochastic state

#### Usage

    MonteCarlo$set_x_stoch_args(
      s = NULL,
      x = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<t

- `x`:

  vector of n states -inf\<x\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `who`:

  object id of sender

#### Returns

list(s,x,y,r,phi)

------------------------------------------------------------------------

### `MonteCarlo$set_t_stoch_args()`

Set t stochastic arguments

#### Usage

    MonteCarlo$set_t_stoch_args(
      t = NULL,
      k = NULL,
      x = NULL,
      omega = NULL,
      Ppct = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `k`:

  threshold -inf\<k\<inf

- `x`:

  initial state -inf\<x\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `Ppct`:

  passage time probability for a percentile 0.01\<=Ppct\<=0.99

- `who`:

  object id of sender

#### Returns

list(t,k,x)

------------------------------------------------------------------------

### `MonteCarlo$set_path_args()`

Set path arguments

#### Usage

    MonteCarlo$set_path_args(paths = NULL, skip = NULL, seed = NULL, method = NULL)

#### Arguments

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time intervals but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

#### Returns

list(paths,skip,seed,method)

------------------------------------------------------------------------

### `MonteCarlo$set_plot_args()`

Set plot arguments

#### Usage

    MonteCarlo$set_plot_args(
      pmax = NULL,
      ptmax = NULL,
      first = NULL,
      last = NULL,
      zbeg = NULL,
      zend = NULL,
      who = NULL
    )

#### Arguments

- `pmax`:

  maximum transition density

- `ptmax`:

  maximum visiting time and first passage time densities

- `first`:

  first path to plot

- `last`:

  last path to plot

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

- `who`:

  object id of caller

#### Returns

list(pmax,ptmax,first,last,zbeg,zend)

------------------------------------------------------------------------

### `MonteCarlo$set_plot_info()`

Set information for plotting

#### Usage

    MonteCarlo$set_plot_info(
      fontfamily = NULL,
      fontsize = NULL,
      fileformat = NULL,
      filewidth = NULL,
      fileheight = NULL,
      theme = NULL,
      opaque = NULL,
      walls = NULL,
      floor = NULL,
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

- `walls`:

  3D walls TRUE or FALSE

- `floor`:

  3D floor TRUE or FALSE

- `labels`:

  title and parameters TRUE or FALSE

- `who`:

  object id of sender

#### Returns

list(font,file,theme,3D,labels)

------------------------------------------------------------------------

### `MonteCarlo$set_plot_type()`

Set plot types by group

#### Usage

    MonteCarlo$set_plot_type(type = NULL, group = NULL)

#### Arguments

- `type`:

  a number identifying the type for a group

- `group`:

  identifier for groups of similar plots

#### Returns

list(type,group)

------------------------------------------------------------------------

### `MonteCarlo$set_flags()`

Set flags for plotting and copying

#### Usage

    MonteCarlo$set_flags(plotit = NULL, copyit = NULL, who = NULL)

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

### `MonteCarlo$get_all()`

Get all arguments and information

#### Usage

    MonteCarlo$get_all()

#### Returns

list(oup_params,y_stoch_args,x_stoch_args,t_stoch_args,path_args,plot_info)

------------------------------------------------------------------------

### `MonteCarlo$get_oup_params()`

Get OUP parameters

#### Usage

    MonteCarlo$get_oup_params()

#### Returns

list(rho,mu,sigma)

------------------------------------------------------------------------

### `MonteCarlo$get_y_stoch_args()`

Get arguments for y as a stochastic state

#### Usage

    MonteCarlo$get_y_stoch_args()

#### Returns

list(t,y,x,psi)

------------------------------------------------------------------------

### `MonteCarlo$get_x_stoch_args()`

Get arguments for x as a stochastic state

#### Usage

    MonteCarlo$get_x_stoch_args()

#### Returns

list(s,y,r,phi)

------------------------------------------------------------------------

### `MonteCarlo$get_t_stoch_args()`

Get arguments for t as a stochastic state

#### Usage

    MonteCarlo$get_t_stoch_args()

#### Returns

list(t,k,x,z,omega,Ppct)

------------------------------------------------------------------------

### `MonteCarlo$get_path_args()`

get path arguments

#### Usage

    MonteCarlo$get_path_args()

#### Returns

list(paths,skip,seed,method)

------------------------------------------------------------------------

### `MonteCarlo$get_plot_args()`

get plot arguments

#### Usage

    MonteCarlo$get_plot_args()

#### Returns

list(pmax,ptmax,first,last,zbeg,zend)

------------------------------------------------------------------------

### `MonteCarlo$get_plot_info()`

Get information for plotting

#### Usage

    MonteCarlo$get_plot_info()

#### Returns

list(font,file,theme,3D,labels)

------------------------------------------------------------------------

### `MonteCarlo$get_plot_colors()`

Get colors for plotting

#### Usage

    MonteCarlo$get_plot_colors()

#### Returns

list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)

------------------------------------------------------------------------

### `MonteCarlo$get_plot_types()`

Get current types for plot routines

#### Usage

    MonteCarlo$get_plot_types()

#### Returns

(list(types,description))

------------------------------------------------------------------------

### `MonteCarlo$get_flags()`

Get flags for plotting and copying

#### Usage

    MonteCarlo$get_flags()

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `MonteCarlo$axes_y_stoch()`

Scale axes for y stochastic paths

#### Usage

    MonteCarlo$axes_y_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `MonteCarlo$axes_x_stoch()`

Scale axes for x stochastic paths

#### Usage

    MonteCarlo$axes_x_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `MonteCarlo$axes_t_stoch()`

Scale axes for t stochastic paths

#### Usage

    MonteCarlo$axes_t_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `MonteCarlo$sync_yxt_stoch()`

Synchronize arguments

#### Usage

    MonteCarlo$sync_yxt_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `MonteCarlo$undo_clear()`

Clear undo list and save current arguments to list

#### Usage

    MonteCarlo$undo_clear()

#### Returns

1

------------------------------------------------------------------------

### `MonteCarlo$undo_save()`

Save current arguments to undo list

#### Usage

    MonteCarlo$undo_save()

#### Returns

number of argument sets

------------------------------------------------------------------------

### `MonteCarlo$undo_undo()`

Replace current arguments from undo list

#### Usage

    MonteCarlo$undo_undo(updn = NULL)

#### Arguments

- `updn`:

  positive to move up, negative to move down

#### Returns

c(index of this argument set, number of argument sets)

------------------------------------------------------------------------

### `MonteCarlo$ForwardPaths()`

Calculate and plot forward paths

#### Usage

    MonteCarlo$ForwardPaths(
      t = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(forward(m,paths))

------------------------------------------------------------------------

### `MonteCarlo$BackwardPaths()`

Calculate and plot backward paths

#### Usage

    MonteCarlo$BackwardPaths(
      s = NULL,
      y = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<inf

- `y`:

  terminal state -inf\<y\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(backward(m,paths))

------------------------------------------------------------------------

### `MonteCarlo$BoundedPaths()`

Calculate and plot bounded paths

#### Usage

    MonteCarlo$BoundedPaths(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `k`:

  threshold or bound -inf\<k\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(bounded(m,paths))

------------------------------------------------------------------------

### `MonteCarlo$Mean()`

Calculate and plot means

#### Usage

    MonteCarlo$Mean(
      t = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(G(m))

------------------------------------------------------------------------

### `MonteCarlo$Variance()`

Calculate and plot variances

#### Usage

    MonteCarlo$Variance(
      t = NULL,
      x = NULL,
      rho = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(H2(m))

------------------------------------------------------------------------

### `MonteCarlo$Density()`

Calculate and plot densities

#### Usage

    MonteCarlo$Density(
      t = NULL,
      y = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(p(m,n))

------------------------------------------------------------------------

### `MonteCarlo$Probability()`

Calculate and plot probabilities

#### Usage

    MonteCarlo$Probability(
      t = NULL,
      y = NULL,
      x = NULL,
      psi = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `x`:

  initial state -inf\<x\<inf

- `psi`:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(P(m,n))

------------------------------------------------------------------------

### `MonteCarlo$DoubleIntegral()`

Calculate and plot double integrals

#### Usage

    MonteCarlo$DoubleIntegral(
      t = NULL,
      y = NULL,
      x = NULL,
      psi = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times -inf\<t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `x`:

  initial state -inf\<x\<inf

- `psi`:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(PP(m,n))

------------------------------------------------------------------------

### `MonteCarlo$Option()`

Calculate and plot option prices

#### Usage

    MonteCarlo$Option(
      s = NULL,
      x = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<t

- `x`:

  vector of n states -inf\<x\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(OO(mxn))

------------------------------------------------------------------------

### `MonteCarlo$VisitingTimeModeMedianMean()`

Calculate and plot visiting time mode, median and mean

#### Usage

    MonteCarlo$VisitingTimeModeMedianMean(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(vtmmm(3x3))

------------------------------------------------------------------------

### `MonteCarlo$VisitingTimePercentiles()`

Calculate and plot visiting time percentiles

#### Usage

    MonteCarlo$VisitingTimePercentiles(
      t = NULL,
      k = NULL,
      x = NULL,
      Ppct = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `Ppct`:

  probability for a percentile 0.01\<=Ppct\<=0.99

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(vtpct(3x3))

------------------------------------------------------------------------

### `MonteCarlo$VisitingTimeDensity()`

Calculate and plot visiting time densities

#### Usage

    MonteCarlo$VisitingTimeDensity(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(pv(m))

------------------------------------------------------------------------

### `MonteCarlo$VisitingTimeProbability()`

Calculate and plot visiting time probabilities

#### Usage

    MonteCarlo$VisitingTimeProbability(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(Pv(m))

------------------------------------------------------------------------

### `MonteCarlo$FirstPassageTimeModeMedianMean()`

Calculate and plot first passage time mode, median and mean

#### Usage

    MonteCarlo$FirstPassageTimeModeMedianMean(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(fptmmm(3x3))

------------------------------------------------------------------------

### `MonteCarlo$FirstPassageTimePercentiles()`

Calculate and plot first passage time percentiles

#### Usage

    MonteCarlo$FirstPassageTimePercentiles(
      t = NULL,
      k = NULL,
      x = NULL,
      Ppct = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `Ppct`:

  probability for a percentile 0.01\<=Ppct\<=0.99

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(fptpct(3x3))

------------------------------------------------------------------------

### `MonteCarlo$FirstPassageTimeDensity()`

Calculate and plot first passage time densities

#### Usage

    MonteCarlo$FirstPassageTimeDensity(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(pf(m))

------------------------------------------------------------------------

### `MonteCarlo$FirstPassageTimeProbability()`

Calculate and plot first passage time probabilities

#### Usage

    MonteCarlo$FirstPassageTimeProbability(
      t = NULL,
      k = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      paths = NULL,
      skip = NULL,
      seed = NULL,
      method = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `paths`:

  number of paths 1\<paths\<1,000,000

- `skip`:

  subdivide time interval but report at times t 1\<=skip\<=50

- `seed`:

  seed for random number generators -inf\<seed\<inf

- `method`:

  4 for 4th order Runge-Kutta, otherwise integral equation

- `who`:

  object id of caller

#### Returns

list(Pf(m))

------------------------------------------------------------------------

### `MonteCarlo$PlotForwardPaths()`

Plot forward paths

#### Usage

    MonteCarlo$PlotForwardPaths(
      type = NULL,
      first = NULL,
      last = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -3,...,0, or 'n','p','d' for next, previous, default

- `first`:

  first path to plot

- `last`:

  last path to plot

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotBackwardPaths()`

Plot backward paths

#### Usage

    MonteCarlo$PlotBackwardPaths(
      type = NULL,
      first = NULL,
      last = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      sbeg = NULL,
      send = NULL
    )

#### Arguments

- `type`:

  = -3,...,0, or 'n','p','d' for next, previous, default

- `first`:

  first path to plot

- `last`:

  last path to plot

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `sbeg`:

  begin value for time axis

- `send`:

  end value for time axis

- `copyit`:

  TRUE or FALSE

#### Returns

plot Plot bounded paths

------------------------------------------------------------------------

### `MonteCarlo$PlotBoundedPaths()`

#### Usage

    MonteCarlo$PlotBoundedPaths(
      type = NULL,
      first = NULL,
      last = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -3,...,0, or 'n','p','d' for next, previous, default

- `first`:

  first path to plot

- `last`:

  last path to plot

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotMean()`

Plot means

#### Usage

    MonteCarlo$PlotMean(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotVariance()`

Plot variances

#### Usage

    MonteCarlo$PlotVariance(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotDensity()`

Plot transition densities

#### Usage

    MonteCarlo$PlotDensity(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      ybeg = NULL,
      yend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `ybeg`:

  begin value for state axis

- `yend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotProbability()`

Plot transition probabilities

#### Usage

    MonteCarlo$PlotProbability(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      ybeg = NULL,
      yend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `ybeg`:

  begin value for state axis

- `yend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotDoubleIntegral()`

Plot double integrals of transition densities

#### Usage

    MonteCarlo$PlotDoubleIntegral(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      ybeg = NULL,
      yend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `ybeg`:

  begin value for state axis

- `yend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotOption()`

Plot options

#### Usage

    MonteCarlo$PlotOption(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      sbeg = NULL,
      send = NULL,
      xbeg = NULL,
      xend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `sbeg`:

  begin value for time axis

- `send`:

  end value for time axis

- `xbeg`:

  begin value for state axis

- `xend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotVisitingTimeModeMedianMean()`

Plot visiting time mode, median and mean

#### Usage

    MonteCarlo$PlotVisitingTimeModeMedianMean(
      type = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotVisitingTimePercentiles()`

Plot visiting time lower, median and upper percentiles

#### Usage

    MonteCarlo$PlotVisitingTimePercentiles(
      type = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotVisitingTimeDensity()`

Plot visiting time densities

#### Usage

    MonteCarlo$PlotVisitingTimeDensity(
      type = NULL,
      pmax = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotVisitingTimeProbability()`

Plot visiting time probabilities

#### Usage

    MonteCarlo$PlotVisitingTimeProbability(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotFirstPassageTimeModeMedianMean()`

Plot first passage time mode, median and mean

#### Usage

    MonteCarlo$PlotFirstPassageTimeModeMedianMean(
      type = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotFirstPassageTimePercentiles()`

Plot first passage time lower, median and upper percentiles

#### Usage

    MonteCarlo$PlotFirstPassageTimePercentiles(
      type = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      tbeg = NULL,
      tend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotFirstPassageTimeDensity()`

Plot first passage time densities

#### Usage

    MonteCarlo$PlotFirstPassageTimeDensity(
      type = NULL,
      pmax = NULL,
      ptmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `ptmax`:

  maximum visiting time and first passage time densities

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `MonteCarlo$PlotFirstPassageTimeProbability()`

Plot first passage time probabilities

#### Usage

    MonteCarlo$PlotFirstPassageTimeProbability(
      type = NULL,
      pmax = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zaxis = NULL,
      tbeg = NULL,
      tend = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = 0, 1, or 'n','p','d' for next, previous, default

- `pmax`:

  maximum transition density

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zaxis`:

  text for z-axis label

- `tbeg`:

  begin value for time axis

- `tend`:

  end value for time axis

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot
