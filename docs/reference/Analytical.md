# R6 class implementing Analytical formulas for the Ornstein-Uhlenbeck Process

A set of functions to calculate probabilities, option prices, visiting
times, first passage times and decision thresholds–everything for simple
exit and entry options in a Real Options Analysis.

## Formulas and Methods:

    z stochastic
      Drift
      Diffusion
    y stochastic
      Mean
      Variance
      Density
      Probability
      DoubleIntegral
    x stochastic
      Option
      OptionEnvelope
      DecisionThreshold
      Obligation
    t stochastic
      PassageTimeModeMedianMean
      PassageTimePercentiles
      PassageTimeDensity
      PassageTimeProbability

## Plots

      PlotDrift
      PlotDiffusion
      PlotMean
      PlotVariance
      PlotDensity
      PlotProbability
      PlotDoubleIntegral
      PlotOption
      PlotOptionEnvelope
      PlotDecisionThreshold
      PlotObligation
      PlotPassageTimeModeMedianMean
      PlotPassageTimePercentiles
      PlotPassageTimeDensity
      PlotPassageTimeProbability

## Arguments of functions:

      All arguments are optional in all functions.
    OUP parameters
      rho:    rate parameter 0<=rho<inf
      mu:     location parameter -inf<mu<inf
      sigma:  scale parameter -inf<sigma<inf
    z stochastic
      z:      vector of states -inf<z<inf
    y stochastic
      t:      vector of times s<=t<inf
      y:      vector of states -inf<y<inf
      s:      initial time -inf<s<inf
      x:      initial state -inf<x<inf
      psi:    <=0 for integral -inf to y,
               >0 for integral y to inf
      eps:    proportion remaining after convergence 0<=eps<=1
    x stochastic
      s:      vector of times -inf<s<t
      x:      vector of states -inf<x<inf
      t:      terminal time -inf<t<inf
      y:      terminal state -inf<y<inf
      r:      discount rate -inf<r<inf
      phi:    <=0 for exit option,
               >0 for entry option
    t stochastic
      t:      vector of times s<=t<inf
      k:      decision threshold -inf<k<int
      s:      initial time -inf<s<inf
      x:      initial state -inf<x<inf
      omega:  degree of irreversibility 0<=omega<=1
      Ppct:   passage time probability for a percentile 0.01<=Ppct<=0.99

## Usage:

The Analytical object must first be instantiated before its methods are
called. There are two ways. The first way instantiates the OUProcess
object and calls a function to get a pointer:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      FD <- OUP$get_FiniteDifference()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()

The Analytical object will coordinate arguments to functions with the
FiniteDifference, MaximumLikelihood and MonteCarlo objects. The second
way instantiates the Analytical object by itself with no coordination:

      A <- Analytical$new()

Once the object is instantiated, its methods can be called, to calculate
and plot a Decision Threshold, for example:

      A$DecisionThreshold()

The plot methods can be used to customize the plots, with a title, for
example:

      A$PlotDecisionThreshold(title="My Decision")

Other functions and methods are called in the same way. To see all the
possibilities, check out the demos below.

## Demos:

Demonstration scripts are in files in the 'demo' directory. The most
convenient way to list and run demos are the commands:

      OUPDemoList()
      OUPDemoRun(<number of demo in list>)

Entering the demos by number in the list saves typing.

## Discussion:

Mathematical formulas and convergence criteria assume exact arithmetic.
Floating-point arithmetic is another matter. It may overflow, underflow
or have extreme cancellation. In principle, everything can be calculated
by Monte Carlo simulation or Finite Difference methods. These solve
continuous equations at discrete nodes, creating another source of
error. Monte Carlo simulation of passage times can be accurate or
biased, sometimes wildly so. The Finite Difference method can be
accurate or give reasonable looking solutions that diverge from the true
answers.

Even if great care is taken, it is a good idea to calculate a solution
in more than one way. This is one reason for analytical formulas.

Another reason is to speed the calculations. Compared with the Finite
Difference method, an analytical formula can calculate option prices 2
to 3 times quicker. Compared with Monte Carlo simulation, an analytical
formula can calculate 800 to 1,000 times quicker.

The bottom line of a real options analysis is the decision threshold and
the time to crossing the decision threshold. The many functions
available are used in the calculations, but they all accumulate to the
bottom line:

      A <- Analytical$new()
      A$DecisionThreshold()
      A$PassageTimeModeMedianMean()
      A$PassageTimePercentiles()

Perhaps the better measure of passage times are the percentiles. The
mean does not exist if rho is zero and the variance does not exist if
rho is small. The mode and median always exist. Percentiles, which
include the median, always exist. The formulas don't calculate the
variance.

The passage time calculations are challenging. For example, the passage
time density is too complicated to explain in this short discussion. But
don't be surprised if it is bi-modal or even negative. Care must be
taken in interpreting passage times.

## Methods

### Public methods

- [`Analytical$new()`](#method-Analytical-initialize)

- [`Analytical$set_oup_params()`](#method-Analytical-set_oup_params)

- [`Analytical$set_z_stoch_args()`](#method-Analytical-set_z_stoch_args)

- [`Analytical$set_y_stoch_args()`](#method-Analytical-set_y_stoch_args)

- [`Analytical$set_x_stoch_args()`](#method-Analytical-set_x_stoch_args)

- [`Analytical$set_t_stoch_args()`](#method-Analytical-set_t_stoch_args)

- [`Analytical$set_plot_args()`](#method-Analytical-set_plot_args)

- [`Analytical$set_plot_info()`](#method-Analytical-set_plot_info)

- [`Analytical$set_plot_type()`](#method-Analytical-set_plot_type)

- [`Analytical$set_flags()`](#method-Analytical-set_flags)

- [`Analytical$get_all()`](#method-Analytical-get_all)

- [`Analytical$get_oup_params()`](#method-Analytical-get_oup_params)

- [`Analytical$get_z_stoch_args()`](#method-Analytical-get_z_stoch_args)

- [`Analytical$get_y_stoch_args()`](#method-Analytical-get_y_stoch_args)

- [`Analytical$get_x_stoch_args()`](#method-Analytical-get_x_stoch_args)

- [`Analytical$get_t_stoch_args()`](#method-Analytical-get_t_stoch_args)

- [`Analytical$get_plot_args()`](#method-Analytical-get_plot_args)

- [`Analytical$get_plot_info()`](#method-Analytical-get_plot_info)

- [`Analytical$get_plot_colors()`](#method-Analytical-get_plot_colors)

- [`Analytical$get_plot_types()`](#method-Analytical-get_plot_types)

- [`Analytical$get_flags()`](#method-Analytical-get_flags)

- [`Analytical$axes_z_stoch()`](#method-Analytical-axes_z_stoch)

- [`Analytical$axes_y_stoch()`](#method-Analytical-axes_y_stoch)

- [`Analytical$axes_x_stoch()`](#method-Analytical-axes_x_stoch)

- [`Analytical$axes_t_stoch()`](#method-Analytical-axes_t_stoch)

- [`Analytical$sync_zyxt_stoch()`](#method-Analytical-sync_zyxt_stoch)

- [`Analytical$undo_clear()`](#method-Analytical-undo_clear)

- [`Analytical$undo_save()`](#method-Analytical-undo_save)

- [`Analytical$undo_undo()`](#method-Analytical-undo_undo)

- [`Analytical$Drift()`](#method-Analytical-Drift)

- [`Analytical$Diffusion()`](#method-Analytical-Diffusion)

- [`Analytical$Mean()`](#method-Analytical-Mean)

- [`Analytical$Variance()`](#method-Analytical-Variance)

- [`Analytical$Density()`](#method-Analytical-Density)

- [`Analytical$Probability()`](#method-Analytical-Probability)

- [`Analytical$DoubleIntegral()`](#method-Analytical-DoubleIntegral)

- [`Analytical$Option()`](#method-Analytical-Option)

- [`Analytical$OptionEnvelope()`](#method-Analytical-OptionEnvelope)

- [`Analytical$DecisionThreshold()`](#method-Analytical-DecisionThreshold)

- [`Analytical$Obligation()`](#method-Analytical-Obligation)

- [`Analytical$PassageTimeModeMedianMean()`](#method-Analytical-PassageTimeModeMedianMean)

- [`Analytical$PassageTimePercentiles()`](#method-Analytical-PassageTimePercentiles)

- [`Analytical$PassageTimeDensity()`](#method-Analytical-PassageTimeDensity)

- [`Analytical$PassageTimeProbability()`](#method-Analytical-PassageTimeProbability)

- [`Analytical$PlotDrift()`](#method-Analytical-PlotDrift)

- [`Analytical$PlotDiffusion()`](#method-Analytical-PlotDiffusion)

- [`Analytical$PlotMean()`](#method-Analytical-PlotMean)

- [`Analytical$PlotVariance()`](#method-Analytical-PlotVariance)

- [`Analytical$PlotDensity()`](#method-Analytical-PlotDensity)

- [`Analytical$PlotProbability()`](#method-Analytical-PlotProbability)

- [`Analytical$PlotDoubleIntegral()`](#method-Analytical-PlotDoubleIntegral)

- [`Analytical$PlotOption()`](#method-Analytical-PlotOption)

- [`Analytical$PlotOptionEnvelope()`](#method-Analytical-PlotOptionEnvelope)

- [`Analytical$PlotDecisionThreshold()`](#method-Analytical-PlotDecisionThreshold)

- [`Analytical$PlotObligation()`](#method-Analytical-PlotObligation)

- [`Analytical$PlotPassageTimeModeMedianMean()`](#method-Analytical-PlotPassageTimeModeMedianMean)

- [`Analytical$PlotPassageTimePercentiles()`](#method-Analytical-PlotPassageTimePercentiles)

- [`Analytical$PlotPassageTimeDensity()`](#method-Analytical-PlotPassageTimeDensity)

- [`Analytical$PlotPassageTimeProbability()`](#method-Analytical-PlotPassageTimeProbability)

------------------------------------------------------------------------

### `Analytical$new()`

Create an Analytical object

#### Usage

    Analytical$new(OUP = NULL)

#### Arguments

- `OUP`:

  pointer set by the OUProcess object

#### Returns

A new Analytical object

------------------------------------------------------------------------

### `Analytical$set_oup_params()`

Set OUP parameters

#### Usage

    Analytical$set_oup_params(rho = NULL, mu = NULL, sigma = NULL, who = NULL)

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(rho,mu,sigma)

------------------------------------------------------------------------

### `Analytical$set_z_stoch_args()`

Set z as a stochastic state

#### Usage

    Analytical$set_z_stoch_args(z = NULL)

#### Arguments

- `z`:

  vector of n states -inf\<z\<inf

#### Returns

list(z)

------------------------------------------------------------------------

### `Analytical$set_y_stoch_args()`

Set y as a stochastic state and its arguments

#### Usage

    Analytical$set_y_stoch_args(
      t = NULL,
      y = NULL,
      s = NULL,
      x = NULL,
      psi = NULL,
      eps = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `psi`:

  \<=0 for integral -inf to y, \>0 for integral y to inf

- `eps`:

  proportion remaining after convergence 0\<=eps\<=1

- `who`:

  object id of caller

#### Returns

list(t,y,s,x,psi,eps)

------------------------------------------------------------------------

### `Analytical$set_x_stoch_args()`

Set x as a stochastic state and its arguments

#### Usage

    Analytical$set_x_stoch_args(
      s = NULL,
      x = NULL,
      t = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      b = NULL,
      c = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<t

- `x`:

  vector of n states -inf\<x\<inf

- `t`:

  terminal time -inf\<t\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `b`:

  lump-sum benefit for entry option

- `c`:

  lump-sum cost for exit option

- `who`:

  object id of caller

#### Returns

list(s,x,t,y,r,phi)

------------------------------------------------------------------------

### `Analytical$set_t_stoch_args()`

Set t stochastic arguments

#### Usage

    Analytical$set_t_stoch_args(
      t = NULL,
      k = NULL,
      s = NULL,
      x = NULL,
      z = NULL,
      omega = NULL,
      Ppct = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  threshold -inf\<k\<inf

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `z`:

  vector of n alternate initial states -inf\<z\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `Ppct`:

  passage time probability for a percentile 0.01\<=Ppct\<=0.99

- `who`:

  object id of caller

#### Returns

list(t,k,s,x,z,omega)

------------------------------------------------------------------------

### `Analytical$set_plot_args()`

Set plot arguments

#### Usage

    Analytical$set_plot_args(pmax = NULL, ptmax = NULL, who = NULL)

#### Arguments

- `pmax`:

  maximum transition density

- `ptmax`:

  maximum visiting time and first passage time densities

- `who`:

  object id of caller

#### Returns

list(pmax,ptmax)

------------------------------------------------------------------------

### `Analytical$set_plot_info()`

Set information for plotting

#### Usage

    Analytical$set_plot_info(
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

  object id of caller

#### Returns

list(type,font,file,theme,3D)

------------------------------------------------------------------------

### `Analytical$set_plot_type()`

Set plot types by group

#### Usage

    Analytical$set_plot_type(type = NULL, group = NULL)

#### Arguments

- `type`:

  a number identifying the type for a group

- `group`:

  identifier for groups of similar plots

#### Returns

list(type,group)

------------------------------------------------------------------------

### `Analytical$set_flags()`

Set flags for plotting and copying

#### Usage

    Analytical$set_flags(plotit = NULL, copyit = NULL, who = NULL)

#### Arguments

- `plotit`:

  automatic plot after calculation TRUE or FALSE

- `copyit`:

  copy to clipboard TRUE or FALSE

- `who`:

  object id of caller

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `Analytical$get_all()`

Get all arguments and information

#### Usage

    Analytical$get_all()

#### Returns

list(oup_params,z_stoch_args,y_stoch_args,x_stoch_args,t_stoch_args,plot_info)

------------------------------------------------------------------------

### `Analytical$get_oup_params()`

Get OUP parameters

#### Usage

    Analytical$get_oup_params()

#### Returns

list(rho,mu,sigma)

------------------------------------------------------------------------

### `Analytical$get_z_stoch_args()`

Get z as a stochastic state

#### Usage

    Analytical$get_z_stoch_args()

#### Returns

list(z)

------------------------------------------------------------------------

### `Analytical$get_y_stoch_args()`

Get y as a stochastic state and its arguments

#### Usage

    Analytical$get_y_stoch_args()

#### Returns

list(t,y,s,x,psi,eps)

------------------------------------------------------------------------

### `Analytical$get_x_stoch_args()`

Get x as a stochastic state and its arguments

#### Usage

    Analytical$get_x_stoch_args()

#### Returns

list(s,x,t,y,r,phi)

------------------------------------------------------------------------

### `Analytical$get_t_stoch_args()`

Get t stochastic arguments

#### Usage

    Analytical$get_t_stoch_args()

#### Returns

list(t,k,s,x,z,omega,Ppct)

------------------------------------------------------------------------

### `Analytical$get_plot_args()`

Get plot arguments

#### Usage

    Analytical$get_plot_args()

#### Returns

list(pmax,ptmax)

------------------------------------------------------------------------

### `Analytical$get_plot_info()`

Get information for plotting

#### Usage

    Analytical$get_plot_info()

#### Returns

list(type,font,file,theme,3D,labels)

------------------------------------------------------------------------

### `Analytical$get_plot_colors()`

Get colors for plotting

#### Usage

    Analytical$get_plot_colors()

#### Returns

list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)

------------------------------------------------------------------------

### `Analytical$get_plot_types()`

Get current types for plot routines

#### Usage

    Analytical$get_plot_types()

#### Returns

(list(types,descriptioon))

------------------------------------------------------------------------

### `Analytical$get_flags()`

Get flags for plotting and copying

#### Usage

    Analytical$get_flags()

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `Analytical$axes_z_stoch()`

Scale axes for z stochastic arguments

#### Usage

    Analytical$axes_z_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `Analytical$axes_y_stoch()`

Scale axes for y stochastic arguments

#### Usage

    Analytical$axes_y_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `Analytical$axes_x_stoch()`

Scale axes for x stochastic arguments

#### Usage

    Analytical$axes_x_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `Analytical$axes_t_stoch()`

Scale axes for t stochastic arguments

#### Usage

    Analytical$axes_t_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `Analytical$sync_zyxt_stoch()`

Synchronize states

#### Usage

    Analytical$sync_zyxt_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `Analytical$undo_clear()`

Clear undo list and save current arguments to list

#### Usage

    Analytical$undo_clear()

#### Returns

1

------------------------------------------------------------------------

### `Analytical$undo_save()`

Save current arguments to undo list

#### Usage

    Analytical$undo_save()

#### Returns

number of argument sets

------------------------------------------------------------------------

### `Analytical$undo_undo()`

Replace current arguments from undo list

#### Usage

    Analytical$undo_undo(updn = NULL)

#### Arguments

- `updn`:

  positive to move up, negative to move down

#### Returns

c(index of this argument set, number of argument sets)

------------------------------------------------------------------------

### `Analytical$Drift()`

Calculate, plot and return drifts

#### Usage

    Analytical$Drift(z = NULL, rho = NULL, mu = NULL, who = NULL)

#### Arguments

- `z`:

  vector of n states

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `who`:

  object id of caller

#### Returns

list(g(1xn))

------------------------------------------------------------------------

### `Analytical$Diffusion()`

Calculate, plot and return diffusions

#### Usage

    Analytical$Diffusion(z = NULL, sigma = NULL, who = NULL)

#### Arguments

- `z`:

  vector of n states

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(h2(1xn))

------------------------------------------------------------------------

### `Analytical$Mean()`

Calculate and plot means and time for means to converge

#### Usage

    Analytical$Mean(
      t = NULL,
      s = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      eps = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `eps`:

  proportion remaining after convergence 0\<=eps\<=1

- `who`:

  object id of caller

#### Returns

list(G(mx1),Gteps)

------------------------------------------------------------------------

### `Analytical$Variance()`

Calculate and plot variances and time for variances to converge

#### Usage

    Analytical$Variance(
      t = NULL,
      s = NULL,
      rho = NULL,
      sigma = NULL,
      eps = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `s`:

  initial time -inf\<s\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `eps`:

  proportion remaining after convergence 0\<=eps\<=1

- `who`:

  object id of caller

#### Returns

list(H2(mx1),H2teps)

------------------------------------------------------------------------

### `Analytical$Density()`

Calculate and plot transition densities

#### Usage

    Analytical$Density(
      t = NULL,
      y = NULL,
      s = NULL,
      x = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(p(mxn))

------------------------------------------------------------------------

### `Analytical$Probability()`

Calculate and plot transition probabilities

#### Usage

    Analytical$Probability(
      t = NULL,
      y = NULL,
      s = NULL,
      x = NULL,
      psi = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `y`:

  vector of nstates -inf\<y\<inf

- `s`:

  initial time -inf\<s\<inf

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

- `who`:

  object id of caller

#### Returns

list(P(mxn))

------------------------------------------------------------------------

### `Analytical$DoubleIntegral()`

Calculate and plot double integrals of transition densities

#### Usage

    Analytical$DoubleIntegral(
      t = NULL,
      y = NULL,
      s = NULL,
      x = NULL,
      psi = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `y`:

  vector of n states -inf\<y\<inf

- `s`:

  initial time -inf\<s\<inf

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

- `who`:

  object id of caller

#### Returns

list(PP(mxn))

------------------------------------------------------------------------

### `Analytical$Option()`

Calculate and plot option prices

#### Usage

    Analytical$Option(
      s = NULL,
      x = NULL,
      t = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      b = NULL,
      c = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<t

- `x`:

  vector of n states -inf\<x\<inf

- `t`:

  terminal time -inf\<t\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `b`:

  lump-sum benefit for entry option

- `c`:

  lump-sum cost for exit option

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(OO(mxn))

------------------------------------------------------------------------

### `Analytical$OptionEnvelope()`

Calculate and plot the envelope of option prices

#### Usage

    Analytical$OptionEnvelope(
      x = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      b = NULL,
      c = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states -inf\<x\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `b`:

  lump-sum benefit for entry option

- `c`:

  lump-sum cost for exit option

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(OOhat(1xn),shat(1xn))

------------------------------------------------------------------------

### `Analytical$DecisionThreshold()`

Calculate and plot the decision threshold

#### Usage

    Analytical$DecisionThreshold(
      x = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      b = NULL,
      c = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states -inf\<x\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for exit option, \>0 for entry option

- `b`:

  lump-sum benefit for entry option

- `c`:

  lump-sum cost for exit option

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(k,OO)

------------------------------------------------------------------------

### `Analytical$Obligation()`

Calculate and plot obligations and prohibitions, ie. benefit/cost
analysis

#### Usage

    Analytical$Obligation(
      s = NULL,
      x = NULL,
      t = NULL,
      y = NULL,
      r = NULL,
      phi = NULL,
      b = NULL,
      c = NULL,
      rho = NULL,
      mu = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times -inf\<s\<t

- `x`:

  vector of n states -inf\<x\<inf

- `t`:

  terminal time -inf\<t\<inf

- `y`:

  terminal state -inf\<y\<inf

- `r`:

  discount rate -inf\<r\<inf

- `phi`:

  \<=0 for obligation, \>0 for prohibition

- `b`:

  lump-sum benefit for entry option

- `c`:

  lump-sum cost for exit option

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `who`:

  object id of caller

#### Returns

list(BC(mxn))

------------------------------------------------------------------------

### `Analytical$PassageTimeModeMedianMean()`

Calculate and plot passage time modes, medians and means

#### Usage

    Analytical$PassageTimeModeMedianMean(
      k = NULL,
      s = NULL,
      x = NULL,
      z = NULL,
      omega = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `k`:

  decision threshold -inf\<k\<int

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `z`:

  vector of alternate initial states -inf\<z\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(tmodemedianmean,tmodesmediansmeans(n),)

------------------------------------------------------------------------

### `Analytical$PassageTimePercentiles()`

Calculate and plot passage time percentiles

#### Usage

    Analytical$PassageTimePercentiles(
      k = NULL,
      s = NULL,
      x = NULL,
      z = NULL,
      omega = NULL,
      Ppct = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `k`:

  decision threshold -inf\<k\<int

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `z`:

  vector of alternate initial states -inf\<z\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `Ppct`:

  passage time probability for a percentile 0.01\<=Ppct\<=0.99

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(tpercentile(3),tpercentiles(3xn))

------------------------------------------------------------------------

### `Analytical$PassageTimeDensity()`

Calculate and plot passage time densities

#### Usage

    Analytical$PassageTimeDensity(
      t = NULL,
      k = NULL,
      s = NULL,
      x = NULL,
      z = NULL,
      omega = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `z`:

  vector of alternate initial states -inf\<z\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(ptx(mx1),pt(mxn))

------------------------------------------------------------------------

### `Analytical$PassageTimeProbability()`

Calculate and plot passage time probabilities

#### Usage

    Analytical$PassageTimeProbability(
      t = NULL,
      k = NULL,
      s = NULL,
      x = NULL,
      z = NULL,
      omega = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `t`:

  vector of m times s\<=t\<inf

- `k`:

  decision threshold -inf\<k\<int

- `s`:

  initial time -inf\<s\<inf

- `x`:

  initial state -inf\<x\<inf

- `z`:

  vector of alternate initial states -inf\<z\<inf

- `omega`:

  degree of irreversibility 0\<=omega\<=1

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(Ptx(1xn),Pt(mxn))

------------------------------------------------------------------------

### `Analytical$PlotDrift()`

Plot drifts

#### Usage

    Analytical$PlotDrift(
      type = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = 0

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `Analytical$PlotDiffusion()`

Plot diffusions

#### Usage

    Analytical$PlotDiffusion(
      type = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      zbeg = NULL,
      zend = NULL
    )

#### Arguments

- `type`:

  = -1, 0, or 'n','p','d' for next, previous, default

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `zbeg`:

  begin value for state axis

- `zend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `Analytical$PlotMean()`

Plot means

#### Usage

    Analytical$PlotMean(
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

  = -1,...,2 or 'n','p','d' for next, previous, default

- `pmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotVariance()`

Plot variances

#### Usage

    Analytical$PlotVariance(
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

  = -1,...,2 or 'n','p','d' for next, previous, default

- `pmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotDensity()`

Plot transition densities

#### Usage

    Analytical$PlotDensity(
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

  = 0, 1 or 'n','p','d' for next, previous, default

- `pmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotProbability()`

Plot transition probabilities

#### Usage

    Analytical$PlotProbability(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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

### `Analytical$PlotDoubleIntegral()`

Plot double integrals of transition densities

#### Usage

    Analytical$PlotDoubleIntegral(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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

### `Analytical$PlotOption()`

Plot options

#### Usage

    Analytical$PlotOption(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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

### `Analytical$PlotOptionEnvelope()`

Plot the option envelope

#### Usage

    Analytical$PlotOptionEnvelope(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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

### `Analytical$PlotDecisionThreshold()`

Plot the decision threshold

#### Usage

    Analytical$PlotDecisionThreshold(
      type = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      xbeg = NULL,
      xend = NULL
    )

#### Arguments

- `type`:

  = 0

- `title`:

  text for plot title

- `xaxis`:

  text for x-axis label

- `yaxis`:

  text for y-axis label

- `xbeg`:

  begin value for state axis

- `xend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `Analytical$PlotObligation()`

Plot obligations

#### Usage

    Analytical$PlotObligation(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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

### `Analytical$PlotPassageTimeModeMedianMean()`

Plot passage time modes, medians and means

#### Usage

    Analytical$PlotPassageTimeModeMedianMean(
      type = NULL,
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

  = -3,...,2 or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotPassageTimePercentiles()`

Plot passage time lower, median and upper percentiles

#### Usage

    Analytical$PlotPassageTimePercentiles(
      type = NULL,
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

  = -3,...,2 or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotPassageTimeDensity()`

Plot passage time densities

#### Usage

    Analytical$PlotPassageTimeDensity(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

- `ptmax`:

  maximum scale for vertical axis

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

### `Analytical$PlotPassageTimeProbability()`

Plot passage time probabilities

#### Usage

    Analytical$PlotPassageTimeProbability(
      type = NULL,
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

  = 0, 1 or 'n','p','d' for next, previous, default

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
