# R6 class implementing a Finite Difference method for the Ornstein-Uhlenbeck Process.

A Finite Difference method–without boundary conditions on the state but
with arbitrary terminal values–for the pricing of index insurance,
perpetual options and sequential options. From the option prices, an
option envelope and a decision threshold can be calculated. Drift,
diffusion and several possible terminal values are pre-programmed.
User-defined terminal values can also be entered.

## Formulas and Methods:

    x stochastic
      Drift
      Diffusion
      TerminalValue_Linear
      TerminalValue_Degenerate
      TerminalValue_Stepped
      TerminalValue_Kinked
      TerminalValue_Butterfly
      TerminalValue_Mitscherlich
      TerminalValue_Gompertz
      TerminalValue_Logistic
      TerminalValue_Transcendental
      TerminalValue_YieldIndex
      TerminalValue
      Option
      OptionEnvelope
      DecisionThreshold

## Plots

      PlotDrift
      PlotDiffusion
      PlotTerminalValue
      PlotOption
      PlotOptionEnvelope
      PlotDecisionThreshold

## Arguments of functions:

      All arguments are optional in all functions.
      s:     vector of times
      x:     vector of states
      V:     vector of terminal values
      r:     discount rate
      phi:   search direction for exit or entry options
      theta: weight of current time in time stepping
      skip:  divide the time interval and report every skip result
      rho:   rate parameter
      mu:    location parameter
      sigma: scale parameter
      xo:    state at the intercept, spike, step or kink
      xi:    state at the inflection point
      xm:    state at the maximum or kink
      vs:    slope or direction of a step
      vr:    rate of change
      Vmax:  maximum terminal value
      Vmin:  minimum terminal value

## Usage:

The FiniteDifference object must first be instantiated before its
methods are called. There are two ways. The first way instantiates the
OUProcess object and calls a function to get a pointer:

      OUP <- OUProcess$new()
      A <- OUP$get_Analytical()
      FD <- OUP$get_FiniteDifference()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()

The FiniteDifference object will coordinate arguments to functions with
the Analytial, MaximumLikelihood and MonteCarlo objects. The second way
instantiates the FiniteDifference object by itself with no coordination:

      FD <- FiniteDifference$new()

Once the object is instantiated, its methods can be called, to calculate
and plot a Decision Threshold, for example:

      FD$DecisionThreshold()

The plot methods can be used to customize the plots, with a title, for
example:

      FD$PlotDecisionThreshold(title="My Decision")

Other functions and methods are called in the same way. To see all the
possibilities, check out the demos below.

## Demos:

Demonstration scripts are in files in the 'demo' directory. The most
convenient way to list and run demos are the commands:

      OUPDemoList()
      OUPDemoRun(<number of demo in list>)

Entering the demos by number in the list saves typing.

## Discussion:

The Finite Difference Method can solve problems with no analytical
solutions. As a numerical procedure, it will have errors. The larger the
scale parameter, the larger the errors. To manage the errors, the domain
for x should be as large as practical. The increment between times is
ds. Parameter skip subdivides ds into smaller intervals. Only
calculations for ds are reported. The increment between states is dx and
ds/skip should be about dx/100. Increasing skip or dx may reduce errors
up to a point. If possible, the Finite Difference Method should be
calibrated against an analytical solution for a simpler problem before
solving more complicated problems.

The functions include many different types of terminal values, but
suppose the bottom line is a simple exit option:

      FD <- FiniteDifference$new()
      FD$DecisionThreshold(mu=10)

Or a simple entry option:

      FD$TerminalValue_Kinked(xo=5,vs=1)
      FD$DecisionThreshold()

(That's x oh, not x zero). And finally a sequence of an entry with an
option to exit, calculating the exit decision first:

      FD$TerminalValue_Kinked(xo=0,vs=-1)
      exit <- FD$OptionEnvelope()$Ohat
      V <- exit+FD$TerminalValue_Kinked(xo=5,vs=1)$V
      FD$DecisionThreshold(V=V,phi=1)

The functions all return named lists. Before being used in calculations,
a list must be stripped of its names. In the above example, the names
are Ohat for the option envelope and V for the terminal value. The
functions try to help by accepting named lists as arguments. For
example:

      V <- FD$TerminalValue_Kinked(xo=5,vs=1)
      FD$DecisionThreshold(V=V)

But the terminal value functions automatically set the terminal value in
the object, so entering V as as argument does nothing useful in this
case.

You may wish to calculate your own terminal value in the console. You
must first get rid of the names. The names to get rid of are listed
under 'Returns' in the help for each function. You can also use double
brackets to eliminate names if you know the position in the list. Both
'OOhat' and 'V' are first in their lists:

      FD$TerminalValue_Kinked(xo=0,vs=-1)
      exit <- FD$OptionEnvelope()[[1]]
      V <- exit+FD$TerminalValue_Kinked(xo=5,vs=1)[[1]]
      FD$DecisionThreshold(V=V,phi=1)

Or you can unlist everything:

      V <- unlist(FD$OptionEnvelope())+unlist(FD$TerminalValue_Kinked(xo=5,vs=1))

Sorry.

## Methods

### Public methods

- [`FiniteDifference$new()`](#method-FiniteDifference-initialize)

- [`FiniteDifference$set_oup_params()`](#method-FiniteDifference-set_oup_params)

- [`FiniteDifference$set_x_stoch_args()`](#method-FiniteDifference-set_x_stoch_args)

- [`FiniteDifference$set_V_linear_args()`](#method-FiniteDifference-set_V_linear_args)

- [`FiniteDifference$set_V_degenerate_args()`](#method-FiniteDifference-set_V_degenerate_args)

- [`FiniteDifference$set_V_stepped_args()`](#method-FiniteDifference-set_V_stepped_args)

- [`FiniteDifference$set_V_kinked_args()`](#method-FiniteDifference-set_V_kinked_args)

- [`FiniteDifference$set_V_butterfly_args()`](#method-FiniteDifference-set_V_butterfly_args)

- [`FiniteDifference$set_V_mitscherlich_args()`](#method-FiniteDifference-set_V_mitscherlich_args)

- [`FiniteDifference$set_V_gompertz_args()`](#method-FiniteDifference-set_V_gompertz_args)

- [`FiniteDifference$set_V_logistic_args()`](#method-FiniteDifference-set_V_logistic_args)

- [`FiniteDifference$set_V_transcendental_args()`](#method-FiniteDifference-set_V_transcendental_args)

- [`FiniteDifference$set_V_yieldindex_args()`](#method-FiniteDifference-set_V_yieldindex_args)

- [`FiniteDifference$set_V_args()`](#method-FiniteDifference-set_V_args)

- [`FiniteDifference$set_V_info()`](#method-FiniteDifference-set_V_info)

- [`FiniteDifference$set_plot_info()`](#method-FiniteDifference-set_plot_info)

- [`FiniteDifference$set_plot_type()`](#method-FiniteDifference-set_plot_type)

- [`FiniteDifference$set_flags()`](#method-FiniteDifference-set_flags)

- [`FiniteDifference$get_all()`](#method-FiniteDifference-get_all)

- [`FiniteDifference$get_oup_params()`](#method-FiniteDifference-get_oup_params)

- [`FiniteDifference$get_x_stoch_args()`](#method-FiniteDifference-get_x_stoch_args)

- [`FiniteDifference$get_V_linear_args()`](#method-FiniteDifference-get_V_linear_args)

- [`FiniteDifference$get_V_degenerate_args()`](#method-FiniteDifference-get_V_degenerate_args)

- [`FiniteDifference$get_V_stepped_args()`](#method-FiniteDifference-get_V_stepped_args)

- [`FiniteDifference$get_V_kinked_args()`](#method-FiniteDifference-get_V_kinked_args)

- [`FiniteDifference$get_V_butterfly_args()`](#method-FiniteDifference-get_V_butterfly_args)

- [`FiniteDifference$get_V_mitscherlich_args()`](#method-FiniteDifference-get_V_mitscherlich_args)

- [`FiniteDifference$get_V_gompertz_args()`](#method-FiniteDifference-get_V_gompertz_args)

- [`FiniteDifference$get_V_logistic_args()`](#method-FiniteDifference-get_V_logistic_args)

- [`FiniteDifference$get_V_transcendental_args()`](#method-FiniteDifference-get_V_transcendental_args)

- [`FiniteDifference$get_V_yieldindex_args()`](#method-FiniteDifference-get_V_yieldindex_args)

- [`FiniteDifference$get_V_args()`](#method-FiniteDifference-get_V_args)

- [`FiniteDifference$get_V_info()`](#method-FiniteDifference-get_V_info)

- [`FiniteDifference$get_plot_info()`](#method-FiniteDifference-get_plot_info)

- [`FiniteDifference$get_plot_colors()`](#method-FiniteDifference-get_plot_colors)

- [`FiniteDifference$get_plot_types()`](#method-FiniteDifference-get_plot_types)

- [`FiniteDifference$get_flags()`](#method-FiniteDifference-get_flags)

- [`FiniteDifference$axes_x_stoch()`](#method-FiniteDifference-axes_x_stoch)

- [`FiniteDifference$undo_clear()`](#method-FiniteDifference-undo_clear)

- [`FiniteDifference$undo_save()`](#method-FiniteDifference-undo_save)

- [`FiniteDifference$undo_undo()`](#method-FiniteDifference-undo_undo)

- [`FiniteDifference$Drift()`](#method-FiniteDifference-Drift)

- [`FiniteDifference$Diffusion()`](#method-FiniteDifference-Diffusion)

- [`FiniteDifference$TerminalValue_Linear()`](#method-FiniteDifference-TerminalValue_Linear)

- [`FiniteDifference$TerminalValue_Degenerate()`](#method-FiniteDifference-TerminalValue_Degenerate)

- [`FiniteDifference$TerminalValue_Stepped()`](#method-FiniteDifference-TerminalValue_Stepped)

- [`FiniteDifference$TerminalValue_Kinked()`](#method-FiniteDifference-TerminalValue_Kinked)

- [`FiniteDifference$TerminalValue_Butterfly()`](#method-FiniteDifference-TerminalValue_Butterfly)

- [`FiniteDifference$TerminalValue_Mitscherlich()`](#method-FiniteDifference-TerminalValue_Mitscherlich)

- [`FiniteDifference$TerminalValue_Gompertz()`](#method-FiniteDifference-TerminalValue_Gompertz)

- [`FiniteDifference$TerminalValue_Logistic()`](#method-FiniteDifference-TerminalValue_Logistic)

- [`FiniteDifference$TerminalValue_Transcendental()`](#method-FiniteDifference-TerminalValue_Transcendental)

- [`FiniteDifference$TerminalValue_YieldIndex()`](#method-FiniteDifference-TerminalValue_YieldIndex)

- [`FiniteDifference$TerminalValue()`](#method-FiniteDifference-TerminalValue)

- [`FiniteDifference$Option()`](#method-FiniteDifference-Option)

- [`FiniteDifference$OptionEnvelope()`](#method-FiniteDifference-OptionEnvelope)

- [`FiniteDifference$DecisionThreshold()`](#method-FiniteDifference-DecisionThreshold)

- [`FiniteDifference$PlotDrift()`](#method-FiniteDifference-PlotDrift)

- [`FiniteDifference$PlotDiffusion()`](#method-FiniteDifference-PlotDiffusion)

- [`FiniteDifference$PlotTerminalValue()`](#method-FiniteDifference-PlotTerminalValue)

- [`FiniteDifference$PlotOption()`](#method-FiniteDifference-PlotOption)

- [`FiniteDifference$PlotOptionEnvelope()`](#method-FiniteDifference-PlotOptionEnvelope)

- [`FiniteDifference$PlotDecisionThreshold()`](#method-FiniteDifference-PlotDecisionThreshold)

------------------------------------------------------------------------

### `FiniteDifference$new()`

Create a new FiniteDifference object

#### Usage

    FiniteDifference$new(OUP = NULL)

#### Arguments

- `OUP`:

  pointer set by the OUProcess object

#### Returns

A new FiniteDifference object

------------------------------------------------------------------------

### `FiniteDifference$set_oup_params()`

Set OUP parameters

#### Usage

    FiniteDifference$set_oup_params(
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

list(rho,mu,sigma)

------------------------------------------------------------------------

### `FiniteDifference$set_x_stoch_args()`

Set x as a stochastic state and its arguments

#### Usage

    FiniteDifference$set_x_stoch_args(
      s = NULL,
      x = NULL,
      V = NULL,
      r = NULL,
      phi = NULL,
      theta = NULL,
      skip = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of times

- `x`:

  vector of states

- `V`:

  vector of terminal values

- `r`:

  discount rate 0\<=r\<inf, scalar

- `phi`:

  search direction for exit or entry options

- `theta`:

  weight of current time in time stepping 0.5\<=theta\<=1

- `skip`:

  subdivide time intervals but report at times s 1\<=skip\<=1000

- `who`:

  object id of sender

#### Returns

list(s,x,V,r,phi,theta,skip,ds,dx)

------------------------------------------------------------------------

### `FiniteDifference$set_V_linear_args()`

Set V linear arguments

#### Usage

    FiniteDifference$set_V_linear_args(xo = NULL, vs = NULL)

#### Arguments

- `xo`:

  state at the intercept

- `vs`:

  slope

#### Returns

list(xo,vs)

------------------------------------------------------------------------

### `FiniteDifference$set_V_degenerate_args()`

Set V degenerate arguments

#### Usage

    FiniteDifference$set_V_degenerate_args(xo = NULL, Vmax = NULL, Vmin = NULL)

#### Arguments

- `xo`:

  state at the step

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_stepped_args()`

Set V stepped arguments

#### Usage

    FiniteDifference$set_V_stepped_args(
      xo = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the step

- `vs`:

  direction of the step

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_kinked_args()`

Set V kinked arguments

#### Usage

    FiniteDifference$set_V_kinked_args(
      xo = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the kink

- `vs`:

  slope

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_butterfly_args()`

Set V butterfly arguments

#### Usage

    FiniteDifference$set_V_butterfly_args(
      xo = NULL,
      xm = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the left wing

- `xm`:

  state at the right wing

- `vs`:

  slope

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,xm,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_mitscherlich_args()`

Set V mitscherlich arguments

#### Usage

    FiniteDifference$set_V_mitscherlich_args(
      xo = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the intercept

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_gompertz_args()`

Set V gompertz arguments

#### Usage

    FiniteDifference$set_V_gompertz_args(
      xi = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xi`:

  state at the inflection point

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xi,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_logistic_args()`

Set V logistic arguments

#### Usage

    FiniteDifference$set_V_logistic_args(
      xi = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xi`:

  state at the inflection point

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xi,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_transcendental_args()`

Set V transcendental arguments

#### Usage

    FiniteDifference$set_V_transcendental_args(
      xo = NULL,
      xi = NULL,
      xm = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the intercept

- `xi`:

  state at the inflection point

- `xm`:

  state at the maximum

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,xi,xm,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_yieldindex_args()`

Set V yield index arguments

#### Usage

    FiniteDifference$set_V_yieldindex_args(
      xo = NULL,
      xi = NULL,
      xm = NULL,
      Vmax = NULL,
      Vmin = NULL
    )

#### Arguments

- `xo`:

  state at the intercept

- `xi`:

  state at the inflection point

- `xm`:

  state at the maximum

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

#### Returns

list(xo,xi,xm,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$set_V_args()`

Set V arguments by position in list

#### Usage

    FiniteDifference$set_V_args(
      Ix = NULL,
      name = NULL,
      v1 = NULL,
      v2 = NULL,
      v3 = NULL,
      v4 = NULL,
      v5 = NULL
    )

#### Arguments

- `Ix`:

  index of terminal value

- `name`:

  name of terminal value

- `v1`:

  first parameter

- `v2`:

  second parameter

- `v3`:

  third parameter

- `v4`:

  fourth parameter

- `v5`:

  fifth parameter

#### Returns

list(Ix,name,v1,v2,v3,v4,v5)

------------------------------------------------------------------------

### `FiniteDifference$set_V_info()`

Set information for terminal value

#### Usage

    FiniteDifference$set_V_info(Ix = NULL, name = NULL)

#### Arguments

- `Ix`:

  index of terminal value

- `name`:

  name of terminal value

#### Returns

list(Ix,name,names,text)

------------------------------------------------------------------------

### `FiniteDifference$set_plot_info()`

Set information for plotting

#### Usage

    FiniteDifference$set_plot_info(
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

### `FiniteDifference$set_plot_type()`

Set plot types by group

#### Usage

    FiniteDifference$set_plot_type(type = NULL, group = NULL)

#### Arguments

- `type`:

  a number identifying the type for a group

- `group`:

  identifier for groups of similar plots

#### Returns

list(types,group)

------------------------------------------------------------------------

### `FiniteDifference$set_flags()`

Set flags for plotting and copying

#### Usage

    FiniteDifference$set_flags(plotit = NULL, copyit = NULL, who = NULL)

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

### `FiniteDifference$get_all()`

Get all arguments and information

#### Usage

    FiniteDifference$get_all()

#### Returns

list(oup_params,x_stoch_args,V_linear_args,V_stepped_args,V_kinked_args,V_butterfly_args,V_mitscherlich_args,V_gompertz_args,V_logistic_args,V_transcendental_args,V_yieldindex_args,plot_info)

------------------------------------------------------------------------

### `FiniteDifference$get_oup_params()`

Get OUP parameters

#### Usage

    FiniteDifference$get_oup_params()

#### Returns

list(rho,mu,sigma)

------------------------------------------------------------------------

### `FiniteDifference$get_x_stoch_args()`

Get x as a stochastic state and its arguments

#### Usage

    FiniteDifference$get_x_stoch_args()

#### Returns

list(s,x,V,r,phi,theta,skip,ds,dx)

------------------------------------------------------------------------

### `FiniteDifference$get_V_linear_args()`

Get V linear arguments

#### Usage

    FiniteDifference$get_V_linear_args()

#### Returns

list(xo,vs)

------------------------------------------------------------------------

### `FiniteDifference$get_V_degenerate_args()`

Get V degenerate arguments

#### Usage

    FiniteDifference$get_V_degenerate_args()

#### Returns

list(xo,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_stepped_args()`

Get V stepped arguments

#### Usage

    FiniteDifference$get_V_stepped_args()

#### Returns

list(xo,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_kinked_args()`

Get V kinked arguments

#### Usage

    FiniteDifference$get_V_kinked_args()

#### Returns

list(xo,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_butterfly_args()`

Get V butterfly arguments

#### Usage

    FiniteDifference$get_V_butterfly_args()

#### Returns

list(xo,xm,vs,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_mitscherlich_args()`

Get V mitscherlich arguments

#### Usage

    FiniteDifference$get_V_mitscherlich_args()

#### Returns

list(xo,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_gompertz_args()`

Get V gompertz arguments

#### Usage

    FiniteDifference$get_V_gompertz_args()

#### Returns

list(xi,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_logistic_args()`

Get V logistic arguments

#### Usage

    FiniteDifference$get_V_logistic_args()

#### Returns

list(xi,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_transcendental_args()`

Get V transcendental arguments

#### Usage

    FiniteDifference$get_V_transcendental_args()

#### Returns

list(xo,xi,xm,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_yieldindex_args()`

Get V yield index arguments

#### Usage

    FiniteDifference$get_V_yieldindex_args()

#### Returns

list(xo,xi,xm,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_args()`

Get V arguments by index or name

#### Usage

    FiniteDifference$get_V_args(Ix = NULL, name = NULL)

#### Arguments

- `Ix`:

  index of terminal value

- `name`:

  name of terminal value

#### Returns

list(xo,xi,xm,vs,vr,Vmax,Vmin)

------------------------------------------------------------------------

### `FiniteDifference$get_V_info()`

Get information for terminal values

#### Usage

    FiniteDifference$get_V_info()

#### Returns

list(Ix,name,names,text)

------------------------------------------------------------------------

### `FiniteDifference$get_plot_info()`

Get information for plotting options

#### Usage

    FiniteDifference$get_plot_info()

#### Returns

list(font,file,theme,3D,labels)

------------------------------------------------------------------------

### `FiniteDifference$get_plot_colors()`

Get colors for plotting

#### Usage

    FiniteDifference$get_plot_colors()

#### Returns

list(red,ylw,grn,cyn,blu,mgn,gry,background,font,reverse)

------------------------------------------------------------------------

### `FiniteDifference$get_plot_types()`

Get current types for plot routines

#### Usage

    FiniteDifference$get_plot_types()

#### Returns

list(types,description)

------------------------------------------------------------------------

### `FiniteDifference$get_flags()`

Get flags for plotting and copying

#### Usage

    FiniteDifference$get_flags()

#### Returns

list(plotit,copyit)

------------------------------------------------------------------------

### `FiniteDifference$axes_x_stoch()`

Scale axes for x stochastic arguments

#### Usage

    FiniteDifference$axes_x_stoch()

#### Returns

NULL

------------------------------------------------------------------------

### `FiniteDifference$undo_clear()`

Clear undo list and save current arguments to list

#### Usage

    FiniteDifference$undo_clear()

#### Returns

1

------------------------------------------------------------------------

### `FiniteDifference$undo_save()`

Save current arguments to undo list

#### Usage

    FiniteDifference$undo_save()

#### Returns

number of argument sets

------------------------------------------------------------------------

### `FiniteDifference$undo_undo()`

Replace current arguments from undo list

#### Usage

    FiniteDifference$undo_undo(updn = NULL)

#### Arguments

- `updn`:

  positive to move up, negative to move down

#### Returns

c(index of this argument set, number of argument sets)

------------------------------------------------------------------------

### `FiniteDifference$Drift()`

Calculate, plot and return drifts

#### Usage

    FiniteDifference$Drift(x = NULL, rho = NULL, mu = NULL, who = NULL)

#### Arguments

- `x`:

  vector of n states

- `rho`:

  rate parameter 0\<=\<rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `who`:

  object id of caller

#### Returns

list(g(1xn))

------------------------------------------------------------------------

### `FiniteDifference$Diffusion()`

Calculate, plot and return diffusions

#### Usage

    FiniteDifference$Diffusion(x = NULL, sigma = NULL, who = NULL)

#### Arguments

- `x`:

  vector of n states

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(h2(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Linear()`

Create and plot linear terminal values

#### Usage

    FiniteDifference$TerminalValue_Linear(
      x = NULL,
      xo = NULL,
      vs = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the intercept

- `vs`:

  slope

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Degenerate()`

Create and plot degenerate terminal values

#### Usage

    FiniteDifference$TerminalValue_Degenerate(
      x = NULL,
      xo = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the spike

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Stepped()`

Create and plot stepped terminal values

#### Usage

    FiniteDifference$TerminalValue_Stepped(
      x = NULL,
      xo = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the step

- `vs`:

  direction of the step

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Kinked()`

Create and plot kinked terminal values

#### Usage

    FiniteDifference$TerminalValue_Kinked(
      x = NULL,
      xo = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the kink

- `vs`:

  slope

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Butterfly()`

Create and plot butterfly terminal values

#### Usage

    FiniteDifference$TerminalValue_Butterfly(
      x = NULL,
      xo = NULL,
      xm = NULL,
      vs = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the left wing

- `xm`:

  state at the right wing

- `vs`:

  slope

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Mitscherlich()`

Create and plot Mitscherlich terminal values

#### Usage

    FiniteDifference$TerminalValue_Mitscherlich(
      x = NULL,
      xo = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the intercept

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Gompertz()`

Create and plot Gompertz terminal values

#### Usage

    FiniteDifference$TerminalValue_Gompertz(
      x = NULL,
      xi = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xi`:

  state at the inflection point

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Logistic()`

Create and plot Logistic terminal values

#### Usage

    FiniteDifference$TerminalValue_Logistic(
      x = NULL,
      xi = NULL,
      vr = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xi`:

  state at the inflection point

- `vr`:

  rate of change

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_Transcendental()`

Create and plot Transcendental terminal values

#### Usage

    FiniteDifference$TerminalValue_Transcendental(
      x = NULL,
      xo = NULL,
      xi = NULL,
      xm = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the intercept

- `xi`:

  state at the inflection point

- `xm`:

  state at the maximum

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue_YieldIndex()`

Create and plot Yield Index terminal values

#### Usage

    FiniteDifference$TerminalValue_YieldIndex(
      x = NULL,
      xo = NULL,
      xi = NULL,
      xm = NULL,
      Vmax = NULL,
      Vmin = NULL,
      who = NULL
    )

#### Arguments

- `x`:

  vector of n states

- `xo`:

  state at the intercept

- `xi`:

  state at the inflection point

- `xm`:

  state at the maximum

- `Vmax`:

  maximum terminal value

- `Vmin`:

  minimum terminal value

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$TerminalValue()`

Retrieves and plots terminal values

#### Usage

    FiniteDifference$TerminalValue(Ix = NULL, name = NULL, who = NULL)

#### Arguments

- `Ix`:

  index number or name of terminal values

- `name`:

  name of terminal values

- `who`:

  object id of caller

#### Returns

list(V(1xn))

------------------------------------------------------------------------

### `FiniteDifference$Option()`

Calculate and plot option prices

#### Usage

    FiniteDifference$Option(
      s = NULL,
      x = NULL,
      V = NULL,
      r = NULL,
      theta = NULL,
      skip = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times

- `x`:

  vector of n states

- `V`:

  vector of n terminal values

- `r`:

  discount rate 0\<=r\<inf

- `theta`:

  weight of current time in time stepping 0.5\<=theta\<=1

- `skip`:

  subdivide time intervals but report at times s 1\<=skip\<=1000

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

### `FiniteDifference$OptionEnvelope()`

Calculate and plot the envelope of option prices

#### Usage

    FiniteDifference$OptionEnvelope(
      s = NULL,
      x = NULL,
      V = NULL,
      r = NULL,
      theta = NULL,
      skip = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times

- `x`:

  vector of n states

- `V`:

  vector of n terminal values

- `r`:

  discount rate 0\<=r\<inf

- `theta`:

  weight of current time in time stepping 0.5\<=theta\<=1

- `skip`:

  subdivide time intervals but report at times s 1\<=skip\<=1000

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(OOhat(1xn))

------------------------------------------------------------------------

### `FiniteDifference$DecisionThreshold()`

Calculate and plot the decision threshold

#### Usage

    FiniteDifference$DecisionThreshold(
      s = NULL,
      x = NULL,
      V = NULL,
      r = NULL,
      phi = NULL,
      theta = NULL,
      skip = NULL,
      rho = NULL,
      mu = NULL,
      sigma = NULL,
      who = NULL
    )

#### Arguments

- `s`:

  vector of m times

- `x`:

  vector of n states

- `V`:

  vector of n terminal values

- `r`:

  discount rate 0\<=r\<inf

- `phi`:

  search direction for exit or entry options

- `theta`:

  weight of current time in time stepping 0.5\<=theta\<=1

- `skip`:

  subdivide time intervals but report at times s 1\<=skip\<=1000

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of caller

#### Returns

list(k,OOhat)

------------------------------------------------------------------------

### `FiniteDifference$PlotDrift()`

Plot drifts

#### Usage

    FiniteDifference$PlotDrift(
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

### `FiniteDifference$PlotDiffusion()`

Plot diffusions

#### Usage

    FiniteDifference$PlotDiffusion(
      type = NULL,
      title = NULL,
      xaxis = NULL,
      yaxis = NULL,
      xbeg = NULL,
      xend = NULL
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

- `xbeg`:

  begin value for state axis

- `xend`:

  end value for state axis

#### Returns

plot

------------------------------------------------------------------------

### `FiniteDifference$PlotTerminalValue()`

Plot terminal values

#### Usage

    FiniteDifference$PlotTerminalValue(
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

### `FiniteDifference$PlotOption()`

Plot options

#### Usage

    FiniteDifference$PlotOption(
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

  = 0, 1, or 'n','p','d' for next, previous, default

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

### `FiniteDifference$PlotOptionEnvelope()`

Plot the option envelope

#### Usage

    FiniteDifference$PlotOptionEnvelope(
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

  = 0, 1, or 'n','p','d' for next, previous, default

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

### `FiniteDifference$PlotDecisionThreshold()`

Plot the decision threshold

#### Usage

    FiniteDifference$PlotDecisionThreshold(
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
