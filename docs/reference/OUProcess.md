# R6 class implementing a container to synchronize OUP objects.

The Analytical, FiniteDifference, MaximumLikelihood and MonteCarlo
objects are a complete set of functions for maximum likelihood
estimation and for the calculation of probabilities, option prices,
visiting times, first passage times, decision thresholds and
more–everything for a Real Options Analysis. Each object can be used by
itself. This object, OUProcess, will instantiate and synchronize the
other objects.

## Usage:

This object, OUProcess, and the Analytical, FiniteDifference,
MaximumLikelihood and MonteCarlo objects can be instantiated together
by:

      OUP <- OUProcess$new()

Then pointers to individual objects can be accessed by:

      A <- OUP$get_Analytical()
      FD <- OUP$get_FiniteDifference()
      ML <- OUP$get_MaximumLikelihood()
      MC <- OUP$get_MonteCarlo()

OUP has other public methods, but these are called by individual objects
and won't do anything if called by a user.

## Discussion

OUP will have pointers to A, FD, ML and MC. A, FD, ML and MC will each
have a pointer to OUP. Thus A, FD, ML and MC can follow the pointers to
synchronize with each other. For example, rho, mu and sigma estimated in
ML will be propagated to the other objects:

      ML$Estimates()
      A$Option()

Alternately, each object can be instantiated individually:

      A <- Analytical$new()
      FD <- FiniteDifference$new()
      ML <- MaximumLikelihood$new()
      MC <- MonteCarlo$new()

In this case, A, FD, ML and MC will not have pointers and will not
synchronize. To estimate rho, mu and sigma and calculate option prices
would require:

      estimates <- ML$Estimates()
      rho <- estimates$rho
      mu <- estimates$mu
      sigma <- estimates$sigma
      A$Option(rho=rho,mu=mu,sigma=sigma)

Upon initialization, OUP will ask if it is being instantiated by an
RShiny app. If so, it will store the session pointer and use it to
coordinate clipboard access between RShiny and the A, FD, ML and MC
objects.

## Methods

### Public methods

- [`OUProcess$new()`](#method-OUProcess-initialize)

- [`OUProcess$get_Analytical()`](#method-OUProcess-get_Analytical)

- [`OUProcess$get_FiniteDifference()`](#method-OUProcess-get_FiniteDifference)

- [`OUProcess$get_MaximumLikelihood()`](#method-OUProcess-get_MaximumLikelihood)

- [`OUProcess$get_MonteCarlo()`](#method-OUProcess-get_MonteCarlo)

- [`OUProcess$send_oup_params()`](#method-OUProcess-send_oup_params)

- [`OUProcess$send_y_stoch_args()`](#method-OUProcess-send_y_stoch_args)

- [`OUProcess$send_x_stoch_args()`](#method-OUProcess-send_x_stoch_args)

- [`OUProcess$send_t_stoch_args()`](#method-OUProcess-send_t_stoch_args)

- [`OUProcess$send_plot_args()`](#method-OUProcess-send_plot_args)

- [`OUProcess$send_plot_info()`](#method-OUProcess-send_plot_info)

- [`OUProcess$send_flags()`](#method-OUProcess-send_flags)

- [`OUProcess$CopyToClipboard()`](#method-OUProcess-CopyToClipboard)

------------------------------------------------------------------------

### `OUProcess$new()`

Create an OUProcess object as a container for other objects

#### Usage

    OUProcess$new(session = NULL)

#### Arguments

- `session`:

  session pointer to RShiny app

#### Returns

A new OUProcess object

------------------------------------------------------------------------

### `OUProcess$get_Analytical()`

Pointer to Analytical object

#### Usage

    OUProcess$get_Analytical()

#### Returns

pointer

------------------------------------------------------------------------

### `OUProcess$get_FiniteDifference()`

Pointer to FiniteDifference object

#### Usage

    OUProcess$get_FiniteDifference()

#### Returns

pointer

------------------------------------------------------------------------

### `OUProcess$get_MaximumLikelihood()`

Pointer to MaximumLikelihood object

#### Usage

    OUProcess$get_MaximumLikelihood()

#### Returns

pointer

------------------------------------------------------------------------

### `OUProcess$get_MonteCarlo()`

Pointer to MonteCarlo object

#### Usage

    OUProcess$get_MonteCarlo()

#### Returns

pointer

------------------------------------------------------------------------

### `OUProcess$send_oup_params()`

Function called by modules to coordinate parameters

#### Usage

    OUProcess$send_oup_params(rho = NULL, mu = NULL, sigma = NULL, who = NULL)

#### Arguments

- `rho`:

  rate parameter 0\<=rho\<inf

- `mu`:

  location parameter -inf\<mu\<inf

- `sigma`:

  scale parameter -inf\<sigma\<inf

- `who`:

  object id of sender

------------------------------------------------------------------------

### `OUProcess$send_y_stoch_args()`

Function called by modules to coordinate arguments for y as a stochastic
state

#### Usage

    OUProcess$send_y_stoch_args(
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

------------------------------------------------------------------------

### `OUProcess$send_x_stoch_args()`

Function called by modules to coordinate x as a stochastic state and its
arguments

#### Usage

    OUProcess$send_x_stoch_args(
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

------------------------------------------------------------------------

### `OUProcess$send_t_stoch_args()`

Function called by modules to coordinate t stochastic arguments

#### Usage

    OUProcess$send_t_stoch_args(
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

------------------------------------------------------------------------

### `OUProcess$send_plot_args()`

Function called by modules to coordinate arguments for plotting

#### Usage

    OUProcess$send_plot_args(pmax = NULL, ptmax = NULL, who = NULL)

#### Arguments

- `pmax`:

  maximum transition density

- `ptmax`:

  maximum visiting time and first passage time densities

- `who`:

  object id of sender

------------------------------------------------------------------------

### `OUProcess$send_plot_info()`

Function called by modules to coordinate information for plotting

#### Usage

    OUProcess$send_plot_info(
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

------------------------------------------------------------------------

### `OUProcess$send_flags()`

Function called by modules to coordinate flags for plotting and copying

#### Usage

    OUProcess$send_flags(plotit = NULL, copyit = NULL, who = NULL)

#### Arguments

- `plotit`:

  automatic plot after calculation TRUE or FALSE

- `copyit`:

  copy to clipboard TRUE or FALSE

- `who`:

  object id of sender

------------------------------------------------------------------------

### `OUProcess$CopyToClipboard()`

Write to clipboard

#### Usage

    OUProcess$CopyToClipboard(clip)

#### Arguments

- `clip`:

  data frame, matrix or vector
