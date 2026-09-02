# GregsOUPR6 package implementing Real Options for Adoption and Resilience.

Global functions to launch an RShiny app, read data and run demos.

## Usage

``` r
OUPShiny()

OUPHelpList()

OUPHelpView(help = "General")

OUPDataList()

OUPDataHelp(help = "MyData")

OUPDataRead(file = "MyData")

OUPDemoList()

OUPDemoRun(demo = "A_Drift")
```

## Arguments

- help:

  name of help file or index in file list

- file:

  name of data file or index in file list

- demo:

  name of demo or index in demo list

## Overview

Many of us dread Matrix Management. It means each of us has lots of
bosses. GregsOUPR6 also has a matrix structure and can be bossed around
in many ways. From top to bottom are levels of software abstraction.
From left to right are solution methods. The same problem can be solved
at different levels by different methods.

This document describes the structure and general use of the package.
Each level of software and each solution method has its own help
documents. There are also help documents for data sets.

There are also functions to get help, view data and run demos.

## Software levels

From top to bottom, the software levels are an RShiny app, R6 Objects
and Rcpp/RcppParallel functions. A user can solve a problem from any
level:

    +-----------------------+
    |  RShiny app           |
    |    user interface     |
    |    server             |\
    +-----------------------+ \
               |               \
    +-----------------------+   \
    |  R6 Objects           |    \
    |    set/get inputs     |     \
    |    set/get outputs    |------user
    |    call calculations  |     /
    |    create plots       |    /
    +-----------------------+   /
               |               /
    +-----------------------+ /
    |  Rcpp / RcppParallel  |/
    |    functions          |
    +-----------------------+

User interactions with the RShiny app and the Rcpp/RcppParallel
functions are simple. Either run the app or call the functions. User
interactions with the R6 Objects are more complicated.

The R6 Objects are reactive. In other words, they calculate if the user
asks but don't calculate again unless an input changes and the user asks
again. There is a cache of inputs and a cache of outputs, managed by
set/get functions. In the diagram to follow, A are arguments, I are
inputs, O are outputs, S is a set function and G is a get function.

     set/get inputs:

                    A
      |    cache    |    cache
      G <--  I  <-> S -->  O
      |             |
      I             I

Arguments enter the set function which validates the arguments and
compares the arguments to the input cache. If there is a change, the
input cache is set and all outputs in the output cache which depend on
that argument are set to NULL. Finally, the set function returns the new
input cache. The get function will return the input cache without
checking or changes.

Once the arguments are set, the user can call other functions. In the
diagram below, C is a calculation function called without arguments.

     set/get outputs:

            |    cache
      G --> C <->  O
            |
            O

The calculation function gets the input cache and checks the output
cache. If the output cache is NULL, functions in Rcpp/RcppParallel are
called for the calculations, the cache is set and then returned. If the
output cache is not NULL, it is returned. Thus a calculation function
calculates once and then acts like a get function.

For convenience, the user can enter arguments directly into the
calculation functions.

     set/get inputs and outputs:

            A
            |    cache
      S <-> C <->  O
            |
            O

A calculation function sends the arguments to the set function which
does its work and returns the input cache. Then the calculation function
checks the output cache, does its work and returns the output cache.

Most calculation functions have a corresponding plot function. There is
a complicated relationship between them. In the R console, the user
probably wants calculations. In RShiny, the user get the plots. In the R
console, the calculation functions may be in charge but in RShiny, the
plot functions are in charge.

In the diagrams below, P is a plot function and T is a flag set to true
and F is the same flag set to false. By default, the flag is set to T. A
set/get function manages the flag, as discussed below. There are two
possibilities if the calculation functions are in charge.

     Calculations in mind:

      F --> C        T --> C --> P
            |              |     |
            O              O    plot

The calculation function gets the flag. If the flag is false, the output
cache is returned. If the flag is true, the output cache is returned and
the plot routine is called which returns a plot.

If the plot routines are in charge, they prevent recursive calls from
the calculation functions by passing a caller id. The caller id tells
the calculation function not to call that plot function.

     Plots in mind:

               |               |
      C <-id-> P      C <-id-> P <-id-> C
               |               |
              plot            plots

The plot function passes its id to a calculation function, gets the
output cache in return, and returns a plot. Most plot functions manage
more than one calculation function and return several different types of
plots.

In summary, the RShiny app allows the user to change inputs and quickly
plot the results. The R6 Objects verify inputs, update outputs and
create plots. A user can choose how to interact with the R6 Objects:

      1. call a calculate function with arguments and get the outputs,
      2. call a set/get function, call a calculate function and get the
         outputs,
      3. call a set/get function, call a plot function and get a plot,
      4. call a calculate function, let it call a plot function and get
         both outputs and a plot.

RShiny uses number three. An R6 object in the console will usually use
number four. For speed and stability, The R6 object only manages lists.
The actual calculations are in Rcpp/RcppParallel functions. These
functions can also be called directly, without the R6 object, or
imported into other apps.

## Solution Methods:

From left to right in the matrix, the R Shiny app, R6 Objects and
Rcpp/RccpParallel levels are modularized into Analytical, Finite
Difference, Maximum Likelihood and Monte Carlo solution methods:

    +------------+------------+------------+------------+
    |            | Finite     | Maximum    | Monte      |
    | Analytical | Difference | Likelihood | Carlo      |
    |            |            |            |            |
    |   RShiny   |   RShiny   |   RShiny   |   RShiny   |
    |   R6       |   R6       |   R6       |   R6       |
    |   Rcpp     |   Rcpp     |   Rcpp     |   Rcpp     |
    +------------+------------+------------+------------+

For calling R6 Objects or Rcpp/RcppParallel functions from the console,
this is enough. Outputs are stored in the global environment and persist
for use amongst modules.

The RShiny app also requires persistence. Usually, RShiny apps are
reactive and recalculate as the user changes inputs. But the
calculations here require many inputs and create large outputs. The
reactivity in RShiny is turned off. Instead, another R6 Object is
created to communicate between modules.

                      +---------------+
           +----------|   OUProcess   |----------+
           |          +---------------+          |
           |            |           |            |
    +------------+------------+------------+------------+
    | Analytical | Finite     | Maximum    | Monte      |
    |            | Difference | Likelihood | Carlo      |
    +------------+------------+------------+------------+

OUProcess instantiates the other modules and stores pointers to them.
Each module has a pointer to OUProcess. If Maximum Likelihood creates
new estimates, for example, it sends the estimates to OUProcess which
calls the set/get functions in the other modules. Every module will have
the new estimates.

RShiny instantiates OUProcess which coordinates the modules to handle
reactivity, validation of inputs and persistence. But OUProcess is also
handy from the console. The user doesn't have to keep track of variables
in the global environment, doesn't have to unlist returns from functions
before passing to other functions and can type fewer commands.

Often, the user may only use one R6 Object to solve a problem. Each
object can be instantiated separately without OUProcess. There is no
inheritance.

Details are in the help for the R6 Objects:

      ?OUProcess
      ?Analytical
      ?FiniteDifference
      ?MaximumLikelihood
      ?MonteCarlo

The calculations are discussed in:

      ?Analytical_Rcpp
      ?FiniteDifference_Rcpp
      ?MaximumLikelihood_Rcpp
      ?MonteCarlo_Rcpp
      ?OptionalPackages

Available data sets for Maximum Likelihood estimation have metadata,
accessed by typing:

      ?Agric_
      ?Climate_
      ?Ecosys_
      ?Finance_
      ?OUP_

and then selecting from the drop-down list.

These and more can viewed by:

      OUPHelpList()
      OUPHelpView(2)

Either a name or an index in the help list can be entered.

The easiest way to query and read the data without changing the working
directory is:

      OUPDataList()
      OUPDataRead(10)

Entering data by number in the list saves typing.

Scripts to demonstrate the methods are in files in the 'demo' directory.
The most convenient way to find and run demos is:

      OUPDemoList()
      OUPDemoRun(2)

The index 2 identifies the second demo in the list.

An RShiny app demonstrates all the methods and data and can be launched
from the console:

      OUPShiny()

Finally, calculation functions print to the console, which is usually
truncated. Plotting doesn't give the numbers behind the plots. As a
convenience, every calculation and plot function copies tab-delimited
text to the clipboard. The text can be pasted as a table in most other
applications. Copying is painless, but you can turn it off from any
module, for example, from the Analytical module:

      A <- Analytical$new()
      A$set_flags(plotit=FALSE,copyit=FALSE)

For good measure, the flag which allows calculation routines to
automatically call plot routines is also set to FALSE.

Plotting is by Plotly. Every plot can be downloaded as an SVG or PNG
file. Optionally, the file format and dimensions can be set.

      A <- Analytical$new()
      A$set_plot_info(fileformat="png",filewidth=640,fileheight=480)

Fonts and colors can also be set. Then clicking a button on the plot
downloads it. SVG files are editable for fonts and colors and are
manageable for small 2D plots. Large 2D and 3D plots are rendered with
WebGL and downloaded by a screen capture to become a PNG file. The PNG
file may we wrapped in an SVG file, but it is still just a PNG file.
