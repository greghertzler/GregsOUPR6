# Default data for the Ornstein-Uhlenbeck Process

Data to estimate parameters, rho, mu and sigma, where rho is the rate of
convergence, mu is the location and sigma is the scale.

## Format

csv file with at least 3 rows and 2 columns

## Details

- tau: time variable

- z: state variable

The data must be in a .csv (Comma Separated Values) file. The first
column should be times and the second column should be states of nature.
There can be more columns for times and states if you wish. There can be
blank entries. The data will be cleaned and sorted by time before it is
used.
