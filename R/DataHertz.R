# roxygen ----
#' Default data for the Ornstein-Uhlenbeck Process
#'
#' Data to estimate parameters, rho, mu and sigma, where rho is the rate of convergence,
#'  mu is the location and sigma is the scale.
#'
#' \itemize{
#'   \item tau: time variable
#'   \item z: state variable
#' }
#'
#' The data must be in a .csv (Comma Separated Values) file.  The first column should
#'  be times and the second column should be states of nature.  There can be more columns
#'  for times and states if you wish.  There can be blank entries.  The data will be
#'  cleaned and sorted by time before it is used.
#'
#' @docType data
#' @keywords datasets
#' @name MyData
#' @format csv file with at least 3 rows and 2 columns
NULL

#' Rates of convergence for the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate different rates of convergence, rho.
#'
#' \itemize{
#'   \item year: time variable in annual increments for all sample paths
#'   \item z1-z5: sample paths in sets of three, each set with the same pseudo-random shocks but different rates of convergence
#' }
#'
#' The rate of convergence, rho, determines the probability distribution of the
#'  estimated parameters and the correlation between two sets of parameters in
#'  hypothesis tests.  Small rho tends toward Brownian Motion, which does not
#'  converge.  Large rho tends toward a stationary or ergodic process which has
#'  converged every time it is observed.  In between is an Ornstein-Uhlenbeck
#'  Process which converges but has not yet converged.
#'
#' Parameters for Browian Motion have an Erlang distribution.  Parameters for
#'  a stationary or ergodic process have a Chi^2 distribution.  These distributions
#'  are special cases of a Gamma distribution.  In general, parameters for the
#'  Ornstein-Uhlenbeck Process have a Gamma distribution.
#'
#' The shape parameter, alpha, identifies the distribution with 0.5 <= alpha <=1. Chi^2
#'  has alpha = 0.5.  Erlang has alpha = 1.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_Convergence
#' @format csv file with 177 rows and 16 columns
NULL

#' Observation intervals of different lengths and the Ornstein-Uhlenbeck Process
#'
#' Monte Carlo simulation to demonstrate that observations are never 'missing'.
#'
#' \itemize{
#'   \item Day: time variable in days
#'   \item Equal: sample paths with equal observation intervals
#'   \item Unequal: sample paths with unequal observation intervals
#'   \item Average: sample paths filled in with the average of states of nature
#'   \item Means: sample paths filled in with Means
#'   \item Year: time variable as decimal year
#' }
#'
#' An observation consists of the initial state, the initial time, the terminal state
#'  and the terminal time.  The terminal time minus the initial time can differ among
#'  observations because each observation has its own mean and variance.  Longer
#'  observation intervals have smaller means and larger variances but, otherwise,
#'  parameters estimated from unequal observation intervals are statistically
#'  equivalent to parameters estimated from equal intervals.
#'
#' The usual time-series methods are special cases of estimating the Ornstein-Uhlenbeck
#'  Process.  They estimate location, mu, and scale, sigma, but not the rate of
#'  convergence, rho.  They eliminate rho by assuming weak stationarity.  The mean of
#'  each observation is assumed to have converged and lost its connection to its
#'  initial state.  The variance of each observation does not depend upon an initial
#'  state and may still be converging.  Observations will have equal variances if
#'  observation intervals are equal.
#'
#' For the Ornstein-Uhlenbeck Process, variances converge twice as fast as means.
#'  If variances are still converging, so are means.  Observations over equal time
#'  intervals will have equal variances but different means.  In other words, weak
#'  stationarity does not exist.  Any assumption of stationarity is strong stationarity
#'  in which both the means and variances have converged.
#'
#' Stationarity is an hypothesis to test, not an assumption to impose.  Goodness
#'  of fit for the Ornstein-Uhlenbeck Process tests for stationarity, and also
#'  for the other extreme of Brownian Motion.
#'
#' If you have a time series with measurements at sporadic times, don't fill in for
#'  observations that aren't missing.  If you must, first estimate the parameters.
#'  Then use the parameters to calculate means for the times of your choosing.
#'  Means are the maximum likelihood estimates of unobserved states of nature.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_NotMissing
#' @format csv file with 366 rows and 6 columns
NULL

#' Observation intervals and the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate the effect of the observation interval on the
#'  rate of convergence, rho, and the scale, sigma.
#'
#' \itemize{
#'   \item Day: time variable in daily increments
#'   \item z: sample paths with daily observations
#'   \item Year: time variable as decimal year
#' }
#'
#' Observation intervals are measured in two different units:  days and years.
#'
#' Both sets of estimates have the same rho(t-s), where rho is the rate of convergence and
#'  t-s is the observation interval.  For observation intervals measured in days, t-s=1 and
#'  rho=1/2.  For observation intervals measured in years, t-s=1/365 and rho=365/2.  Location
#'   parameter, mu, is the same in both sets of estimates.  Scale parameter, sigma, is larger
#'   if rho is larger because the asymptotic variance, sigma^2/2rho, is the same in both sets
#'   of estimates.  If rho is 365 times bigger, sigma is 365^0.5 = 19.105 times bigger.
#'
#' Because rho(t-s) is the same, alpha = 0.5(1+exp(-rho(t-s))) is also the same, where alpha
#'  identifies the probability distribution of the estimates.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_ObservationInterval
#' @format csv file with 366 rows and 3 columns
NULL

#' Sample sizes for the Ornstein-Uhlenbeck Process
#'
#' Monte-Carlo simulation to demonstrate the effect of sample size on discovering the
#'  parameters of the Ornstein-Uhlenbeck Process.
#'
#' \itemize{
#'   \item year: time variable in annual increments for all sample paths
#'   \item small: 5 observations
#'   \item medium: 50 observations
#'   \item large: 500 observations
#' }
#'
#' The first 5 and the first 50 observations are the same as those of 500 simulated
#'  observations.  More observations get closer to the true parameters of rho=0.5,
#'  mu=15 and sigma=50.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_SampleSize
#' @format csv file with 501 rows and 4 columns
NULL

#' Smoothed sample paths for the Ornstein-Uhlenbeck Process
#'
#' Monte Carlo simulation to demonstrate the effect of smoothing on parameter estimates.
#'
#' \itemize{
#'   \item year: time variable in annual increments
#'   \item Data: simulated sample path
#'   \item G uno-G nueve: nine successively smoother paths
#' }
#'
#' Simulated data is smoothed ten times.  First, the parameters are estimated
#'  and the means calculated for each observation. Then the calculated means
#'  are used, as if they are data, in a subsequent estimation. Means are
#'  calculated again and used in the next estimation and so on ten times.
#'  The true rate of convergence, rho, and location, mu, are recovered from
#'  each estimation, but the scale, sigma, goes toward zero.
#'
#' By the eighth smoothing, the log of the likelihood becomes positive. Hence,
#'  the likelihood, as the anti-log, becomes greater than one. In other words,
#'  if the eighth and ninth smoothings were real samples, the probability of
#'  observing them would be greater than 100%. Further smoothings would increase
#'  this probability.  A small sigma is a tell-tale sign the data has been smoothed.
#'  A positive log likelihood is a sure sign.
#'
#' This is the best possible smoothing method. The model used for the smoothing
#'  is consistent with the data. Surprisingly, hypothesis tests and decision
#'  thresholds are not greatly affected.  However, passage times are calculated
#'  to be much larger and are completely unreliable.
#'
#' The best possible smoothing is unlikely. Any model used to smooth the data is
#'  probably not the Ornstein-Uhlenbeck Process.  The estimates will be wrong
#'  and the actual system will be much more uncertain. How much more uncertain
#'  is uncertain.
#'
#' Always use the raw data.
#'
#' @docType data
#' @keywords datasets
#' @name OUP_SmoothedData
#' @format csv file with 177 rows and 11 columns
NULL

#' Water in Farm Dams in the Riverina of New South Wales
#'
#' Stock water management to improve drought resilience in intensive-use grazing landscapes.
#'
#' \itemize{
#'   \item Run Days: time variable in daily increments
#'   \item Baseline: water volumes in cubic metres without a windbreak
#'   \item TwentyPct scenario: water volumes in cubic metres with a windbreak to suppress evaporation by 20%
#'   \item FortyPct scenario: water volumes in cubic metres with a windbreak to suppress evaporation by 40%
#' }
#'
#' MatLab simulated water volumes for a 4,000 cubic metre dam.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_FarmDamsRiverina
#' @format csv file with 19313 rows and 4 columns
#' @author Helena Clayton, Sally Thompson, Tim Capon, Greg Hertzler, Philip Graham, and David Lindenmayer, 2026
#' @source corresponding author: Tim Capon \email{tim.capon@csiro.au}
NULL

#' Farm Adaptation at Cootamundra, New South Wales
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S).
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMCootamundra
#' @format csv file with 50 rows and 5 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Narrendera, New South Wales
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S).
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMNarrendera
#' @format csv file with 50 rows and 5 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Temora, New South Wales
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item CWWCWW GM: Gross margins per hectare for CWWCWW rotation
#'   \item CWWSSS GM: Gross margins per hectare for CWWSSS rotation
#'   \item SSSSC GM: Gross margins per hectare for SSSSC rotation
#'   \item S GM: Gross margins per hectare for continuous S
#' }
#'
#' APSIM generated gross margins for wheat (W), canola (C) and sheep (S).
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_GMTemora
#' @format csv file with 50 rows and 5 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Soil Health with Stubble Management
#'
#' CSIRO Harden Long-Term Tillage Experiment. v3. CSIRO. Data Collection.
#'
#' \itemize{
#'   \item Year: time variable in sporadic increments
#'   \item N Burn: kilograms per hectare of mineral soil nitrogen in the soil profile from 0 to 160 millimetres with stubble burning
#'   \item N Bash: kilograms per hectare of mineral soil nitrogen in the soil profile from 0 to 160 millimetres with stubble bashing
#'   \item W Burn: millimetres of water in the soil profile from 0 to 160 millimetres with stubble burning
#'   \item W Bash: millimetres of water in the soil profile from 0 to 160 millimetres with stubble bashing
#' }
#'
#' Soil nitrogen and water in the top 160 millimetres.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_NSW_SoilHealthHarden
#' @format csv file with 30 rows and 5 columns
#' @author Kirkegaard, John; & Lilley, Julianne, 2023
#' @source https://doi.org/10.25919/2jqe-mz45
NULL

#' Farm Adaptation at Clare, South Australia
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMClare
#' @format csv file with 108 rows and 7 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Hawker, South Australia
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMHawker
#' @format csv file with 108 rows and 7 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Farm Adaptation at Orroroo, South Australia
#'
#' Will primary producers continue to adjust practices and technologies, change production systems or transform their industry – An application of real options.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Annual rain: Annual rainfall in millimetres
#'   \item Apr-Oct rain: Growing season rainfall in millimetres
#'   \item Wheat Yld: Wheat yield in tonnes per hectare
#'   \item Sheep DSE: Dry-sheep equivalents per hectare
#'   \item Wheat GM: Wheat gross margins per hectare
#'   \item Sheep GM: Sheep gross margins per hectare
#' }
#'
#' APSIM generated Wheat yields and manually simulated dry-sheep equivalents.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_GMOrroroo
#' @format csv file with 108 rows and 7 columns
#' @author Greg Hertzler, Todd Sanderson, Tim Capon, Peter Hayman, Ross Kingwell, Anthea McClintock, Jason Crean and Alan Randall, 2013
#' @source https://nccarf.edu.au/will-primary-producers-continue-adjust-practices-and-technologies-change-production/
NULL

#' Long Term Crop Rotations, South Australia
#'
#' Waite Permanent Rotation Trial. v4. CSIRO. Data Collection.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Plt 1-2 WF W: Wheat Fallow rotation, Wheat yields in kilograms per hectare
#'   \item Plt 3-4 WPe W: Wheat Peas rotation, Wheat yields in kilograms her hectare
#'   \item Plt 3-4 WPe Pe: Wheat Peas rotation, Pea yield in kilograms per hectare
#'   \item Plt 5-7 WPF W: Wheat Pasture Fallow rotation, Wheat yield in kilograms per hectare
#'   \item Plt 5-7 WPF P: Wheat Pasture Fallow rotation, Pasture dry matter in kilograms per hectare
#'   \item Plt 8-10 WOF W: Wheat Oats Fallow rotation, Wheat yield in kilograms per hectare
#'   \item Plt 8-10 WOF O: Wheat Oats Fallow rotation, Oat yield in kilograms per hectare
#'   \item Plt 11-16 WWPPPP W: Wheatx2 Pasturex4 rotation, Wheat yield in kilograms per hectare
#'   \item Plt 11-16 WWPPPP P: Wheatx2 Pasturex4 rotation, Pasture dry matter in kilograms per hectare
#'   \item Plt 17 W W: continuous Wheat, Wheat yield in kilograms per hectare
#'   \item Plt 18-20 WPP W: Wheat Pasturex2 rotation, Wheat yield in kilograms per hectare
#'   \item Plt 18-20 WPP P: Wheat Pasturex2 rotation, Pasture dry matter in kilograms per hectare
#'   \item Plt 21-23 WBPe W: Wheat Barley Peas rotation, Wheat yield in kilograms per hectare
#'   \item Plt 21-23 WBPe B: Wheat Barley Peas rotation, Barley yield in kilograms per hectare
#'   \item Plt 21-23 WBPe Pe: Wheat Barley Peas rotation, Pea yield in kilograms per hectare
#'   \item Plt 24-27 WOPF W: Wheat Oats Pasture Fallow rotation, Wheat yield in kilograms per hectare
#'   \item Plt 24-27 WOPF O: Wheat Oats Pasture Fallow rotation, Oat yield in kilograms per hectare
#'   \item Plt 24-27 WOPF P: Wheat Oats Pasture Fallow rotation, Pasture dry matter in kilograms per hectare
#'   \item Plt 28-29 P P: continuous Pasture, Pasture dry matter in kilograms per hectare
#'   \item Plt 30-33 WPPF W: Wheat Pasturex2 Fallow rotation, Wheat yield in kilograms per hectare
#'   \item Plt 30-33 WPPF W: Wheat Pasturex2 Fallow rotation, Pasture dry matter in kilograms per hectare
#'   \item Plt 34-35 WF W: Wheat Fallow rotation, Wheat yield in kilograms per hectare
#'   \item Annual Rain: Annual rainfall in millimetres
#'   \item Apr-Oct Rain: April through October rainfall in millimetres
#'   \item Plt 1 WF Carbon: Wheat Fallow rotation, Carbon in milligrams per hectare
#'   \item Plt 10 WOF Carbon: Wheat Oats Fallow rotation, Carbon in milligrams per hectare
#'   \item Plt 13 WWPPPP Carbon: Wheatx2 Pasturex4 rotation, Carbon in milligrams per hectare
#'   \item Plt 17 W Carbon: continuous Wheat, Carbon in milligrams per hectare
#'   \item Plt 29 P Carbon: continuous Pasture, Carbon in milligrams per hectare
#' }
#'
#' Experimental crop yields and pasture dry matter for permanent rotations from 1925 to 1993.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_SA_WaiteRotationTrial
#' @format csv file with 70 rows and 30 columns
#' @author Sanderman, Jonathan; David, Rakesh; Moore, Andrew; Keith, Heather; & Farquharson, Ryan, 2015
#' @source https://doi.org/10.4225/08/55E5165EC0D29
NULL

#' Tree shelter belts in Tasmania
#'
#' CSIRO Perennial Prosperity Project.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item With trees: sheep gross margins per paddock with tree shelter belt
#'   \item Without trees: sheep gross margins per paddock without tree shelter belt
#' }
#'
#' Simulated effects of shelter belts with prices and costs for sheep production.
#'
#' @docType data
#' @keywords datasets
#' @name Agric_Tas_TreeShelterBelts
#' @format csv file with 52 rows and 3 columns
#' @author Tim Capon, Daniel Mendham, 2024
#' @source corresponding author, Tim Capon \email{tim.capon@csiro.au}
NULL

#' Price of Carbon
#'
#' European Union Emission Trading System
#'
#' \itemize{
#'   \item Year: time variable as decimal year
#'   \item Open: Price at market open in EUR/tonne
#'   \item High: Daily high price in EUR/tonne
#'   \item Low: Daily low price in EUR/tonne
#'   \item Close: Price at market close in EUR/tonne
#'   \item Day: time variable in days since 1 January 2021
#' }
#'
#' Daily price of Carbon Permits from 1 January 2021 to 31 December 2025
#'
#' @docType data
#' @keywords datasets
#' @name Climate_CarbonCredits_EECXM
#' @format csv file with 1519 rows and 6 columns
#' @source https://seekingalpha.com/symbol/EECXM
NULL

#' Sea Level at Port Kembla
#'
#' Australian Baseline Sea Level Monitoring Project.
#'
#' \itemize{
#'   \item Day in 2023: time variable in hourly increments
#'   \item Sea Level: metres above Tide Gauge Zero
#'   \item Water Temperature: degrees Celsius
#'   \item Air Temperature: degrees Celsius
#'   \item Barometric Pressure: hPa
#' }
#'
#' Sea levels above Tide Guage Zero for 2023.
#'
#' This is radically smoothed data.  Do not use this data for decisions under uncertainty.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_SeaLevel_PortKembla
#' @format csv file with 7985 rows and 5 columns
#' @author Bureau of Meteorology
#' @source http://www.bom.gov.au/oceanography/projects/abslmp/data/data.shtml
NULL

#' Sunspot Numbers
#'
#' International Sunspot Number V2.0.
#'
#' \itemize{
#'   \item DecimalYear: time variable in daily increments
#'   \item Average: Simple average of the daily total sunspot number over all days of each calendar month.
#' }
#'
#' Sunspot numbers and groups of sunspot numbers, averaged by month.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_Sunspots
#' @format csv file with 3316 rows and 2 columns
#' @author Royal Observatory of Belgium
#' @source https://doi.org/10.24414/qnza-ac80
NULL

#' Maximum daily temperatures at Cape Otway
#'
#' Climate Data Online
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' Raw (unhomogenized) data.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_CapeOtway
#' @format csv file with 8888 rows and 2 columns
#' @author Bureau of Meteorology
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 090015
NULL

#' Maximum daily temperatures at Darwin
#'
#' Climate Data Online
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' Raw (Unhomogenized) data.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_Darwin
#' @format csv file with 9127 rows and 2 columns
#' @author Bureau of Meteorology
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 014040
NULL

#' Maximum daily temperatures at Tennant Creek
#'
#' Climate Data Online
#'
#' \itemize{
#'   \item Year: time variable in daily increments
#'   \item Max: maximum temperature in degrees centigrade
#' }
#'
#' Raw (Unhomogenized) data.
#'
#' @docType data
#' @keywords datasets
#' @name Climate_TempsMax_TennantCreek
#' @format csv file with 9115 rows and 2 columns
#' @author Bureau of Meteorology
#' @source http://www.bom.gov.au/climate/data/?ref=ftr, Station number 015135
NULL

#' Albatross Egg Counts
#'
#' The Australian Threatened Species Index, 2024.
#'
#' \itemize{
#'   \item Year: time variable in annual increments
#'   \item Wandering Albatross: egg count
#'   \item Black-browed Albatross: egg count
#'   \item Grey-headed Albatross: egg count
#' }
#'
#' Nesting egg counts of three Albatross species in Tasmania.
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_Albatross
#' @format csv file with 60 rows and 4 columns
#' @author Australian Government, Department of Climate Change, Energy, the Environment and Water
#' @source https://tsx.org.au/tsx2024
NULL

#' Water Supply for Irrigated Agriculture in Eastern Australia
#'
#' National data for Australia's ecosystem services: irrigated agricultural water supply (250 m resolution): 2000-01 to 2022-23. v7. CSIRO. Data Collection.
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item WaterSupply: irrigation water supplied in Australia in megalitres
#' }
#'
#' Regions are aggregated for Eastern Australia.
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_IrrigationWaterSupply
#' @format csv file with 13 rows and 2 columns
#' @author Liu, Ning; Smith, Greg; Evans, David; Tetreault Campbell, Sally; Pascoe, Sean; & Schmidt, Becky, 2024
#' @source https://doi.org/10.25919/7jj7-8826
NULL

#' Kangaroo Population and Harvest in South Australia
#'
#' Population surveys and meat processor records of harvest.
#'
#' \itemize{
#'   \item Year: time variable
#'   \item Red Pop: Estimated population of Red Kangaroos in 1000 head
#'   \item Red Harvest: Number of Red Kangaroos commercially harvested in 1000 head
#'   \item Grey Pop: Estimated population of Western Grey Kangaroos in 1000 head
#'   \item Grey Harvest: Number of Western Grey Kangaroos commercially harvested in 1000 head
#'   \item Euro Pop: Estimated population of European Kangaroos in 1000 head
#'   \item Euro Harvest: Number of European Kangaroos commercially harvested in 1000 head
#' }
#'
#' Regions are aggregated for all of South Australia.
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_Kangaroos
#' @format csv file with 46 rows and 7 columns
#' @author Government of South Australia, Department of Environment
#' @source https://www.environment.sa.gov.au/topics/animals-and-plants/sustainable-use-of-animals-and-plants/kangaroo-conservation-and-management/quotas-harvest-data
NULL

#' Southern Bluefin Tuna in Australian Waters
#'
#' National data for Australia's ecosystem services: fisheries biomass provisioning services (2000-01 to 2020-21). v9. CSIRO. Data Collection.
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item Catch: kilograms caught per year
#'   \item GVP: Gross value product in dollars as market price times catch
#'   \item EV: Exchange value in dollars as the value of ecosystem services as if a market existed
#' }
#'
#' Regions are aggregated for all of Australia.
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_SouthernBluefinTuna
#' @format csv file with 22 rows and 4 columns
#' @author Pascoe, Sean; Liu, Ning; Scheufele, Gabriela; Tetreault Campbell, Sally; & Schmidt, Becky, 2024
#' @source https://doi.org/10.25919/nzrp-6702
NULL

#' Sydney Drinking Water Catchment
#'
#' WaterNSW, WaterInsight.
#'
#' \itemize{
#'   \item Day: time variable for days from August 2015
#'   \item All: All Storage volume in gigalitres
#'   \item Blue Mtns: Blue Mountains Dams volume in gigalitres
#'   \item Nepean: Nepean Dam volume in gigalitres
#'   \item Avon: Avon Dam volume in gigalitres
#'   \item Wingecarribe: Wingecarribe Reservoir volume in gigalitres
#'   \item Cordeaux: Cordeaux Dam volume in gigalitres
#'   \item Cataract: Cataract Dam volume in gigalitres
#'   \item Warragamba: Warragamba Dam volume in gigalitres
#'   \item Woronora: Woronora Dam volume in gigalitres
#'   \item Prospect: Prospect Reservoir volume in gigalitres
#'   \item Tallowa: Tallowa Dam volume in gigalitres
#'   \item Fitzroy: Fitzroy Falls Dam volume in gigalitres
#'   \item Year: time variable as decimal year
#' }
#'
#' Drinking Water Storage from August 2015 to July 2025.
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_SydneyWater
#' @format csv file with 120 rows and 14 columns
#' @author WaterNSW
#' @source https://waterinsights.waternsw.com.au/12964-sydney-drinking-water-catchment/storage
NULL

#' Tropical Rock Lobster in Australian Waters
#'
#' National data for Australia's ecosystem services: fisheries biomass provisioning services (2000-01 to 2020-21). v9. CSIRO. Data Collection.
#'
#' \itemize{
#'   \item Fin Year: time variable for financial year 1 July to 30 June
#'   \item Catch: kilograms caught per year
#'   \item GVP: Gross value product in dollars as market price times catch
#'   \item EV: Exchange value in dollars as the value of ecosystem services as if a market existed
#' }
#'
#' Regions are aggregated for all of Australia
#'
#' @docType data
#' @keywords datasets
#' @name Ecosys_TropicalRockLobsters
#' @format csv file with 22 rows and 4 columns
#' @author Pascoe, Sean; Liu, Ning; Scheufele, Gabriela; Tetreault Campbell, Sally; & Schmidt, Becky, 2024
#' @source https://doi.org/10.25919/nzrp-6702
NULL

#' Metals and Energy Commodities
#'
#' London Bullion Market Association, World Gold Council, GoldHub, Perth Mint and International Monetary Fund prices.
#'
#' \itemize{
#'   \item Day: time variable in days since 1 January 2009
#'   \item Gold: London Bullion Market fix at the Perth (Aus) Mint in USD/troy ounce
#'   \item Silver: London Bullion Market fix in USD/troy ounce
#'   \item Copper: London Metals Exchange spot price, Copper, Grade A, USD/tonne
#'   \item Iron Ore: Cleared for Export Tianjin port, Iron Ore Fines 62% FE spot price in USD/tonne
#'   \item WTI: West Texas Intermediate in USD/bbl
#'   \item Brent: Brent Crude in USD/bbl
#'   \item Year: time variable as decimal year
#' }
#'
#' Closing prices on the first day of each month from 1 January 2009 to 1 May 2025.
#'
#' @docType data
#' @keywords datasets
#' @name Finance_Commodities
#' @format csv file with 197 rows and 8 columns
#' @source https://files.marketindex.com.au/files/workbooks/commodities-workbook.xlsx
NULL

#' Futures Prices
#'
#' Kansas City Hard Red Wheat Futures
#'
#' \itemize{
#'   \item Day: time variable in days since 1 July 2024
#'   \item Open: Price at market open in USD/ton
#'   \item High: Daily high price in USD/ton
#'   \item Low: Daily low price in USD/ton
#'   \item Close: Price at market close in USD/ton
#'   \item Adj Close: Same as Close for futures
#'   \item Volume: Daily number of trades
#'   \item Year: time variable as decimal year
#' }
#'
#' Daily futures for September 2025 maturity from 1 July 2024 to 30 June 2025
#'
#' @docType data
#' @keywords datasets
#' @name Finance_KansasCity_WheatFutures
#' @format csv file with 250 rows and 8 columns
#' @source https://finance.yahoo.com/quote (search for KE=F)
NULL

#' Exchange Rates
#'
#' US Dollars per Australian Dollar
#'
#' \itemize{
#'   \item Day: time variable in days since 1 July 2024
#'   \item Open: Rate at session open in USD/AUD
#'   \item High: Daily high rate in USD/AUD
#'   \item Low: Daily low rate in USD/AUD
#'   \item Close: Rate at session close in USD/AUD
#'   \item Adj Close: Same as Close for exchange rates
#'   \item Year: time variable as decimal year
#' }
#'
#' Daily exchange rates from 1 July 2024 to 30 June 2025, with AUD as the base currency and USD as the quote currency.
#'
#' @docType data
#' @keywords datasets
#' @name Finance_USDAUD
#' @format csv file with 258 rows and 7 columns
#' @source https://finance.yahoo.com/quote (search for AUDUSD=X)
NULL
