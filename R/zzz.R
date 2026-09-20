#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib GregsOUPR6, .registration = TRUE
#' @importFrom Rcpp  evalCpp
## usethis namespace: end
NULL

#' @rdname OptionalPackages
#' @usage  RcppParallelThreads()
#' @return int <- RcppParallelThreads()
#' @export
RcppParallelThreads = function()
{
  if(RcppParallelInstalled()) { return(RcppParallel::defaultNumThreads()) }
  else { return(1) }
}

.onAttach <- function(libname,pkgname) {
    msg <- paste0("\nGreetings!\n",
        "\u2139\ufe0f OUPShiny() to launch the RShiny app,\n",
        "\u2139\ufe0f OUPHelpList() of help topics,\n",
        "\u2139\ufe0f OUPDataList() of data sets,\n",
        "\u2139\ufe0f OUPDemoList() of demos.\n\n",
        "\u2705 Compiled with Rcpp, ")
    if (RcppParallelInstalled()) {  parallel <- paste0("with RcppParallel (",RcppParallelThreads()," threads), ") }
    else { parallel <- "without RcppParallel, "}
    if (RcppdqrngInstalled()) { dqrng <- "with dqrng, " }
    else { dqrng <- "without dqrng, " }
    if (RcppsitmoInstalled()) { sitmo <- "and with sitmo.\n" }
    else { sitmo <- "and without sitmo.\n" }
    packageStartupMessage(paste0(msg,parallel,dqrng,sitmo))
}
