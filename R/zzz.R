#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib GregsOUPR6, .registration = TRUE
#' @importFrom Rcpp  evalCpp
## usethis namespace: end
NULL

.onAttach <- function(libname,pkgname) {
  if (RcppParallelInstalled()) {
    if (RcppsitmoInstalled()) {
      packageStartupMessage(
        "\nGreetings!\n",
        "\u2139\ufe0f OUPShiny() to launch the RShiny app,\n",
        "\u2139\ufe0f OUPHelpList() of help topics,\n",
        "\u2139\ufe0f OUPDataList() of data sets,\n",
        "\u2139\ufe0f OUPDemoList() of demos.\n\n",
        "\u2705 Compiled with Rcpp and RcppParallel with sitmo.\n"
      )
    } else {
      packageStartupMessage(
        "\nGreetings!\n",
        "\u2139\ufe0f OUPShiny() to launch the RShiny app,\n",
        "\u2139\ufe0f OUPHelpList() of help topics,\n",
        "\u2139\ufe0f OUPDataList() of data sets,\n",
        "\u2139\ufe0f OUPDemoList() of demos.\n\n",
        "\u2705 Compiled with Rcpp and RcppParallel without sitmo.\n"
      )
    }
  } else {
    packageStartupMessage(
      "\nGreetings!\n",
      "\u2139\ufe0f OUPShiny() to launch the RShiny app,\n",
      "\u2139\ufe0f OUPHelpList() of help topics,\n",
      "\u2139\ufe0f OUPDataList() of data sets,\n",
      "\u2139\ufe0f OUPDemoList() of demos.\n\n",
      "\u2714\ufe0f Compiled with Rcpp but without RcppParallel.\n"
    )
  }
}
