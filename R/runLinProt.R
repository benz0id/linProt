#' Launch the shiny app for the linProt package
#'
#'
#' The shiny app trains a linear model on a user supplied alignment and
#' accompanying functional labels. The app allows for the specification of
#' hyper parameters and encoding techniques used. The trained model is then used
#' to provide insights into how different residues and thier properties
#' contribute to protein function. Code adapted from (Silva A., 2022).
#'
#' @return No return value but opens up a shiny page.
#'
#' @references
#' Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A,
#' Borges B (2022). _shiny: Web Application Framework for R_. R package version 1.7.3,
#' <https://CRAN.R-project.org/package=shiny>.
#'
#' Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
#' Unpublished. URL https://github.com/anjalisilva/TestingPackage."
#'
#'
#' @export
#' @importFrom shiny runApp
runLinProt <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "linProt")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}


# [END]
