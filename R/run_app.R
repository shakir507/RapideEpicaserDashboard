run_app <- function(...) {

  app_dir <- system.file(
    "shiny",
    package = "RapideEpicaserDashboard"
  )

  if (!nzchar(app_dir)) {
    stop("Dashboard files were not found in the installed package.")
  }

  shiny::runApp(appDir = app_dir, ...)
}
