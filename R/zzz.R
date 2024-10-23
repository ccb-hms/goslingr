#' Gos object
#'
#' Uses the reticulate framework to access the Gos API.
#'
#' The Gos Python package is exposed through the `gos` object.
#' You can create and add to chart using its methods and classes,
#' as outlined in the
#' [Gos Python documentation](https://gosling-lang.github.io/gos).
#'
#' In this package, use the `$` operator wherever you see the `.` operator
#' used in Python.
#'
#' @seealso [Gos Python documentation](https://gosling-lang.github.io/gos)
#' @export gos
#'
gos <- NULL

# =============================================================================
# Note to maintainers:
#
# To change the supported Python version, set the option in .onLoad
# =============================================================================
.onLoad <- function(libname, pkgname) {
  reticulate::use_virtualenv("r-gosling")
  gos <<- reticulate::import("gosling", delay_load = yield_r_to_python)
}

# setup reticulate to continually yield R <> Python
# needed to enable local file paths for gosling data
yield_r_to_python <- function() {
  
  # check if packages for local higlass server were installed
  py_packages <- reticulate::py_list_packages()
  higlass_server <- all(c('clodius', 'servir') %in% py_packages$package)
  
  if (higlass_server) {
    reticulate::py_run_string("from time import sleep")
    py_yield_and_register_next_yield <- function() {
      reticulate::py_eval("sleep(0.001)")
      later::later(py_yield_and_register_next_yield, .001)
      invisible()
    }
    later::later(py_yield_and_register_next_yield)
  }
}
