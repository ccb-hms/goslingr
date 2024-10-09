#' Install Gos Python package
#'
#' This function wraps installation functions from [reticulate][reticulate::reticulate]
#' to install the Python package **gos**.
#'
#' This package uses the [reticulate][reticulate::reticulate] package
#' to make an interface with the [Gos](https://github.com/gosling-lang/gos)
#' Python package. To promote consistency in usage of **reticulate** among
#' different R packages, it is
#' [recommended](https://rstudio.github.io/reticulate/articles/package.html#installing-python-dependencies)
#' to use a common Python environment, called `"r-reticulate"`.
#'
#' Depending on your setup, you can create this environment using
#' `reticulate::conda_create()` or [reticulate::virtualenv_create()],
#' as described in this
#' [reticulate article](https://rstudio.github.io/reticulate/articles/python_packages.html#conda-installation).
#'
#' @param envname `character`, name of python environment into which to install (default is `'r-gosling'`)
#' @param version `character`, version of Gos to install. For general use of this package,
#' @param new_env `boolean`, should the existing python environment recreated before installation.
#' @param ... other arguments sent to `reticulate::py_install()`
#'
#' @return invisible `NULL`, called for side-effects
#'
#' @export
#'
install_gosling <- function(envname = "r-gosling",
                            version = NULL,
                            new_env = identical(envname, "r-gosling"),
                            ...) {
  
  if(new_env && reticulate::virtualenv_exists(envname))
    reticulate::virtualenv_remove(envname)
  
  gosling_pkg <- "gosling"
  
  if (!is.null(version)) {
    gosling_pkg <- paste(gosling_pkg, version, sep = "==")
  }
  
  reticulate::py_install(gosling_pkg, envname = envname, ...)
  invisible(NULL)
}