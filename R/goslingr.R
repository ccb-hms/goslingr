#' Title TODO
#'
#' Description TODO
#' 
#' @param spec Gosling specification either as an R list or as a json string.
#'
#' @export
goslingr <- function(spec, width = NULL, height = NULL) {
  
  # is spec is json string, lets make it an R list
  if (is.character(spec)) {
    spec <- jsonlite::parse_json(spec)
  }
  
  # serve/update local data urls
  spec <- serve_and_update_local_data(spec)
  
  x <- list(spec = spec)
  
  # create widget
  htmlwidgets::createWidget(
    name = 'goslingr',
    x = x,
    width = width,
    height = height,
    package = 'goslingr'
  )
}


# recursive function to traverse and update local data urls in nested lists
serve_and_update_local_data <- function(x) {
  if (is.list(x)) {
    for (i in seq_along(x)) {
      if (!is.null(names(x)) && names(x)[i] == 'data') {
        x[[i]] <- local_data_to_gos(x[[i]])
      } else {
        x[[i]] <- serve_and_update_local_data(x[[i]])
      }
    }
  }
  return(x)
}

# use gos (e.g. gos$multivec(...)) for local data url so that it gets served
local_data_to_gos <- function(data_entry) {
  
  # if not local url, return data unchanged
  data_url <- data_entry$url
  is_local_url <- !is.null(data_url) && file.exists(data_url)
  if (!is_local_url) return(data_entry)
  
  # otherwise we recreate data with call to
  # e.g. gos$multivec(...) which generates localhost urls
  data_type <- data_entry$type
  
  # ensure we have gos loaded from reticulate
  load_gos()
  
  data_entry$type <- NULL
  new_data <- do.call(gos[[data_type]], data_entry)
  return(new_data)
}

# need if first call to gos is gos[[data_type]](...) but not e.g. gos$multivec(...)
load_gos <- function() {
  invisible(gos$data)
}

#' Shiny bindings for goslingr
#'
#' Output and render functions for using goslingr within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a goslingr
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name goslingr-shiny
#'
#' @export
goslingrOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'goslingr', width, height, package = 'goslingr')
}

#' @rdname goslingr-shiny
#' @export
renderGoslingr <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, goslingrOutput, env, quoted = TRUE)
}