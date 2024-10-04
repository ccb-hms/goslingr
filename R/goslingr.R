#' Title TODO
#'
#' Description TODO
#'
#' @import htmlwidgets
#'
#' @export
goslingr <- function(spec, width = NULL, height = NULL, elementId = NULL) {
  
  # describe a React component to send to the browser for rendering.
  content <- reactR::component(
    "GoslingComponent", 
    list(spec = spec)
  )
  
  params <- reactR::reactMarkup(content)
  
  # arguments used by jsonlite::toJSON
  # see https://www.htmlwidgets.org/develop_advanced.html#custom-json-serializer
  # 'rows' restores default from 'columns' in htmlwidgets
  attr(params, 'TOJSON_ARGS') <- list(dataframe = 'rows')
  
  # create widget
  htmlwidgets::createWidget(
    name = 'goslingr',
    x = params,
    width = width,
    height = height,
    package = 'goslingr',
    elementId = elementId
  )
}

#' Called by HTMLWidgets to produce the widget's root element.
#' @noRd
widget_html.goslingr <- function(id, style, class, ...) {
  htmltools::tagList(
    # Necessary for RStudio viewer version < 1.2
    reactR::html_dependency_corejs(),
    reactR::html_dependency_react(),
    reactR::html_dependency_reacttools(),
    htmltools::tags$div(id = id, class = class, style = style)
  )
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