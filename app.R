library(shiny)
library(goslingr)

ui <- fluidPage(
  titlePanel("Lollipop, lollipop, oh lolli lollipop"),
  goslingrOutput('lollipop')
)

server <- function(input, output, session) {
  output$lollipop <- renderGoslingr({
    lollipop_spec_fpath <- system.file('extdata', 'lollipop_plots_spec.json', package = 'goslingr')
    lollipop_spec <- jsonlite::read_json(lollipop_spec_fpath)
    
    goslingr(lollipop_spec)
  })
}

shinyApp(ui, server)