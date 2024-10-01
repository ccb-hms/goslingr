library(shiny)
library(goslingr)

ui <- fluidPage(
  titlePanel("reactR HTMLWidget Example"),
  goslingrOutput('widgetOutput')
)

server <- function(input, output, session) {
  output$widgetOutput <- renderGoslingr(
    goslingr("Hello world!")
  )
}

shinyApp(ui, server)