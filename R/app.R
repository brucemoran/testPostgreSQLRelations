# app.R

#' testPostgreSQLRelations function
#'

testPostgreSQLRelations <- function(){

  ui <- source(system.file(package = "testPostgreSQLRelations",
                           "app/ui.R"))
  server <- source(system.file(package = "testPostgreSQLRelations",
                                                    "app/server.R"))$value
  shiny::shinyApp(ui, server)
}
