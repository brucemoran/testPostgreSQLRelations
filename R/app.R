# app.R

#' Run testPostgreSQLRelations
#' @return running shiny app
#' @rdname testPostgreSQLRelations
#' @export

testPostgreSQLRelations <- function(){

  ui <- source(system.file(package = "testPostgreSQLRelations",
                           "app/ui.R"))$value
  server <- source(system.file(package = "testPostgreSQLRelations",
                               "app/server.R"))$value
  shiny::shinyApp(ui, server)
}
