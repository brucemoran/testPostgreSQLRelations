
shiny::shinyUI(fluidPage(

    # Application title
    shiny::titlePanel("Test PostgreSQL Relations"),

    shiny::tabsetPanel(
        shiny::tabPanel("Patient Data",
            shiny::mainPanel(fluid = TRUE, DT::dataTableOutput("maintable1"))
        ),
        shiny::tabPanel("Result Data",
            shiny::actionButton("click_res", icon = icon("poll"), label = "Result Tables"),
            shiny::mainPanel(fluid = TRUE, DT::dataTableOutput("maintable2"))
        )
    )
))
