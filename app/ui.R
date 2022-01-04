#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Test PostgreSQL Relations"),

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
