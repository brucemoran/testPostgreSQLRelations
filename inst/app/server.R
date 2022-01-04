# server.R

server <- function(input, output, session) {

    ##read in CSV data
    data <- shiny::reactiveValues()

    shiny::observe({

        data$patient <- readr::read_csv(system.file(package = "testPostgreSQLRelations", "data/test_patient_data.csv"))
        data$results <- readr::read_csv(system.file(package = "testPostgreSQLRelations", "data/test_result_data.csv"))
        data$anno <- readr::read_csv(system.file(package = "testPostgreSQLRelations", "data/test_anno_data.csv"))

    })

    output$maintable1 <- DT::renderDataTable ({
        data$patient
    })

    shiny::observeEvent(input$click_res, {

        ##join Lab ID in paient and result data
        data$datres <- dplyr::left_join(data$patient, data$results)
        data$datresanno <- dplyr::left_join(data$datres, data$anno, by = c("Gene" = "Gene", "HGVSc" = "HGVSc"))
        output$maintable2 <- DT::renderDataTable ({
            data$datresanno
        })
    })
}
