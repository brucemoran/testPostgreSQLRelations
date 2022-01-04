# **testPostgreSQLRelations**

## **Overview**

Testing use of multiple 'tables' (here CSV files but in reality would be read from Postgres DB) for resulting and annotation of panel testing on patient data in a diagnostic lab setting.

## **Usage**

Install package: `devtools::install_github("brucemoran/testPostgreSQLRelations")`

Launch Shiny app: `testPostgreSQLRelations::testPostgreSQLRelations()`

In `Results Data` tab, click `Result Tables` to join patient, result and annotation data.

## **Future Work**

* Create a mock report that is generated on result, annotation data.
* ?Actually use a Postgres interface, with login etc? (see [shinySetupPostgreSQL](https://github.com/brucemoran/shinySetupPostgreSQL/) for more on how this is achieved)
