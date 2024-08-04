##
# Shiny App - Dataanalysis Spox Assay 
# Analyse Spox Data - 2 Plates per assay
# Version Beta
# Daniel Stern RKI
# 2023-07-18
##

# 01: Load libraries
library("shinydashboard") #"shiny",
library("shinycssloaders")
library("DT")
library("readxl")
library("rio")
library("tidyverse")
library("lubridate")
library("ggthemes")
library("caret")
library("glmnet")
library("docxtractr")
library("officer")
library("reshape") 
library("Hmisc")
library("chron")
library("gdata") 
library("plyr")
library("stringr")
library("ggplot2")
library("minpack.lm")
library("stats")
library("irr")
library("rootSolve")
library("msm")




# Delete Workspace
#rm(list = ls(all.names = TRUE))


##
# 1) Load models, data and metadata ####
# Load models and data
source("./functions/function_read_multiplex.R", encoding = "UTF-8") ## Function for data import
source("./functions/functions_quantify.R", encoding = "UTF-8") ## Function for data quantification
source("./functions/function_normalize.R", encoding = "UTF-8") ## Add QC data
source("./functions/function_classification.R", encoding = "UTF-8") ## Function for classification


# Source all  drLumi Files to test, if it works without loading the package
files.sources = list.files("./functions/drLumi/")
sapply(paste("./functions/drLumi/",files.sources, sep = ""), source)

ecdata <- import("./data/EC_input.xlsx") %>%  # Load ec-file for quantification
  pivot_longer(c(-sample), names_to = "analyte", values_to = "ec")
load("./data/threshold.Rdata") # Load cutoff values for IgG and IgM 
threshold <- 
  threshold %>% 
  filter(assaytype == "Multiplex")

load("./data/LDAmodels.Rdata") # Load LDA models 
load("./data/NTLassoReg.Rdata") # Loat NT-Lasso Regressin data
load("./data/GLMmodels.RData") # Load GLM Models with all data
load("./data/threshold_MaxSpec.RData") # Load Cutoff-Values for 100% specificity

ATI_N_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "ATI.N" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()

Combined_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "Combined" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()



# 01_01:    Define dashboard header ####
header <- dashboardHeader(
  title = "Spox Analyse",
  titleWidth = 235
)

# 01_02:    Define dashboard sidebar ####
sidebar <- dashboardSidebar(
  width = 235,
  sidebarMenu(
    p(img(src = "rkiLogo.png", height = 90, width = 235)),
    menuItem("Dateneingabe", tabName = "dataInput", icon = icon("hdd")),
    menuItem("Tabellen Eingabe", tabName = "outputTable", icon = icon("table")),
    menuItem("Tabelle Zusammenfassung", tabName = "summaryTable", icon = icon("table")),
    menuItem("Ergebnisse Quantifizierung", tabName = "summaryQuantification", icon = icon("table")),
    menuItem("Plot Standardkurven", tabName ="standardPlot", icon = icon("chart-line")),
    #    menuItem("Graphen Testreichweite", tabName ="reachPlot", icon = icon("chart-line")),
    #    menuItem("Graphen Lagebericht", tabName ="lageberichtPlot", icon = icon("chart-line")),
    #    menuItem("Graphen QC", tabName = "qcPlot", icon = icon("chart-line")),
    menuItem("Tabellen Klassifizierung", tabName = "classifiedTable", icon = icon("table")),
    menuItem("Download Ergebnisse", tabName = "dataOutput", icon = icon("download"))
    #    menuItemOutput("selectCutoff")
  )
)

# 01_03:    Define dashboard body ####
body <- dashboardBody(
  tabItems(
    # First tab: data input
    tabItem(tabName = "dataInput",
            h2("Eingabedatei laden"),
            fluidRow(
              # First tab box: load files
              tabBox(
                title = "Dateneingabe",
                id = "dataInput",
                # First panel: Load report sheet file
                tabPanel("Ergebnisse Platte 1 laden",
                         fileInput("dataFile1", "Auswahl Ergebnisdatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx"))
                ),
                tabPanel("Ergebnisse Platte 2 laden",
                         fileInput("dataFile2", "Auswahl Eingabedatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx"))
                ),
                tabPanel("Experimentbeschreibung laden",
                         fileInput("metaWord", "Auswahl Eingabedatei [docx]",
                                   accept = c("application/msword", ".docx"))
                )
              )
            )
    ),
    
    # Second tab: output of input tables (leave out or at least tidy input for indiviual labs)
    tabItem(tabName = "outputTable",
            h2("Eingabedateien Kontrolle"),
            fluidRow(
              tabBox(
                title = "Daten Eingabe",
                id = "outputTable",
                tabPanel("Daten Bioplex Platte 1",
                         DTOutput("inputData1_data_output")
                ),
                tabPanel("Daten Bioplex Platte 2",
                         DTOutput("inputData2_data_output")
                ),
                tabPanel("Metadaten",
                         DTOutput("inputmetadata_data_output")
                ),
                width = 12
              ) 
            )
    ),
    # Third tab: output of summary tables
    tabItem(tabName = "summaryTable",
            h2("Zusammenfassung der Daten"),
            fluidRow(
              tabBox(
                title = "Zusammenfassung",
                id = "summaryTable",
                tabPanel("Daten mit Metadaten",
                         DTOutput("inputData_data_output")
                ),
                width = 12
              )
            )
    ),
    tabItem(tabName = "summaryQuantification",
            h2("Ergebnisse der Quantifizierung"),
            fluidRow(
              tabBox(
                title = "Quantifizierung",
                id = "summaryQuantification",
                tabPanel("Quantifizierte Daten",
                         DTOutput("outputData_data_output")
                ),
                width = 12
              )
            )
    ),
    tabItem(tabName = "classifiedTable",
            h2("Ergebnisse der Klassifizierung"),
            fluidRow(
              tabBox(
                title = "Klassifizierung",
                id = "classifiedTable",
                tabPanel("Klassifizierte Daten",
                         DTOutput("outputData_classifiedTable")
                ),
                width = 12
              )
            )
    ),
    tabItem(
      tabName = "standardPlot",
      h2("Plot der Standardkurven"),
      fluidPage(
        tabBox(
          id = "plotStandards",
          title = "Plot der Standardkurven",
          tabPanel(
            "Platte 1",
            plotOutput("plotStandardPlate1", width = "100%", height = "600px")
          ),
          tabPanel(
            "Platte 2",
            plotOutput("plotStandardPlate2", width = "100%", height = "600px")
          )
        )
      )
    ),
    tabItem(tabName = "dataOutput",
            h2("Download der Ergebnisse"),
            fluidPage(
              tabBox(
                title = "Download der Ergebnisse",
                id = "download",
                tabPanel(
                  "Zusammenfassung der Ergebnisse",
                  p(downloadLink("downloadOutput", "Zusammenfassung Ergebnisse"))
                ),
                # downloadQuantified
                tabPanel(
                  "Quantifizierter Input",
                  p(downloadLink("downloadQuantified", "Zusammenfassung Input"))
                ),
                #downloadStandard
                tabPanel(
                  "Standardkurven",
                  p(downloadLink("downloadStandard", "Standardkurven"))
                )
              )
            )
    )
  )
)

#             tabPanel("Zusammenfassung der Daten Serology",
#                      DTOutput("summary_data_serology_output")
#             ), 
#             tabPanel("Zusammenfassung der Daten Antigentests",
#                      DTOutput("summary_data_antigen_output")
#             ), 
#             tabPanel("Zusammenfassung Testreichweite Zusammen",
#                      DTOutput("summaryTidyDataReachCombined_output")
#             ),
#             tabPanel("Zusammenfassung PCR-Ergebnisse S-Gen",
#                      DTOutput("fileVOXCOSGenPCR_SingleLabsReactive_output")
#             ),
#             tabPanel("Zusammenfassung Sequenzierung",
#                      DTOutput("summary_sequencing_output")
#             ),
#             tabPanel("Zusammenfassung Testkapazität Tage Summenstatistik",
#                      DTOutput("testReachWeekFunctionReactive_output")
#             ),
#             tabPanel("Zusammenfassung Rückstau",
#                      DTOutput("summaryTidyDataBacklogReactive_output")
#             ),
#             tabPanel("Summe aller PCR Tests",
#                      DTOutput("summariseInputCombinedReactive_output")
#             ),
#             width = 12
#           ) 
#         )
# ),
# # Fourth tab: Plot output from raw data curves
# tabItem(
#   tabName = "outPlot",
#   h2("Graphen Zusammenfassung"),
#   fluidPage(
#     tabBox(
#       id = "plotOutputRatio",
#       title = "Ratio Positive (%)",
#       tabPanel(
#         "PCR",
#         plotOutput("summaryQCRatioPlotOutput", width = "100%", height = "300px")
#       ),
#       tabPanel(
#         "Serology",
#         plotOutput("summaryQCRatioPlotSerology", width = "100%", height = "300px")
#       )
#     ),
#     tabBox(
#       id = "plotOutputLabs",
#       title = "Anzahl der Labore",
#       tabPanel(
#         "PCR",
#         plotOutput("summaryLabsPlotOutput", width = "100%", height = "300px")
#       ),
#       tabPanel(
#         "Serology",
#         plotOutput("summaryLabsPlotSerologyOutput", width = "100%", height = "300px")
#       )
#     ),
#     tabBox(
#       id = "plotOutputTests",
#       title = "Anzahl der Tests",
#       tabPanel(
#         "PCR",
#         plotOutput("summaryTestsPlotOutput", width = "100%", height = "300px")
#       ), 
#       tabPanel(
#         "Serology",
#         plotOutput("summaryTestsPlotSerologyOutput", width = "100%", height = "300px")
#       )
#     ),
#     tabBox(
#       id = "plotOutputPosives",
#       title = "Anzahl positiver Tests",
#       tabPanel(
#         "PCR",
#         plotOutput("summaryPositivePlotOutput", width = "100%", height = "300px")
#       ), 
#       tabPanel(
#         "Serology",
#         plotOutput("summaryPositivePlotSerologyOutput", width = "100%", height = "300px")
#       )
#     )
#   )
# ),
# 
# #Fifth tab: QC plots
# tabItem(tabName = "qcPlot",
#         h2("Graphen QC"),
#         fluidRow(
#           tabBox(
#             title = "QC ALM VOXCO",
#             id = "outputTable",
#             tabPanel("Anteil positver Ergebnisse (%)",
#                      plotOutput("individualQCRatioPlotOutput", height = 300)
#             ),
#             tabPanel("Interaktiv: Anteil positver Ergebnisse (%)",
#                      plotlyOutput("individualQCRatioPlotlyOutput",  height = "600px") %>% withSpinner(color="#f4b943")
#             ),
#             tabPanel("VOXCO ALM All",
#                      plotOutput("numberTestPlotAllOutput")
#             ),
#             tabPanel("VOXCO ALM 1",
#                      plotOutput("numberTestPlot50Output")
#             ),
#             tabPanel("VOXCO ALM 2",
#                      plotOutput("numberTestPlot75Output")
#             ),
#             tabPanel("VOXCO ALM 3",
#                      plotOutput("numberTestPlot104Output")
#             ),
#             tabPanel("VOXCO ALM 4",
#                      plotOutput("numberTestPlot133Output")
#             ),
#             tabPanel("VOXCO ALM 5",
#                      plotOutput("numberTestPlot157Output")
#             ),
#             tabPanel("VOXCO ALM 6",
#                      plotOutput("numberTestPlot191Output")
#             ),
#             tabPanel("VOXCO ALM 7",
#                      plotOutput("numberTestPlot230Output")
#             ),width = 12
#           )
#        )
#),

#Fifth tab: QC plots
# tabItem(tabName = "lageberichtPlot",
#         h2("Lagebericht Plot"),
#         fluidRow(
#           tabBox(
#             title = "Lagebericht Abbildung",
#             id = "outputLagebericht",
#             tabPanel("",
#                      plotOutput("LabsizeTablePlotOutput", height = 500)
#             ),width = 12
#           )
#         )
# ),
# 
# 
# 
# #Sixth tab: Reach plots
# tabItem(tabName = "reachPlot",
#         h2("Graphen Reichweite"),
#         fluidRow(
#           tabBox(
#             title = "Testreichweite",
#             id = "reachPlot",
#             tabPanel("Reichweite Tests pro Woche",
#                      plotOutput("reachPlotOutput", height = 300)
#             ),
#             width = 12
#           )
#         )
# ),
# 
# # Seventh tab: output of qc tables
# tabItem(tabName = "qcTable",
#         h2("Qualitätskontrolle der Daten"),
#         fluidRow(
#           tabBox(
#             title = "Qualitätskontrollen",
#             id = "summaryTable",
#             tabPanel("QC der Einzellabore",
#                      DTOutput("QCfailedSingleOutput"),
#                      p("Kontrolliert die folgenden Qualitätskriteren"),
#                      p("- QCpatients: Anzahl Patienten < Anzahl Tests (TRUE)"),
#                      p("- QCPositiveTests: Anzahl positiver Tests < Anzahl Tests (TRUE)"),
#                      p("- QCNegativeTests: Anzahl negativer Tests < Anzahl Tests (TRUE)"),
#                      p("- QCSumTest: Summe positiver und negativer Tests < Anzahl Tests (TRUE)")
#             ),
#             tabPanel("Abweichungen bei Positivrate",
#                      DTOutput("outliersRatioPositiveOutput"), 
#                      helpText("Zeigt Labore an, bei denen der Anteil der positiv
#                               getesteten über 50% liegt (unerwartet)")
#             ),
#             width = 12
#           ) 
#         )
# ),
# 
# # Last Tab: Download
# tabItem(tabName = "dataOutput",
#         h2("Download der Ergebnisse"),
#         fluidPage(
#           tabBox(
#             title = "Download der Ergebnisse",
#             id = "download",
#             tabPanel(
#               "Testzahlerfassung",
#               p(downloadLink("downloadTestzahlerfassung", "Testzahlerfassung final"))
#             ),
#             tabPanel(
#               "Testzahlerfassung Open Data",
#               p(downloadLink("downloadTestzahlerfassungOD", "Testzahlerfassung final Open Data")),
#               p(downloadLink("downloadTestzahlerfassungODxlsx", "Testzahlerfassung final Open Data Excel"))
#             ),
#             tabPanel(
#               "QC Daten",
#               p(downloadLink("downloadQCSummary", "QC Daten Zusammenfassung"))
#             ),
#             tabPanel(
#               "Zusammenfassung der Ergebnisse",
#               p(downloadLink("downloadSummary", "Zusammenfassung Ergebnisse")),
#               p(downloadLink("downloadSequenzierung", "Zusammenfassung Sequenzierung"))
#             ),
#             tabPanel(
#               "Ergebnisse Einzellabore",
#               p(downloadLink("downloadSingleLabs", "Ergebnisse Einzellabore aktuell")),
#               p(downloadLink("downloadSingleLabs1", "Ergebnisse Einzellabore 2020 KW 5-14")),
#               p(downloadLink("downloadSingleLabs2", "Ergebnisse Einzellabore 2020 KW 15-23")),
#               p(downloadLink("downloadSingleLabs3", "Ergebnisse Einzellabore 2020 KW 24-33")),
#               p(downloadLink("downloadSingleLabs4", "Ergebnisse Einzellabore 2020 KW 34-42")),
#               p(downloadLink("downloadSingleLabs5", "Ergebnisse Einzellabore 2020 KW 43-53"))
#             ),
#             tabPanel(
#               "Abbildungen",
#               p(downloadLink("downloadFigLagebericht", "Abbildungen Lagebericht")),
#               p(downloadLink("downloadFigBacklog", "Abbildungen Rückstau")),
#               p(downloadLink("downloadFigPCR", "Abbildungen PCR")),
#               p(downloadLink("downloadFigSerology", "Abbildungen Serologie")),
#               p(downloadLink("downloadFigBoxplotQC", "Abbildungen Boxplot QC"))
#             )
#           )
#         )
#  )

#)


# 01_04:    Combine header, sidebar, and body to ui ####
ui <- dashboardPage(header, sidebar, body)

# 02:       Server function ####
server <- function(input, output, session) {
  ##
  # File input
  # Import results from both files and word document with metadata
  # inputData <- reactive({
  #   req(input$dataFile1, input$dataFile2, input$metaWord)
  #   read_multiplex(input$dataFile1, input$dataFile2, input$metaWord)
  # })
  
  inputDataPlate1 <- reactive({
    req(input$dataFile1)
    read_plate_function(input$dataFile1$datapath, 1)
  })
  
  inputDataPlate2 <- reactive({
    req(input$dataFile2)
    read_plate_function(input$dataFile2$datapath, 2)
  })
  
  inputmetadata <- reactive({
    req(input$metaWord)
    read_metadata_function(input$metaWord$datapath, input$dataFile1$datapath, input$dataFile2$datapath)
  })
  
  inputData <- reactive({
    req(inputDataPlate1(), 
        inputDataPlate2(),
        inputmetadata())
    join_data_function(inputDataPlate1(), 
                       inputDataPlate2(),
                       inputmetadata() )
  })
  
  
  # Add QC flags
  inputDataQC <- reactive({
    req(inputData())
    addQC_function(inputData())
  })
  
  # Generate quantified results as list
  quantifiedList <- reactive({
    req(inputDataQC())
    generate_quantify_output_function(ecdata, inputDataQC(), threshold)
  })
  
  
  classifiedTable <- reactive({
    req(quantifiedList()[[3]])
    classify_sero_function(quantifiedList()[[3]],
                           lDAPanelsAll,
                           lDAPanelsDelta,
                           lDAPanelsHigh,
                           Combined_spec,
                           ATI_N_spec,
                           aucCombinedIgGPos)
  })
  
  
  # Output table to check file input
  # # 1.3) Output table for ALM dataInput
  output$inputData1_data_output <- renderDT({
    req(inputDataPlate1())
    inputDataPlate1()
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  output$inputData2_data_output <- renderDT({
    req(inputDataPlate2())
    inputDataPlate2()
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  output$inputmetadata_data_output <- renderDT({
    req(inputmetadata())
    inputmetadata()
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  
  output$inputData_data_output <- renderDT({
    req(inputmetadata())
    inputDataQC()
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  
  output$outputData_data_output <- renderDT({
    req(quantifiedList()[[3]])
    quantifiedList()[[3]]
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  
  
  output$outputData_classifiedTable <- renderDT({
    req(classifiedTable())
    classifiedTable()
  }, options = list(scrollX = TRUE, pageLength = 10), filter = "top")
  
  # Geneerate Plots for standards
  plotStandardPlate1In <- reactive({
    req(quantifiedList()[[1]])
    plot(quantifiedList()[[1]])
  })
  
  output$plotStandardPlate1 <- renderPlot(plotStandardPlate1In())
  
  plotStandardPlate2In <- reactive({
    req(quantifiedList()[[2]])
    plot(quantifiedList()[[2]])
  })
  
  output$plotStandardPlate2 <- renderPlot(plotStandardPlate2In())
  
  #  quantifiedList
  
  output$downloadQuantified <- downloadHandler(
    filename = function(){
      paste("Quantifizierter Input", ".xlsx", sep="")
    },
    content = function(file) {
      export(as.data.frame(quantifiedList()[[3]]),file)
    }
  )
  
  output$downloadStandard <- downloadHandler(
    filename = function(){
      paste("Standards", ".xlsx", sep="")
    },
    content = function(file) {
      export(as.data.frame(quantifiedList()[[4]]),file)
    }
  )
  
  
  output$downloadOutput <- downloadHandler(
    filename = function() {
      paste("Zusammenfassung Ergebnisse", ".xlsx", sep="")
    },
    content = function(file) {
      export(as.data.frame(classifiedTable()),file)
    }
  )
  
}

# 03:       Run shiny app ####

shinyApp(ui = ui, server = server)