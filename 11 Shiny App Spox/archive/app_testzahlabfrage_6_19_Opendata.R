##
# Shiny App - Testzahlabfrage 
# Check VOXCO Data Input fuer Janna
# Version 6.18
# Daniel Stern RKI
# 2022-09-20
# Verbesserungen zur Vorgängerversion
# 1) Entfernung doppelter Einträge in VOXCO und durch VOXCO hervorgerufen in ALM
# 2) Korrektur des Fehlers bei der Berechnung der Testreichweite aus ALM-Tabelle: Prospektive Woche wurde mit aktueller Woche aus VOXCO kombiniert
# 3) Entfernung von RKI_IDs ohne Eintrag aus ALM-Tabelle
# 4) Aufnahme und Kompatibilität mit den neuen Mustertabellen für die Abfrage der Antigentestungen
# 5) Dokumentation und Kommentierung
# 6) Aufnahme der Abfrage von S-Gen PCRs
# 7) Verbesserung der abgegriffenen Variablen aus VOXCO (Tests vs. Patienten)
# 8) Aufnahme der Abfrage der Sequenzierung
# 9) Erstellung eines neuen Plots für den Lagebericht
# 10) Erstellung neuer Tabellen für den Lagebericht -Automatisierung des Outputs
# 11) Umsetzen der Wünsche von Sindy für den automatisierten Bericht
# 12) Update der Testzahlabfrage für die VOC ab einschließlich KW09 zur getrennten Darstellung nach Sequenzierungs- und PCR-Ergebnissen
# 13) Updaten beim Zusammenfügen der VOC-Daten und Verbesserung der Ausgabe
# 14) Einbau der neuen Abbildung im Lagebericht ohne die Tabelle
# 15) Anpassung des Cutoff-Wertes für falsch positive Werte auf 75%
# 16) Einfügen der Opendata Auswertung
# 17) Erhöhung des Cutoff-Wertes auf 100% am 20.09.2022
##

# 01: Load libraries
packageList <- c("shiny",
                 "shinydashboard",
                 "shinycssloaders",
                 "DT",
                 "readxl",
                 "tidyverse",
                 "xlsx",
                 "ggsci",
                 "ggpubr",
                 "plotly",
                 "viridis",
                 "ggridges",
                 "jsonlite",
                 "lubridate",
                 "scales",
                 "gridExtra",
                 "grid",
                 "janitor",
                 "stringr")


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packageList)

rm(packageList, ipak)


# Delete Workspace
rm(list = ls(all.names = TRUE))


# Source functions
source("functions/functions_fileALM_data_input_year.R", encoding = "utf-8", local = knitr::knit_global())
source("functions/functions_fileVOXCO_data_input.R")
source("functions/functions_dataALMVOXCO_combine.R")
source("functions/functions_dataALMVOXCO_tidy.R")
source("functions/functions_summarise.R")
source("functions/functions_dataALMVOXCO_combine_serology.R")
source("functions/functions_dataALMVOXCO_tidy_serology.R")
source("functions/functions_dataALMVOXCO_combine_antigen.R")
source("functions/functions_dataALMVOXCO_tidy_antigen.R")
source("functions/functions_summarise_serology.R")
source("functions/functions_summarise_antigen.R")
source("functions/functions_plotting_Reach.R")
source("functions/functions_plotting_QC.R")
source("functions/functions_plotting_summary.R")
source("functions/functions_summarise_testreach.R")
source("functions/functions_summarise_QC.R")
source("functions/functions_summarise_S-Gen_PCR.R")
# NeU: Neue Funktionen ab KW9, um die Ergebnisse der VOCs getrennt voneinander einlesen zu können
# source("functions/functions_dataALMVOXCO_combine_sequencing.R")
#source("functions/functions_dataALMVOXCO_combine_sequencing_2.R")
#source("functions/functions_dataALMVOXCO_combine_sequencing_3_1_Indien.R", encoding = "utf-8", local = knitr::knit_global())
source("functions/functions_dataALMVOXCO_combine_sequencing_omicron.R", encoding = "utf-8", local = knitr::knit_global())
#source("functions/functions_plotting_Lagebericht.R")

# 01_01:    Define dashboard header ####
header <- dashboardHeader(
  title = "Testzahlabfrage",
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
    menuItem("Graphen Zusammenfassung", tabName ="outPlot", icon = icon("chart-line")),
    menuItem("Graphen Testreichweite", tabName ="reachPlot", icon = icon("chart-line")),
    menuItem("Graphen Lagebericht", tabName ="lageberichtPlot", icon = icon("chart-line")),
    menuItem("Graphen QC", tabName = "qcPlot", icon = icon("chart-line")),
    menuItem("Tabellen QC", tabName = "qcTable", icon = icon("table")),
    menuItem("Download Ergebnisse", tabName = "dataOutput", icon = icon("download")),
    menuItemOutput("selectCutoff")
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
                tabPanel("BUND ALM laden",
                         fileInput("ALM_data_input", "Auswahl Eingabedatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx")),
                         helpText("Vor Import alle Zellen als Zahl formatieren. 
                                  Reparatur PCR-Testung Ueberschrift in KW 10.")
                ),
                tabPanel("VOXCO Daten laden",
                         fileInput("VOXCO_data_input", "Auswahl Eingabedatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx"))
                ),
                tabPanel("ARS/RespVir Daten laden",
                         fileInput("ARSRespVir_data_input", "Auswahl Eingabedatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx"))
                ),
                tabPanel("ALM Sequenzierdaten laden",
                         fileInput("ALM_seq_data_input", "Auswahl Eingabedatei [xlsx]",
                                   accept = c("application/msexcel", ".xlsx"))
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
                tabPanel("ALM Daten Eingabe",
                         DTOutput("fileALM_data_output")
                ),
                tabPanel("VOXCO Daten Eingabe",
                         DTOutput("fileVOXCO_data_output")
                ),
                tabPanel("RespVir ARS Daten PCR",
                         DTOutput("fileARSRespVirPCR_data_output")
                ),
                tabPanel("RespVir ARS Daten Serologie",
                         DTOutput("fileARSRespVirSerology_data_output")
                ),
                tabPanel("ALM Sequenzierdaten",
                         DTOutput("fileALM_seq_data_output")
                ),
                tabPanel("VOXCO Sequenzierdaten",
                         DTOutput("fileVOXCO_seq_data_output")
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
                tabPanel("Zusammenfassung der Daten PCR",
                         DTOutput("summary_data_output")
                ),
                tabPanel("Zusammenfassung der Daten Serology",
                         DTOutput("summary_data_serology_output")
                ), 
                tabPanel("Zusammenfassung der Daten Antigentests",
                         DTOutput("summary_data_antigen_output")
                ), 
                tabPanel("Zusammenfassung Testreichweite Zusammen",
                         DTOutput("summaryTidyDataReachCombined_output")
                ),
                tabPanel("Zusammenfassung PCR-Ergebnisse S-Gen",
                         DTOutput("fileVOXCOSGenPCR_SingleLabsReactive_output")
                ),
                tabPanel("Zusammenfassung Sequenzierung",
                         DTOutput("summary_sequencing_output")
                ),
                tabPanel("Zusammenfassung Testkapazität Tage Summenstatistik",
                         DTOutput("testReachWeekFunctionReactive_output")
                ),
                tabPanel("Zusammenfassung Rückstau",
                         DTOutput("summaryTidyDataBacklogReactive_output")
                ),
                tabPanel("Summe aller PCR Tests",
                         DTOutput("summariseInputCombinedReactive_output")
                ),
                width = 12
              ) 
            )
    ),
    # Fourth tab: Plot output from raw data curves
    tabItem(
      tabName = "outPlot",
      h2("Graphen Zusammenfassung"),
      fluidPage(
        tabBox(
          id = "plotOutputRatio",
          title = "Ratio Positive (%)",
          tabPanel(
            "PCR",
            plotOutput("summaryQCRatioPlotOutput", width = "100%", height = "300px")
          ),
          tabPanel(
            "Serology",
            plotOutput("summaryQCRatioPlotSerology", width = "100%", height = "300px")
          )
        ),
        tabBox(
          id = "plotOutputLabs",
          title = "Anzahl der Labore",
          tabPanel(
            "PCR",
            plotOutput("summaryLabsPlotOutput", width = "100%", height = "300px")
          ),
          tabPanel(
            "Serology",
            plotOutput("summaryLabsPlotSerologyOutput", width = "100%", height = "300px")
          )
        ),
        tabBox(
          id = "plotOutputTests",
          title = "Anzahl der Tests",
          tabPanel(
            "PCR",
            plotOutput("summaryTestsPlotOutput", width = "100%", height = "300px")
          ), 
          tabPanel(
            "Serology",
            plotOutput("summaryTestsPlotSerologyOutput", width = "100%", height = "300px")
          )
        ),
        tabBox(
          id = "plotOutputPosives",
          title = "Anzahl positiver Tests",
          tabPanel(
            "PCR",
            plotOutput("summaryPositivePlotOutput", width = "100%", height = "300px")
          ), 
          tabPanel(
            "Serology",
            plotOutput("summaryPositivePlotSerologyOutput", width = "100%", height = "300px")
          )
        )
      )
    ),
    
    #Fifth tab: QC plots
    tabItem(tabName = "qcPlot",
            h2("Graphen QC"),
            fluidRow(
              tabBox(
                title = "QC ALM VOXCO",
                id = "outputTable",
                tabPanel("Anteil positver Ergebnisse (%)",
                         plotOutput("individualQCRatioPlotOutput", height = 300)
                ),
                tabPanel("Interaktiv: Anteil positver Ergebnisse (%)",
                         plotlyOutput("individualQCRatioPlotlyOutput",  height = "600px") %>% withSpinner(color="#f4b943")
                ),
                tabPanel("VOXCO ALM All",
                         plotOutput("numberTestPlotAllOutput")
                ),
                tabPanel("VOXCO ALM 1",
                         plotOutput("numberTestPlot50Output")
                ),
                tabPanel("VOXCO ALM 2",
                         plotOutput("numberTestPlot75Output")
                ),
                tabPanel("VOXCO ALM 3",
                         plotOutput("numberTestPlot104Output")
                ),
                tabPanel("VOXCO ALM 4",
                         plotOutput("numberTestPlot133Output")
                ),
                tabPanel("VOXCO ALM 5",
                         plotOutput("numberTestPlot157Output")
                ),
                tabPanel("VOXCO ALM 6",
                         plotOutput("numberTestPlot191Output")
                ),
                tabPanel("VOXCO ALM 7",
                         plotOutput("numberTestPlot230Output")
                ),width = 12
              )
            )
    ),
    
    #Fifth tab: QC plots
    tabItem(tabName = "lageberichtPlot",
            h2("Lagebericht Plot"),
            fluidRow(
              tabBox(
                title = "Lagebericht Abbildung",
                id = "outputLagebericht",
                tabPanel("",
                         plotOutput("LabsizeTablePlotOutput", height = 500)
                ),width = 12
              )
            )
    ),
    
    
    
    #Sixth tab: Reach plots
    tabItem(tabName = "reachPlot",
            h2("Graphen Reichweite"),
            fluidRow(
              tabBox(
                title = "Testreichweite",
                id = "reachPlot",
                tabPanel("Reichweite Tests pro Woche",
                         plotOutput("reachPlotOutput", height = 300)
                ),
                width = 12
              )
            )
    ),
    
    # Seventh tab: output of qc tables
    tabItem(tabName = "qcTable",
            h2("Qualitätskontrolle der Daten"),
            fluidRow(
              tabBox(
                title = "Qualitätskontrollen",
                id = "summaryTable",
                tabPanel("QC der Einzellabore",
                         DTOutput("QCfailedSingleOutput"),
                         p("Kontrolliert die folgenden Qualitätskriteren"),
                         p("- QCpatients: Anzahl Patienten < Anzahl Tests (TRUE)"),
                         p("- QCPositiveTests: Anzahl positiver Tests < Anzahl Tests (TRUE)"),
                         p("- QCNegativeTests: Anzahl negativer Tests < Anzahl Tests (TRUE)"),
                         p("- QCSumTest: Summe positiver und negativer Tests < Anzahl Tests (TRUE)")
                ),
                tabPanel("Abweichungen bei Positivrate",
                         DTOutput("outliersRatioPositiveOutput"), 
                         helpText("Zeigt Labore an, bei denen der Anteil der positiv
                                  getesteten über 50% liegt (unerwartet)")
                ),
                width = 12
              ) 
            )
    ),
    
    # Last Tab: Download
    tabItem(tabName = "dataOutput",
            h2("Download der Ergebnisse"),
            fluidPage(
              tabBox(
                title = "Download der Ergebnisse",
                id = "download",
                tabPanel(
                  "Testzahlerfassung",
                  p(downloadLink("downloadTestzahlerfassung", "Testzahlerfassung final"))
                ),
                tabPanel(
                  "Testzahlerfassung Open Data",
                  p(downloadLink("downloadTestzahlerfassungOD", "Testzahlerfassung final Open Data")),
                  p(downloadLink("downloadTestzahlerfassungODxlsx", "Testzahlerfassung final Open Data Excel"))
                ),
                tabPanel(
                  "QC Daten",
                  p(downloadLink("downloadQCSummary", "QC Daten Zusammenfassung"))
                ),
                tabPanel(
                  "Zusammenfassung der Ergebnisse",
                  p(downloadLink("downloadSummary", "Zusammenfassung Ergebnisse")),
                  p(downloadLink("downloadSequenzierung", "Zusammenfassung Sequenzierung"))
                ),
                tabPanel(
                  "Ergebnisse Einzellabore",
                  p(downloadLink("downloadSingleLabs", "Ergebnisse Einzellabore aktuell")),
                  p(downloadLink("downloadSingleLabs1", "Ergebnisse Einzellabore 2020 KW 5-14")),
                  p(downloadLink("downloadSingleLabs2", "Ergebnisse Einzellabore 2020 KW 15-23")),
                  p(downloadLink("downloadSingleLabs3", "Ergebnisse Einzellabore 2020 KW 24-33")),
                  p(downloadLink("downloadSingleLabs4", "Ergebnisse Einzellabore 2020 KW 34-42")),
                  p(downloadLink("downloadSingleLabs5", "Ergebnisse Einzellabore 2020 KW 43-53"))
                ),
                tabPanel(
                  "Abbildungen",
                  p(downloadLink("downloadFigLagebericht", "Abbildungen Lagebericht")),
                  p(downloadLink("downloadFigBacklog", "Abbildungen Rückstau")),
                  p(downloadLink("downloadFigPCR", "Abbildungen PCR")),
                  p(downloadLink("downloadFigSerology", "Abbildungen Serologie")),
                  p(downloadLink("downloadFigBoxplotQC", "Abbildungen Boxplot QC"))
                )
              )
            )
    )
  )
)


# 01_04:    Combine header, sidebar, and body to ui ####
ui <- dashboardPage(header, sidebar, body)

# 02:       Server function ####
server <- function(input, output, session) {
  
  
  
  
  
  #########################################################################################
  #########################################################################################
  # 02_01:    File input ####
  # 02_01_01: Generate responsive file inputs for report point table ####
  # Generate output tables to check for correct data input
  
  #########################################################################################
  # 1.1) Load responsive Input ALM data BUND
  
  fileALM_data_input_raw <- reactive({
    fileALM_data_input_raw_app_function(input$ALM_data_input)  
  })
  
  
  # 1.2) Tidy Datainput mit Funktion in functions_fileALM_data_input
  fileALM_data_input <- reactive({
    req(fileALM_data_input_raw())
    fileALM_data_input_function(fileALM_data_input_raw())
  })
  
  # 1.3) Output table for ALM dataInput
  output$fileALM_data_output <- renderDT({
    req(fileALM_data_input())
    fileALM_data_input()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  #########################################################################################
  
  
  
  #########################################################################################
  # 2.1) Load responsive Input VOXCO data
  fileVOXCO_data_input_raw <- reactive({
    inFileVOXCO_data_input <- input$VOXCO_data_input
    if (is.null(inFileVOXCO_data_input)){
      return(NULL)   
    }
    dataInput <- read_xlsx(inFileVOXCO_data_input$datapath, col_names = TRUE, col_types = "text") %>%
      rename(last_connection_date = "{Last Connection Date}", start_time_last_connection = "{Start time of last connection}") %>%
      select(!starts_with("{"))
    return(dataInput)
  })
  
  # 2.2) Tidy Datainput mit Funktion in functions_fileVOXCO_data_input
  # Backup: Alte Importfunktion
  #fileVOXCO_data_input <- reactive({
  #  req(fileVOXCO_data_input_raw())
  #  fileVOXCO_data_input_function(fileVOXCO_data_input_raw())
  #})
  
  # Neu: Import mit Verbesserung bei frühen Wochen und positive Tests statt Patienten am 29.01.2021
  fileVOXCO_data_input <- reactive({
    req(fileVOXCO_data_input_raw())
    fileVOXCO_data_input_function_tests(fileVOXCO_data_input_raw())
  })
  
  # 2.3) Output table for VOXCO dataInput
  output$fileVOXCO_data_output <- renderDT({
    req(fileVOXCO_data_input())
    fileVOXCO_data_input()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ##########################################################################################
  
  
  
  ##########################################################################################
  # 3.1 a)Responsive Input ARSRespVir PCR
  fileARSRespVirPCR_data_input <- reactive({
    fileARSRespVirPCR_data_input <- input$ARSRespVir_data_input
    if (is.null(fileARSRespVirPCR_data_input)){
      return(NULL)}
    dataInput <- read_excel(fileARSRespVirPCR_data_input$datapath, sheet = "ARSRespVir_PCR", col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
    dataInput$week <- as.integer(dataInput$week)
    dataInput <- dataInput %>%
      # New 05.01.2020 Fix year
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week))
    return(dataInput)
  })
  
  # 3.3.a) Output table for ARSRespVir PCR
  output$fileARSRespVirPCR_data_output <- renderDT({
    req(fileARSRespVirPCR_data_input())
    fileARSRespVirPCR_data_input()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  
  # 3.1. b) Responsive Input ARSRespVir Serology
  fileARSRespVirSerology_data_input <- reactive({
    fileARSRespVirSerology_data_input <- input$ARSRespVir_data_input
    if (is.null(fileARSRespVirSerology_data_input)){
      return(NULL)}
    dataInput <- read_excel(fileARSRespVirSerology_data_input$datapath, sheet = "ARSRespVir_Serology", col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
    dataInput$week <- as.integer(dataInput$week)
    dataInput <- dataInput %>%
      # New 05.01.2020 Fix year
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week))
    return(dataInput)
  })
  
  # 3.3. b) Responsive Input ARSRespVir Serology
  output$fileARSRespVirSerology_data_output <- renderDT({
    req(fileARSRespVirSerology_data_input())
    fileARSRespVirSerology_data_input()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  #########################################################################################
  
  
  
  #########################################################################################
  # Cutoff to determine positive ratio. Default 75% positive Rate treated as mistake
  # Can be changed via slider input
  output$selectCutoff <- renderMenu({
    req(fileALM_data_input(), fileVOXCO_data_input(), fileARSRespVirPCR_data_input())
    menuItem("Cutoff falsch positive Verhältnis", icon = icon("sliders-h"), sliderInput("cutoff", "", 0,100,95)
    )
  })
  ##########################################################################################
  
  
  
  
  
  ###########################################################################################
  ###########################################################################################
  # 02_02: Data manipulation and calculation ####
  # Based on functions definend in extral R scripts
  
  
  
  ###########################################################################################
  # 02_02_01: Define labels ####
  cutoff <- reactive(input$cutoff)
  # Define variables, that are included in either ARS or RespVir dataset
  vars <- c(institute = "RKI_ID",
            federal_state = "federal_state",
            tests = "tests_combined",
            patients = "patients_combined",
            negative = "negative_combined",
            positive = "positive_combined",
            week = "week",
            year = "year",
            year_week = "year_week",
            source = "source_combined")
  varsSerology <- c(institute = "RKI_ID",
                    federal_state = "federal_state",
                    tests = "tests_combined",
                    positive = "positive_combined",
                    week = "week",
                    year = "year",
                    year_week = "year_week",
                    source = "source_combined")
  # Neu 11.12.2020 Define Vars Antigentestung
  varsAntigen <- c(institute = "RKI_ID",
                   federal_state = "federal_state",
                   tests = "tests_combined",
                   positive = "positive_combined",
                   PCRconfirmed = "fracPCR_combined",
                   week = "week",
                   year = "year",
                   year_week = "year_week",
                   source = "source_combined")
  # Vector with RKI_IDs, for which duplicated values are included in either ARS or RespVir dataset
  duplicateRKI_IDs <- c(25, 30, 34, 160, 161, 162, 165, 168)
  ##########################################################################################
  
  
  
  ##########################################################################################
  # 02_02_02_A: Data wrangling for PCR data to generate File ready for QC Input
  # Select PCR Date and combine ALM data with VOXCO data
  fileALMVOXCOPCR_data_inputRaw <- reactive({
    req(fileALM_data_input(), fileVOXCO_data_input(), vars, duplicateRKI_IDs, cutoff())
    fileALMVOXCOPCR_data_inputRaw_function(fileALM_data_input(), fileVOXCO_data_input(), vars, duplicateRKI_IDs, cutoff())
  })
  
  # Tidy data, remove double data sets and missing data sets for positive results
  fileALMVOXCOPCR_data_input <- reactive({
    req(fileALMVOXCOPCR_data_inputRaw())
    fileALMVOXCOPCR_data_input_tidy_function(fileALMVOXCOPCR_data_inputRaw())
  })
  
  # Summarise data from ALM, VOXCO and ARS/Respvir
  summariseInputDataReactive <- reactive({
    req(fileALMVOXCOPCR_data_input(), fileARSRespVirPCR_data_input())
    summariseInputData(fileALMVOXCOPCR_data_input(), fileARSRespVirPCR_data_input())
  })
  
  # Output summary data  as table
  output$summary_data_output <- renderDT({
    req(summariseInputDataReactive())
    summariseInputDataReactive() %>%
      filter(labsClass == "Combined") %>%
      select("Jahr" = year,
             "KW* 2020" = week, 
             "Anzahl Testungen" = sumTests,
             "Positiv getestet" = sumPositive,
             "Positivquote (%)" = ratioPositive,
             "Anzahl übermittelnde Labore" = sumLabs)
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  
  # Calculate total summary for PCR results
  summariseInputCombinedReactive <- reactive({
    req(summariseInputDataReactive())
    sumTestTotal20 <- summariseInputDataReactive() %>%
      filter(year == 2020) %>%
      filter(week > 10) %>%
      filter(labsClass == "Combined") %>%
      summarise(totalTest = (sum(sumTests) + 124716),
                totalPositive = sum(sumPositive) + 3892)
    
    sumTestTotal21 <- summariseInputDataReactive() %>%
      filter(year == 2021) %>%
      filter(labsClass == "Combined") %>%
      summarise(totalTest = (sum(sumTests)),
                totalPositive = sum(sumPositive))
    
    sumTestTotal <- sumTestTotal20 %>%
      add_row(sumTestTotal21) %>%
      summarise(totalTest = sum(totalTest),
                totalPositive = sum(totalPositive))
  })
  
  # Output total output for PCR results
  output$summariseInputCombinedReactive_output <- renderDT({
    req(summariseInputCombinedReactive())
    summariseInputCombinedReactive()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ###############################################################################################
  
  
  
  ##############################################################################################
  # Select Serology data and combine ALM data with VOXCO data
  fileALMVOXCOSerology_data_inputRaw <- reactive({
    req(fileALM_data_input(), fileVOXCO_data_input(), varsSerology)
    fileALMVOXCOSerology_data_input_function(fileALM_data_input(), fileVOXCO_data_input(), varsSerology)
  })
  
  # Tidy data, remove double data sets and missing data sets for positive results
  fileALMVOXCOSerology_data_input <- reactive({
    req(fileALMVOXCOSerology_data_inputRaw())
    fileALMVOXCOSerology_data_input_tidy_function(fileALMVOXCOSerology_data_inputRaw())
  })
  
  # Summarise serology data
  summariseInputDataSerologyReactive <- reactive({
    req(fileALMVOXCOSerology_data_input(), fileARSRespVirSerology_data_input())
    summariseInputDataSerology(fileALMVOXCOSerology_data_input(), fileARSRespVirSerology_data_input())
  })
  
  # Output summary for serology data
  output$summary_data_serology_output <- renderDT({
    req(summariseInputDataSerologyReactive())
    summariseInputDataSerologyReactive() %>%
      select("Jahr" = year,
             "KW* 2020" = week, 
             "Anzahl Testungen" = sumTests,
             "Positiv getestet" = sumPositive,
             "Positivquote (%)" = ratioPositive,
             "Anzahl übermittelnde Labore" = sumLabs)
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ##################################################################################################
  
  
  
  ##############################################################################################
  # New: 14.12.2020 Select Antigen data and combine ALM data with VOXCO data
  fileALMVOXCOAntigen_data_inputRaw <- reactive({
    req(fileALM_data_input(), fileVOXCO_data_input(), varsAntigen)
    fileALMVOXCOAntigen_data_input_function(fileALM_data_input(), fileVOXCO_data_input(), varsAntigen)
  })
  
  
  # Tidy data, remove double data sets and missing data sets for positive results
  fileALMVOXCOAntigen_data_input <- reactive({
    req(fileALMVOXCOAntigen_data_inputRaw())
    fileALMVOXCOAntigen_data_input_tidy_function(fileALMVOXCOAntigen_data_inputRaw())
  })
  
  # Summarise antigen data
  summariseInputDataAntigenReactive <- reactive({
    req(fileALMVOXCOAntigen_data_input())
    summariseInputDataAntigen(fileALMVOXCOAntigen_data_input())
  })
  
  # Output summary for antigen data
  output$summary_data_antigen_output <-   renderDT({
    req(summariseInputDataSerologyReactive())
    summariseInputDataAntigenReactive() %>%
      select("Jahr" = year,
             "KW* 2020" = week, 
             "Anzahl Testungen" = sumTests,
             "Positiv getestet" = sumPositive,
             "Davon per PCR bestätigt" = sumPCRconfirmed,
             "Positivquote (%)" = ratioPositive,
             "Anzahl übermittelnde Labore" = sumLabs)
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ##################################################################################################
  # Auswertung der sequenzierungsdaten
  fileALM_data_input_seq_raw <- reactive({
    fileALM_data_input_seq_raw_app_function(input$ALM_seq_data_input)
  })
  
  
  # Tidy Datainput mit Funktion in functions_fileALM_data_input
  fileALMSeqData <- reactive({
    req(fileALM_data_input_seq_raw())
    fileALMSeqdata_function(fileALM_data_input_seq_raw())
  })
  
  # Output table for ALM dataInput
  output$fileALM_seq_data_output <- renderDT({
    req(fileALMSeqData())
    fileALMSeqData()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  
  # Datainput der Sequenzdaten aus den VOXCO Daten 
  fileVOXCOSeqData <- reactive({
    req(fileVOXCO_data_input())
    sequencing_VOXCO_input_function(fileVOXCO_data_input())
  })
  
  # Output table für VOXCO Daten Sequenzierung
  output$fileVOXCO_seq_data_output <- renderDT({
    req(fileVOXCOSeqData())
    fileVOXCOSeqData()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  
  
  # Zusammenführen der Sequenzierdaten aus ALM und VOXCO
  sequencinq_ALM_VOXCO_final <- reactive({
    req(fileALMSeqData(), fileVOXCOSeqData(), fileALMVOXCOPCR_data_input())
    combine_sequencing_function(fileALMSeqData(), fileVOXCOSeqData(), fileALMVOXCOPCR_data_input())
  })
  
  # Zusammenfassen nach Woche
  summary_sequencing <- reactive({
    req(sequencinq_ALM_VOXCO_final())
    summarise_sequencing_function(sequencinq_ALM_VOXCO_final())
  })
  
  
  # Zusammenfassen nach Woche
  summary_sequencing_ALM <- reactive({
    req(sequencinq_ALM_VOXCO_final())
    summarise_sequencing_function_ALM(sequencinq_ALM_VOXCO_final())
  })
  
  
  # Zusammenfassen nach Woche
  summary_sequencing_VOXCO <- reactive({
    req(sequencinq_ALM_VOXCO_final())
    summarise_sequencing_function_VOXCO(sequencinq_ALM_VOXCO_final())
  })
  
  # Output table für Zusammenfassung der Ergebnisse
  output$summary_sequencing_output <- renderDT({
    req(summary_sequencing())
    summary_sequencing()
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ##################################################################################################  
  
  
  
  
  
  
  
  
  ##################################################################################################
  # Generate Data summary for lab capacities
  fileALMVOXCOPCR_allLabs <- reactive({
    req(fileALM_data_input(), fileVOXCO_data_input())
    fileALMVOXCOPCR_allLabs_function(fileALM_data_input(), fileVOXCO_data_input())
  })
  ################################################################################################
  
  
  
  #################################################################################################
  # Combined table for test reach
  summaryTidyDataReachCombined <- reactive({
    req(fileALMVOXCOPCR_allLabs())
    summaryTidyDataReachCombined_function(fileALMVOXCOPCR_allLabs())
  })
  
  # Output test reach based combined
  output$summaryTidyDataReachCombined_output <- renderDT({
    req(summaryTidyDataReachCombined())
    summaryTidyDataReachCombined() %>%
      select("Jahr" = year,
             "KW*, für die die Angabe prognostisch erfolgt ist" = week,
             "Anzahl übermittelnder Labore (Arbeitstage)" = sumLabsWorkweek,
             "Anzahl übermittelnder Labore (Reichweite)" = sumLabsCapacity,
             "Testkapazität pro Tag" = sumDay,
             "Neu ab KW15: wöchentliche Kapazität anhand von Wochenarbeitstagen" = sumWeekWorkdays,
             "Reale Kapazität anhand von Testreichweite" = sumWeekCapacity,
             "Reale Kapazität anhand von Testreichweite und Wochenarbeitstagen (ALM Kriterien)" = sumWeekCapacityALM)
    #    summarise(sumLabsWorkweek = length(which(!is.na(testsCapacity_combined))),
    #          sumLabsCapacity = length(which(!is.na(testsCapacity_combined))),
    #          sumDay = sum(testsCapacity_combined, na.rm = TRUE),
    #          sumWeekWorkdays = sum(testCapacityWeek, na.rm = TRUE),
    #          sumWeekCapacity = sum(testCapacityWeekReach, na.rm = TRUE),
    #          sumWeekCapacityALM = sum(testCapacityWeekReachALM, na.rm = TRUE))
  }, options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  #################################################################################################
  
  
  
  
  
  
  
  #################################################################################################
  # Calculate backlog of tests
  summaryTidyDataBacklogReactive <- reactive({
    req(fileALMVOXCOPCR_allLabs())
    summaryTidyDataBacklogReactive_function(fileALMVOXCOPCR_allLabs())
  })
  
  # Output of test backlog
  output$summaryTidyDataBacklogReactive_output <- renderDT({
    req(summaryTidyDataBacklogReactive())
    summaryTidyDataBacklogReactive() %>%
      select("Jahr" = year,
             "KW* 2020" = week,
             "Anzahl übermittelnde Labore" = sumLabszero,
             "Rückstau" = sumBacklog)
  }, options = list(scrollX = TRUE, pageLength = 50, pageLength = 50), filter = "top")
  ##################################################################################################
  
  
  
  ##################################################################################################
  # Dev: 28.12.2020: Summary VOXCO S-Gen PCR-Ergebnisse und Zusammenfassung
  fileVOXCOSGenPCR_SingleLabsReactiveInput <- reactive({
    req(fileALMVOXCOPCR_allLabs())
    fileVOXCOSGenPCR_SingleLabs_function(fileALMVOXCOPCR_allLabs())
  })
  
  fileVOXCOSGenPCR_SingleLabsReactive <- reactive({
    req(fileVOXCOSGenPCR_SingleLabsReactiveInput())
    fileVOXCOSGenPCR_SingleLabsReactive_output_function(fileVOXCOSGenPCR_SingleLabsReactiveInput())
  })
  
  
  output$fileVOXCOSGenPCR_SingleLabsReactive_output <- renderDT({
    req(fileVOXCOSGenPCR_SingleLabsReactive())
    fileVOXCOSGenPCR_SingleLabsReactive()
  }, options = list(scrollX = TRUE, pageLength = 50, pageLength = 50), filter = "top")
  
  
  ##################################################################################################
  
  
  
  
  
  ##################################################################################################
  # Output summary statics for test reach per lab
  testReachWeekFunctionReactive <- reactive({
    req(fileALMVOXCOPCR_allLabs())
    testReachWeekFunction(fileALMVOXCOPCR_allLabs())
  })
  
  output$testReachWeekFunctionReactive_output <- renderDT({
    req(testReachWeekFunctionReactive())
    testReachWeekFunctionReactive()
  }, options = list(scrollX = TRUE, pageLength = 50, pageLength = 50), filter = "top")
  
  
  
  ##################################################################################################
  ##################################################################################################
  # Plot Test Reach Plots (Neu 01.10.2020)
  reachPlotFunctionReactive <- reactive({
    req(fileALMVOXCOPCR_allLabs())
    reachPlotFunction(fileALMVOXCOPCR_allLabs())
  })
  
  output$reachPlotOutput <- renderPlot(reachPlotFunctionReactive())
  
  
  
  # Plots summary comparison positive ratio
  summaryQCRatioPlotReactive <- reactive({
    req(summariseInputDataReactive())
    summaryQCRatioPlot(summariseInputDataReactive())
  })
  
  output$summaryQCRatioPlotOutput <- renderPlot(summaryQCRatioPlotReactive())
  
  # Plot summary comparison positive ratio serology
  summaryQCRatioPlotSerologyReactive <- reactive({
    req(summariseInputDataSerologyReactive())
    summaryQCRatioPlotSerology(summariseInputDataSerologyReactive())
  })
  
  output$summaryQCRatioPlotSerology <- renderPlot(summaryQCRatioPlotSerologyReactive())
  
  
  
  
  # Plot QC plots for labs summary for week > 4
  summaryLabsPlotReactive <- reactive({
    req(summariseInputDataReactive())
    summaryLabsPlot(summariseInputDataReactive())
  })
  
  output$summaryLabsPlotOutput <- renderPlot(summaryLabsPlotReactive())
  
  # Plot QC plots for labs summary serology for week > 4
  summaryLabsPlotSerologReactive <- reactive({
    req(summariseInputDataSerologyReactive())
    summaryLabsPlotSerology(summariseInputDataSerologyReactive())
  })
  
  output$summaryLabsPlotSerologyOutput <- renderPlot(summaryLabsPlotSerologReactive())
  
  
  
  # Plot summary of tests performed
  summaryTestsPlotReactive <- reactive({
    req(summariseInputDataReactive())
    summaryTestsPlot(summariseInputDataReactive())
  })
  
  output$summaryTestsPlotOutput <- renderPlot(summaryTestsPlotReactive())
  
  # Plot summary of tests performed serology
  summaryTestsPlotSerologyReactive <- reactive({
    req(summariseInputDataSerologyReactive())
    summaryTestsSerologyPlot(summariseInputDataSerologyReactive())
  })
  
  output$summaryTestsPlotSerologyOutput <- renderPlot(summaryTestsPlotSerologyReactive())
  
  
  
  # Plot summary of positive plots
  summaryPositivePlotReactive <- reactive({
    req(summariseInputDataReactive())
    summaryPositivePlot(summariseInputDataReactive())
  })
  
  output$summaryPositivePlotOutput <- renderPlot(summaryPositivePlotReactive())
  
  # Plot summary of positive plots serology
  summaryPositivePlotSerologyReactive <- reactive({
    req(summariseInputDataSerologyReactive())
    summaryPositiveSerologyPlot(summariseInputDataSerologyReactive())
  })
  
  output$summaryPositivePlotSerologyOutput <- renderPlot(summaryPositivePlotSerologyReactive())
  
  ##################################################################################################################
  # Neu 05.02.2021: Plot von neuer Abbildung für Lagebericht
  
  # Funktionen hier definieren, da ansonsten die Umlaute nicht funktionieren
  lageberichtPlotFunction <- function(summariseInputDataReactive, summaryTidyDataReachCombined, lowerdate, breakdate){
    if (is.null(summariseInputDataReactive) | is.null(summaryTidyDataReachCombined) ){
      return(NULL)}
    lageberichtFigData <-
      summariseInputDataReactive %>%
      filter(labsClass == "Combined") %>%
      full_join(summaryTidyDataReachCombined, by = c("week" = "week", "year" = "year")) %>%
      select(year, week, sumTests, sumLabs, ratioPositive, sumWeekCapacity, sumLabsCapacity) %>%
      # End 05.01.2020: Fix year bug
      mutate(year_week =if_else(week != 53, 
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53, 
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      mutate(kalenderwoche = paste(week, year, sep = "/")) %>%
      arrange(year_week) %>%
      # Neu Plot mit Darstellung nur jeder zweiten KW
      mutate(kalenderwoche=factor(kalenderwoche, levels=c(kalenderwoche, "\t"))) %>%
      filter(year_week >= lowerdate)
    
    # Erstellen des Plots
    plotLagebericht <-
      lageberichtFigData %>%
      #arrange(year_week) %>%
      #mutate(kalenderwoche=factor(kalenderwoche, levels=kalenderwoche)) %>%
      pivot_longer(cols = c("sumTests", "sumWeekCapacity")) %>%
      #filter(year_week >= lowerdate) %>%
      ggplot(mapping = aes(x = (kalenderwoche), y = value, fill = name)) + 
      geom_col(position = position_dodge()) +
      geom_point(aes(y = ratioPositive*65000), fill = "red",color = "red", size = 4, shape = 18) +
      #  geom_line(aes(y = ratioPositive*130000), fill = "red",color = "red") +
      scale_fill_manual(values = c("#0072bc",
                                   "lightgrey"),
                        labels = c("Durchgeführte Tests",
                                   "Testkapazitäten"))+
      labs(fill = "") +
      #  guides(fill = FALSE, size = FALSE) +
      #  scale_x_datetime(name = "KW", date_breaks = "1 week", 
      #                   expand = c(0,0), date_labels = "%W/%Y") +
      scale_x_discrete(name = "KW", expand = c(0,0), 
                       breaks = lageberichtFigData$kalenderwoche[seq(1, length(lageberichtFigData$kalenderwoche), by = 4)])+
      # scale_y_continuous(name = "Anzahl Tests", expand = c(0,10000),labels = comma_format(big.mark = ".",
      #                                                                                     decimal.mark = ","),
      #                    sec.axis = sec_axis(~./130000, name="Positivenanteil (%)", breaks = c(0,2,4,6,8,10,12,14,16))) +
      # Neu: 23.11.2021 - Anpassung der y-Achse
      scale_y_continuous(name = "Anzahl Tests", expand = c(0,50000), breaks = c(0,500000, 1000000, 1500000, 2000000, 2500000, 3000000),
                         labels = comma_format(big.mark = ".",
                                               decimal.mark = ","),
                         sec.axis = sec_axis(~./65000, name="Positivenanteil (%)", breaks = c(0,5,10,15,20,25,30,35,40,45,50,55))) +
      theme_classic() +
      theme(text = element_text(size = 18))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.title.y.right = element_text(color = "red"),
            axis.text.y.right = element_text(color = "red")) +
      theme(legend.position='top', 
            legend.justification='left',
            legend.direction='horizontal')
    
    return(plotLagebericht)
  }
  
  
  
  lageberichtPlotExportFunction <- function(summariseInputDataReactive, summaryTidyDataReachCombined, lowerdate, breakdate){
    if (is.null(summariseInputDataReactive) | is.null(summaryTidyDataReachCombined) ){
      return(NULL)}
    lageberichtFigData <-
      summariseInputDataReactive %>%
      filter(labsClass == "Combined") %>%
      full_join(summaryTidyDataReachCombined, by = c("week" = "week", "year" = "year")) %>%
      select(year, week, sumTests, sumLabs, ratioPositive, sumWeekCapacity, sumLabsCapacity) %>%
      # End 05.01.2020: Fix year bug
      mutate(year_week =if_else(week != 53, 
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53, 
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      mutate(kalenderwoche = paste(week, year, sep = "/")) %>%
      arrange(year_week) %>%
      mutate(kalenderwoche=factor(kalenderwoche, levels=kalenderwoche)) %>%
      #  pivot_longer(cols = c("sumTests", "sumWeekCapacity")) %>%
      filter(year_week >= lowerdate)
    
    # Erstellen des Plots
    plotLagebericht <-
      lageberichtFigData %>%
      #  arrange(year_week) %>%
      #  mutate(kalenderwoche=factor(kalenderwoche, levels=kalenderwoche)) %>%
      pivot_longer(cols = c("sumTests", "sumWeekCapacity")) %>%
      #  filter(year_week >= lowerdate) %>%
      ggplot(mapping = aes(x = (kalenderwoche), y = value, fill = name)) + 
      geom_col(position = position_dodge()) +
      geom_point(aes(y = ratioPositive*65000), fill = "red",color = "red", size = 4, shape = 18) +
      #  geom_line(aes(y = ratioPositive*130000), fill = "red",color = "red") +
      scale_fill_manual(values = c("#0072bc",
                                   "lightgrey"),
                        labels = c("Durchgeführte Tests",
                                   "Testkapazitäten"))+
      labs(fill = "") +
      #  guides(fill = FALSE, size = FALSE) +
      #  scale_x_datetime(name = "KW", date_breaks = "1 week", 
      #                   expand = c(0,0), date_labels = "%W/%Y") +
      scale_x_discrete(name = "KW", expand = c(0,0),
                       breaks = lageberichtFigData$kalenderwoche[seq(1, length(lageberichtFigData$kalenderwoche), by = 4)])+
      # scale_y_continuous(name = "Anzahl Tests", expand = c(0,10000),labels = comma_format(big.mark = ".",
      #                                                                                     decimal.mark = ","),
      #                    sec.axis = sec_axis(~./130000, name="Positivenanteil (%)", breaks = c(0,2,4,6,8,10,12,14,16))) +
      # Neu: 23.11.2021 - Anpassung der y-Achse
      scale_y_continuous(name = "Anzahl Tests", expand = c(0,50000), breaks = c(0,500000, 1000000, 1500000, 2000000, 2500000, 3000000),
                         labels = comma_format(big.mark = ".",
                                               decimal.mark = ","),
                         sec.axis = sec_axis(~./65000, name="Positivenanteil (%)", breaks = c(0,5,10,15,20,25,30,35,40,45,50,55))) +
      theme_classic() +
      theme(text = element_text(size = 22))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.title.y.right = element_text(color = "red"),
            axis.text.y.right = element_text(color = "red")) +
      theme(legend.position='top', 
            legend.justification='left',
            legend.direction='horizontal',
            legend.text=element_text(size=22))
    
    
    # Erzeugen der Tabelle mit den Meldungen Unterhalb des Plots, aufgeteilt in zwei Teiltabellen
    # Hier dann das Datum anpassen, sonst wird die untere Tabelle immer länger (könnte man auch über length regeln)
    tableLabssumLabs1 <-
      lageberichtFigData %>%
      filter(year_week >= lowerdate & year_week <= breakdate) %>%
      select(kalenderwoche, sumLabs) %>%
      pivot_wider(names_from = kalenderwoche, values_from = c("sumLabs")) 
    
    tableLabssumLabsCapacity1 <-
      lageberichtFigData %>%
      filter(year_week >= lowerdate & year_week <= breakdate) %>%
      select(kalenderwoche, sumLabsCapacity) %>%
      pivot_wider(names_from = kalenderwoche, values_from = c("sumLabsCapacity")) 
    
    tableLabs1 <-
      tableLabssumLabs1 %>%
      add_row(tableLabssumLabsCapacity1)
    
    tablethemelagebericht <- ttheme_default(colhead=list(fg_params = list(rot=90)),
                                            base_size = 12)
    tbl1lagebericht <- tableGrob(tableLabs1, rows=c("Testzahlen", "Testkapazitäten"), theme=tablethemelagebericht)
    # Add Title
    h1lagebericht <- grobHeight(tbl1lagebericht)
    titlelagebericht <- textGrob("Anzahl der Labore, die Daten übermittelt haben", y=unit(0.5,"npc") + 0.5*h1lagebericht, 
                                 vjust=-2.5, gp=gpar(fontsize=18))
    gt1lagebericht <- gTree(children=gList(tbl1lagebericht, titlelagebericht))
    #grid.draw(gt1lagebericht)
    
    
    # Erzeugen des unteren Teils der Tabelle mit den späteren Daten
    tableLabssumLabs2 <-
      lageberichtFigData %>%
      filter(year_week > breakdate) %>%
      select(kalenderwoche, sumLabs) %>%
      pivot_wider(names_from = kalenderwoche, values_from = c("sumLabs")) 
    
    tableLabssumLabsCapacity2 <-
      lageberichtFigData %>%
      filter(year_week > breakdate) %>%
      select(kalenderwoche, sumLabsCapacity) %>%
      pivot_wider(names_from = kalenderwoche, values_from = c("sumLabsCapacity")) 
    
    tableLabs2 <-
      tableLabssumLabs2 %>%
      add_row(tableLabssumLabsCapacity2)
    tbl2lagebericht <- tableGrob(tableLabs2, rows=c("Testzahlen", "Testkapazitäten"), theme=tablethemelagebericht)
    
    # Plot chart and table into one object
    LabsizeTablePlot <-
      arrangeGrob(plotLagebericht , gt1lagebericht, tbl2lagebericht,
                  nrow=3,
                  as.table=TRUE,
                  heights=c(4,1.5,1))
    return(plotLagebericht)
  }
  
  
  
  
  LabsizeTablePlot <- reactive({
    req(summariseInputDataReactive(),summaryTidyDataReachCombined)
    lageberichtPlotFunction(summariseInputDataReactive = summariseInputDataReactive(),
                            summaryTidyDataReachCombined = summaryTidyDataReachCombined(),
                            lowerdate = "2020-02-26 01:00:00",
                            breakdate = "2020-10-21 02:00:00 CEST")
  })
  
  output$LabsizeTablePlotOutput <- renderPlot(LabsizeTablePlot())
  
  
  
  LabsizeTablePlotExport <- reactive({
    req(summariseInputDataReactive(),summaryTidyDataReachCombined)
    lageberichtPlotExportFunction(summariseInputDataReactive = summariseInputDataReactive(),
                                  summaryTidyDataReachCombined = summaryTidyDataReachCombined(),
                                  lowerdate = "2020-02-26 01:00:00",
                                  breakdate = "2020-10-21 02:00:00 CEST")
  })
  ##################################################################################################################
  # Neu 05.02.2021 Plot Probenrückstau
  plotBacklogFunction <- function(tabelle3){
    plotTabelle3 <-
      tabelle3 %>%
      arrange(year_week) %>%
      filter(year_week > as.POSIXct("2020-06-03 02:00:00")) %>%
      mutate(KW=factor(KW, levels=KW)) %>%
      ggplot(mapping = aes(x = KW, y = `Probenrückstau`)) + 
      geom_col(position = position_dodge(), fill = "#0072bc") +
      #  geom_line(aes(y = ratioPositive*130000), fill = "red",color = "red") +
      # scale_fill_manual(values = c("#0072bc"))+
      #  labs(fill = "") +
      #  guides(fill = FALSE, size = FALSE) +
      scale_x_discrete(name = "KW", expand = c(0,0))+
      scale_y_continuous(name = "Anzahl Proben im Rückstau", expand = c(0,1000),labels = comma_format(big.mark = ".",
                                                                                                      decimal.mark = ",")) +
      theme_classic() +
      theme(text = element_text(size = 12))+
      #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Probenrückstau")
    
    tableLabs1backlog <-
      tabelle3 %>%
      filter(year_week > as.POSIXct("2020-06-03 02:00:00")) %>%
      # filter(year_week >= "2020-02-26 01:00:00") %>% #& year_week <= "2020-07-22 02:00:00 CEST") %>%
      select(KW, `Labore mit Rückstau`) %>%
      pivot_wider(names_from = KW, values_from = c("Labore mit Rückstau"))
    
    tt <- ttheme_default(colhead=list(fg_params = list(rot=90)),
                         base_size = 10)
    tbl1 <- tableGrob(tableLabs1backlog, rows=c(""), theme=tt)
    # Add Title
    table1 <-tbl1
    grid.newpage()
    h1 <- grobHeight(table1)
    w1 <- grobWidth(table1)
    title <- textGrob("Anzahl der Labore, die Daten übermittelt haben", y=unit(0.5,"npc") + 0.5*h1, 
                      vjust=-2.5, gp=gpar(fontsize=18))
    #footnote <- textGrob("footnote", 
    #                     x=unit(0.5,"npc") - 0.5*w,
    #                     y=unit(0.5,"npc") - 0.5*h, 
    #                     vjust=1, hjust=0,gp=gpar( fontface="italic"))
    gt1 <- gTree(children=gList(table1, title))
    grid.draw(gt1)
    
    Table3Plot <-
      arrangeGrob(plotTabelle3 , gt1,
                  nrow=2,
                  as.table=TRUE,
                  heights=c(4,1.5))
    return(Table3Plot)
  }
  
  plotBacklogExport <- reactive({
    req(tabelle3())
    plotBacklogFunction(tabelle3())
  })
  
  
  
  
  
  ##################################################################################################################
  # Analysis for individual labs from VOXCO and ALM
  # Calculate QC for VOXCO_ALM single labs
  dataInputQCReactive <- reactive({
    req(fileALMVOXCOPCR_data_inputRaw())
    dataInputQC(fileALMVOXCOPCR_data_inputRaw())
  })
  # Für den ersten Plot: Die bereinigten Daten verwenden, um die Auswirkung des Cutoffs zu sehen
  dataInputQCReactiveCutoff <- reactive({
    req(fileALMVOXCOPCR_data_input())
    dataInputQC(fileALMVOXCOPCR_data_input())
  })
  
  # Plot QC plot of number of tests performed per Lab
  # Plot single labs positive ratio
  individualQCRatioPlotReactive <- reactive({
    req(dataInputQCReactiveCutoff())
    individualQCRatioPlot(dataInputQCReactiveCutoff())
  })
  
  output$individualQCRatioPlotOutput <- renderPlot(individualQCRatioPlotReactive())
  
  ## New 30.09.2020 Interaktiver Plotly Plot für die Darstellung der Positivraten der Labore
  individualQCRatioPlotlyReactive <- reactive({
    req(dataInputQCReactive())
    individualQCRatioPlotly(dataInputQCReactive())
  })
  
  output$individualQCRatioPlotlyOutput <- renderPlotly(individualQCRatioPlotlyReactive())
  
  
  # Plot weekly number of tests by instituions: Several Plots needed for visual inspection of plausibility
  # To Do: find better suited method for plotting
  numberTestPlotAllReactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlotAll(dataInputQCReactive())
  })
  output$numberTestPlotAllOutput <- renderPlot(numberTestPlotAllReactive())
  
  
  numberTestPlot50Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 1, 50)
  })
  output$numberTestPlot50Output <- renderPlot(numberTestPlot50Reactive())
  
  
  numberTestPlot75Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 50, 76)
  })
  output$numberTestPlot75Output <- renderPlot(numberTestPlot75Reactive())
  
  
  numberTestPlot104Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 76, 104)
  })
  output$numberTestPlot104Output <- renderPlot(numberTestPlot104Reactive())
  
  
  numberTestPlot133Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 104, 133)
  })
  output$numberTestPlot133Output <- renderPlot(numberTestPlot133Reactive())
  
  
  numberTestPlot157Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 133, 157)
  })
  output$numberTestPlot157Output <- renderPlot(numberTestPlot157Reactive())
  
  numberTestPlot191Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 157, 191)
  })
  output$numberTestPlot191Output <- renderPlot(numberTestPlot191Reactive())
  
  
  numberTestPlot230Reactive <- reactive({
    req(dataInputQCReactive())
    numberTestPlot(dataInputQCReactive(), 191, 400)
  })
  output$numberTestPlot230Output <- renderPlot(numberTestPlot230Reactive())
  ########################################################################################################
  
  
  
  
  
  ########################################################################################################
  ########################################################################################################
  # Output of QC datatables
  QCfailedSingleReactive <- reactive({
    req(dataInputQCReactive())
    QCfailedSingle(dataInputQCReactive())
  })
  
  output$QCfailedSingleOutput <- renderDT(QCfailedSingleReactive() %>%
                                            mutate(QCratioPositive = if_else(QCratioPositive == TRUE, "Ja", "Nein")) %>%
                                            mutate(QCpatients = if_else(QCpatients == TRUE, "Ja", "Nein")) %>%
                                            mutate(QCpositiveTestReported = if_else(QCpositiveTestReported == TRUE, "Ja", "Nein")) %>%
                                            select("Jahr" = year,
                                                   "KW" = week,
                                                   "RKI-ID" = institute,
                                                   "Datenquelle" = source,
                                                   "Verhältnis positiver Test Ok?" = QCratioPositive,
                                                   "Mehr Tests als Patienten?" = QCpatients,
                                                   "Positive Testergebnisse gemeldet?" = QCpositiveTestReported,
                                                   "Anzahl Tests" = tests,
                                                   "Anzahl Patienten" = patients,
                                                   "Anzahl negativer Tests" = negative,
                                                   "Anzahl positiver Tests" = positive,
                                                   "Verhältnis positiver Tests" = ratioPositive),
                                          options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  
  outliersRatioPositiveReactive <- reactive({
    req(dataInputQCReactive())
    outliersRatioPositive(dataInputQCReactive(), cutoff())
  })
  
  output$outliersRatioPositiveOutput <- renderDT(outliersRatioPositiveReactive() %>%
                                                   mutate(QCratioPositive = if_else(QCratioPositive == TRUE, "Ja", "Nein")) %>%
                                                   mutate(QCpatients = if_else(QCpatients == TRUE, "Ja", "Nein")) %>%
                                                   mutate(QCpositiveTestReported = if_else(QCpositiveTestReported == TRUE, "Ja", "Nein")) %>%
                                                   select("Jahr" = year,
                                                          "KW" = week,
                                                          "RKI-ID" = institute,
                                                          "Datenquelle" = source,
                                                          "Verhältnis positiver Test Ok?" = QCratioPositive,
                                                          "Mehr Tests als Patienten?" = QCpatients,
                                                          "Positive Testergebnisse gemeldet?" = QCpositiveTestReported,
                                                          "Anzahl Tests" = tests,
                                                          "Anzahl Patienten" = patients,
                                                          "Anzahl negativer Tests" = negative,
                                                          "Anzahl positiver Tests" = positive,
                                                          "Verhältnis positiver Tests" = ratioPositive),
                                                 options = list(scrollX = TRUE, pageLength = 50), filter = "top")
  ########################################################################################################
  # Neu 05.02.2021: Ausgabe einer finalen Exceltabelle mit allen Formaten 
  # Dev Neu 04.02.2021: Output der Exeltabelle im finalen Format
  tabelle1 <- reactive({
    #  if (is.null(summariseInputDataReactive())){
    #    return(NULL)}
    summariseInputDataReactive() %>%
      filter(labsClass == "Combined") %>%
      mutate(Kalenderwoche = case_when((week <= 10 & year == 2020) ~ "Bis einschließlich KW10, 2020",
                                       (week > 10 & week <= 45 & year == 2020) ~ paste(week, year, sep = "/"),
                                       TRUE ~ paste(week, "/", year, sep = "") ))%>%
      group_by(Kalenderwoche) %>%
      summarise("Anzahl Testungen" = sum(sumTests),
                "Positiv getestet" = sum(sumPositive),
                "Positivenanteil (%)" = (sum(sumPositive)/sum(sumTests))*100,
                "Anzahl übermittelnder Labore" = sum(sumLabs),
                year_week = year_week) %>%
      arrange(year_week) %>%
      mutate(`Positivenanteil (%)` = if_else(Kalenderwoche == "Bis einschließlich KW10, 2020", NA_real_, `Positivenanteil (%)`),
             `Anzahl übermittelnder Labore` = if_else(Kalenderwoche == "Bis einschließlich KW10, 2020", NA_real_, `Anzahl übermittelnder Labore`),
             year_week = if_else(Kalenderwoche == "Bis einschließlich KW10, 2020", min(year_week), year_week)) %>%
      unique() %>%
      adorn_totals("row", name = "Summe") %>%
      mutate(`Positivenanteil (%)` = if_else(Kalenderwoche == "Summe", NA_real_, `Positivenanteil (%)`),
             `Anzahl übermittelnder Labore` = if_else(Kalenderwoche == "Summe", NA_real_, `Anzahl übermittelnder Labore`))
  })
  
  
  tabelle2 <- reactive({
    # if (is.null(summaryTidyDataReachCombined()) ){
    #    return(NULL)}
    summaryTidyDataReachCombined() %>%
      ungroup() %>%
      select(year, week, sumLabsCapacity, sumDay, sumWeekWorkdays, sumWeekCapacity) %>%
      # End 05.01.2020: Fix year bug
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      arrange(year_week) %>%
      filter(year_week >= "2020-03-11 01:00:00") %>%
      mutate(`KW, für die die Angabe prognostisch erfolgt ist` = paste(week,year,sep = "/")) %>%
      select(`KW, für die die Angabe prognostisch erfolgt ist`,
             `Anzahl übermittelnde Labore` = sumLabsCapacity,
             `Testkapazität pro Tag` = sumDay,
             `Theoretische wöchentliche Kapazität anhand von Wochenarbeitstagen` = sumWeekWorkdays,
             `Reale Testkapazität zum Zeitpunkt der Abfrage` = sumWeekCapacity)
  })
  
  
  tabelle3 <- reactive({
    # if (is.null(summaryTidyDataBacklogReactive()) ){
    #    return(NULL)}
    summaryTidyDataBacklogReactive() %>%
      ungroup()%>%
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      arrange(year_week) %>%
      filter(year_week >= "2020-04-08 02:00:00") %>%
      select(`Labore mit Rückstau` = sumLabszero,
             KW =  week,
             `Probenrückstau` = sumBacklog,
             year_week)
  })
  
  
  shortenTabelle1Function <- function(tabelle1In){
    tabelle1neu <- tabelle1In
    tabelle1neu$ID <- seq.int(nrow(tabelle1neu))
    tabelle1neu %>%
      ungroup() %>%
      filter(Kalenderwoche != "Summe") %>%
      mutate(Kalenderwoche = case_when((ID <= nrow(tabelle1neu)-11 ) ~ paste("Bis einschließlich KW", tabelle1neu$Kalenderwoche[nrow(tabelle1neu)-11],sep = ""),
                                       TRUE ~ Kalenderwoche))%>%
      group_by(Kalenderwoche) %>%
      summarise("Anzahl Testungen" = sum(`Anzahl Testungen`),
                "Positiv getestet" = sum(`Positiv getestet`),
                "Positivenanteil (%)" = (sum(`Positiv getestet`)/sum(`Anzahl Testungen`))*100,
                "Anzahl übermittelnder Labore" = sum(`Anzahl übermittelnder Labore`),
                year_week = year_week) %>%
      arrange(year_week) %>%
      mutate(`Positivenanteil (%)` = if_else(Kalenderwoche == paste("Bis einschließlich KW", tabelle1neu$Kalenderwoche[nrow(tabelle1neu)-11],sep = ""), NA_real_, `Positivenanteil (%)`),
             `Anzahl übermittelnder Labore` = if_else(Kalenderwoche == paste("Bis einschließlich KW", tabelle1neu$Kalenderwoche[nrow(tabelle1neu)-11],sep = ""), NA_real_, `Anzahl übermittelnder Labore`),
             year_week = if_else(Kalenderwoche == paste("Bis einschließlich KW", tabelle1neu$Kalenderwoche[nrow(tabelle1neu)-11],sep = ""), min(year_week), year_week)) %>%
      unique() %>%
      adorn_totals("row", name = "Summe") %>%
      mutate(`Positivenanteil (%)` = if_else(Kalenderwoche == "Summe", NA_real_, `Positivenanteil (%)`),
             `Anzahl übermittelnder Labore` = if_else(Kalenderwoche == "Summe", NA_real_, `Anzahl übermittelnder Labore`))
  }
  
  
  tabelle1kurz <- reactive({
    req(tabelle1())
    shortenTabelle1Function(tabelle1())
  })
  
  tabelle2kurz <- reactive({
    tabelle2() %>%
      tail(12)
  })
  
  erlaeuterungen <- c("Das RKI erfasst wöchentlich die SARS-CoV-2-Testzahlen. Hierfür werden deutschlandweit Daten von Universitätskliniken, Forschungseinrichtungen sowie klinischen und ambulanten Laboren zusammengeführt. Die Erfassung basiert auf einer freiwilligen Mitteilung der Labore und erfolgt über eine webbasierte Plattform (VOXCO, RKI-Testlaborabfrage) oder in Zusammenarbeit mit der am RKI etablierten, laborbasierten SARS-CoV-2-Surveillance (eine Erweiterung der Antibiotika-Resistenz-Surveillance, ARS), dem Netzwerk für respiratorische Viren (RespVir) sowie der Abfrage eines labormedizinischen Berufsverbands. Die Erfassung liefert Hinweise zur aktuellen Situation in den Laboren, erlaubt aber keine detaillierten Auswertungen oder Vergleiche mit den gemeldeten Fallzahlen.") %>%
    as_tibble()
  ########################################################################################################
  
  
  #### New 2022-09-09 Generate open date output ####
  generateOpenData <- function(summariseInputDataReactive, summaryTidyDataReachCombined,
                               summaryTidyDataBacklogReactive){
    tabelle1csv <-
      summariseInputDataReactive %>%
      filter(labsClass == "Combined") %>%
      mutate(date = case_when((week <= 10 & year == 2020) ~ "2020-W10",
                              TRUE ~ paste(year, "-W", week, sep = "")))%>%
      group_by(date) %>%
      summarise(tests_total = sum(sumTests),
                tests_positive = sum(sumPositive),
                tests_positive_ratio = round(sum(sumPositive)/sum(sumTests), digits = 4),
                laboratories_tests = sum(sumLabs),
                year_week = year_week) %>%
      arrange(year_week) %>%
      mutate(laboratories_tests = if_else(date == "2020-W10", NA_real_, laboratories_tests),
             year_week = if_else(date == "2020-W10", min(year_week), year_week)) %>%
      unique() %>%
      ungroup() %>%
      mutate(tests_positive_accumulated = cumsum(ifelse(is.na(tests_positive), 0, tests_positive)),
             tests_positive_accumulated = if_else(tests_positive_accumulated == 0, NA_real_,
                                                  tests_positive_accumulated),
             tests_total_accumulated = cumsum(tests_total)) %>%
      dplyr::select(date, tests_total, tests_total_accumulated, tests_positive,
                    tests_positive_accumulated,
                    tests_positive_ratio, laboratories_tests)
    
    
    tabelle2csv <- summaryTidyDataReachCombined %>%
      ungroup() %>%
      dplyr::select(year, week, sumLabsCapacity, sumDay, sumWeekWorkdays, sumWeekCapacity) %>%
      # End 05.01.2020: Fix year bug
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      arrange(year_week) %>%
      filter(year_week >= "2020-03-11 01:00:00") %>%
      mutate(date = case_when(week-1 != 0 ~ paste(year, "-W", week-1, sep = ""),
                              week-1 == 0 & year == 2021 ~ paste(year-1, "-W", 53, sep = ""),
                              week-1 == 0 & year == 2022 ~ paste(year-1, "-W", 52, sep = ""))) %>%
      dplyr::select(date,
                    capacities_daily = sumDay,
                    capacities_weekly_theoretically = sumWeekWorkdays,
                    capacities_weeklyweek_actually = sumWeekCapacity,
                    laboratories_capacities = sumLabsCapacity)
    
    tabelle3csv <- summaryTidyDataBacklogReactive %>%
      ungroup()%>%
      mutate(year_week =if_else(week != 53,
                                paste(year, week, 2, sep="/"),
                                paste((year + 1), 1, 2, sep="/"))) %>%
      mutate(year_week = parse_date_time(year_week,'Y/W/u'))%>%
      # new 12.01.2021: Einfügen von KW53 für das Plotting 
      mutate(year_week = if_else(week == 53,
                                 as.POSIXct("2020-12-31  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-31  01:00:00"), as.POSIXct("2020-12-29  01:00:00"),
                                 if_else(year(year_week) == 2020, 
                                         year_week -  as.difftime(6, unit="days"), year_week))) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-30 01:00:00 CET"), as.POSIXct("2021-01-06  01:00:00"),
                                 year_week)) %>%
      mutate(year_week = if_else(year_week == as.POSIXct("2020-12-29 01:00:00 CET"), as.POSIXct("2020-12-30 01:00:00 CET"),
                                 year_week)) %>%
      arrange(year_week) %>%
      filter(year_week >= "2020-04-08 02:00:00") %>%
      mutate(date = paste(year, "-W", week, sep = "")) %>%
      dplyr::select(date,
                    laboratories_samplebacklog = sumLabszero,
                    samplebacklog = sumBacklog) %>%
      mutate(laboratories_samplebacklog = if_else(laboratories_samplebacklog == 0,
                                                  NA_integer_, as.integer(laboratories_samplebacklog)),
             samplebacklog = if_else(samplebacklog == 0, NA_integer_, as.integer(samplebacklog)))
    
    tabellecsv <-
      tabelle1csv %>%
      left_join(tabelle2csv, by = "date") %>%
      left_join(tabelle3csv, by = "date")
    return(tabellecsv)
  }
  
  
  tabellecsv <- reactive({
    req(summariseInputDataReactive(),
        summaryTidyDataReachCombined(),
        summaryTidyDataBacklogReactive())
    generateOpenData(summariseInputDataReactive(),
                     summaryTidyDataReachCombined(),
                     summaryTidyDataBacklogReactive())
  })
  
  
  
  
  
  
  
  ########################################################################################################
  ########################################################################################################
  # Generate downloadable result files
  # Neu 05.02.2021: Finale Exceltabelle
  output$downloadTestzahlerfassung <- downloadHandler(
    filename = function() {
      paste("Testzahlerfassung_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(as.data.frame(erlaeuterungen), file,
                 sheetName = "0_Erläuterungen", append = FALSE, row.names= FALSE)
      write.xlsx(as.data.frame(tabelle1()) %>%
                   select(-year_week), file,
                 sheetName = "1_Testzahlerfassung", append = TRUE, row.names= FALSE)
      write.xlsx(as.data.frame(tabelle1kurz()) %>%
                   select(-year_week) %>%
                   mutate(Kalenderwoche = str_remove(Kalenderwoche, "\\*")), file,
                 sheetName = "1_Testzahlerfassung kurz", append = TRUE, row.names= FALSE)
      write.xlsx(as.data.frame(tabelle2()) %>%
                   mutate(`Testkapazität pro Tag` = as.integer(`Testkapazität pro Tag`)) %>%
                   mutate(`Reale Testkapazität zum Zeitpunkt der Abfrage` = as.integer(`Reale Testkapazität zum Zeitpunkt der Abfrage`)),
                 file,
                 sheetName = "2_Testkapazitäten", append = TRUE, row.names= FALSE)
      write.xlsx(as.data.frame(tabelle2kurz()) %>%
                   mutate(`Testkapazität pro Tag` = as.integer(`Testkapazität pro Tag`)) %>%
                   mutate(`Reale Testkapazität zum Zeitpunkt der Abfrage` = as.integer(`Reale Testkapazität zum Zeitpunkt der Abfrage`)),
                 file,
                 sheetName = "2_Testkapazitäten kurz", append = TRUE, row.names= FALSE)
      write.xlsx(as.data.frame(tabelle3()) %>%
                   select(-year_week) %>%
                   mutate(`Probenrückstau` = as.integer(`Probenrückstau`)), file,
                 sheetName = "3_Probenrückstau", append = TRUE, row.names= FALSE)
    }
  )
  
  
  # Neu 09.09.2022: Open Data
  output$downloadTestzahlerfassungOD <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2-PCR-Testungen_in_Deutschland", ".csv", sep="")
    },
    content = function(file) {
      write.csv(tabellecsv(),file,
                fileEncoding = "UTF-8",
                row.names = FALSE,
                quote = FALSE)
    }
  )
  
  output$downloadTestzahlerfassungODxlsx <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2-PCR-Testungen_in_Deutschland", ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(as.data.frame(tabellecsv()),file,
                sheetName = "Testzahlen", append = FALSE, row.names= FALSE)
    }
  )
  
  
  # QC results
  output$downloadQCSummary <- downloadHandler(
    filename = function() {
      #week(max(fileARSRespVirSerology_data_input$year_week))
      #year(max(fileARSRespVirSerology_data_input$year_week))
      paste("QCfailedSummary_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(QCfailedSingleReactive() %>%
                   mutate(QCratioPositive = if_else(QCratioPositive == TRUE, "Ja", "Nein")) %>%
                   mutate(QCpatients = if_else(QCpatients == TRUE, "Ja", "Nein")) %>%
                   mutate(QCpositiveTestReported = if_else(QCpositiveTestReported == TRUE, "Ja", "Nein")) %>%
                   select("Jahr" = year,
                          "KW" = week,
                          "RKI-ID" = institute,
                          "Datenquelle" = source,
                          "Verhältnis positiver Test Ok?" = QCratioPositive,
                          "Mehr Tests als Patienten?" = QCpatients,
                          "Positive Testergebnisse gemeldet?" = QCpositiveTestReported,
                          "Anzahl Tests" = tests,
                          "Anzahl Patienten" = patients,
                          "Anzahl negativer Tests" = negative,
                          "Anzahl positiver Tests" = positive,
                          "Verhältnis positiver Tests" = ratioPositive)
                 , file,
                 sheetName = "QCfailedIndividualReactive", append = TRUE)
      write.xlsx(outliersRatioPositiveReactive() %>%
                   mutate(QCratioPositive = if_else(QCratioPositive == TRUE, "Ja", "Nein")) %>%
                   mutate(QCpatients = if_else(QCpatients == TRUE, "Ja", "Nein")) %>%
                   mutate(QCpositiveTestReported = if_else(QCpositiveTestReported == TRUE, "Ja", "Nein")) %>%
                   select("Jahr" = year,
                          "KW" = week,
                          "RKI-ID" = institute,
                          "Datenquelle" = source,
                          "Verhältnis positiver Test Ok?" = QCratioPositive,
                          "Mehr Tests als Patienten?" = QCpatients,
                          "Positive Testergebnisse gemeldet?" = QCpositiveTestReported,
                          "Anzahl Tests" = tests,
                          "Anzahl Patienten" = patients,
                          "Anzahl negativer Tests" = negative,
                          "Anzahl positiver Tests" = positive,
                          "Verhältnis positiver Tests" = ratioPositive), file,
                 sheetName = "QCfailedRatio", append = TRUE)
    }
  )
  
  # Summary of results
  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste("SummaryResults_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(as.data.frame(summariseInputDataReactive()) %>%
                   filter(labsClass == "Combined") %>%
                   select("Jahr" = year,
                          "KW" = week, 
                          "Anzahl Testungen" = sumTests,
                          "Positiv getestet" = sumPositive,
                          "Positivquote (%)" = ratioPositive,
                          "Anzahl übermittelnde Labore" = sumLabs), file,
                 sheetName = "SummaryPCRReactive", append = FALSE)
      write.xlsx(as.data.frame(summariseInputDataSerologyReactive()) %>%
                   mutate(sumTests = as.integer(sumTests),
                          sumPositive = as.integer(sumPositive)) %>%
                   select("Jahr" = year,
                          "KW" = week, 
                          "Anzahl Testungen" = sumTests,
                          "Positiv getestet" = sumPositive,
                          "Positivquote (%)" = ratioPositive,
                          "Anzahl übermittelnde Labore" = sumLabs), file,
                 sheetName = "SummarySerology", append = TRUE)
      write.xlsx(as.data.frame(summariseInputDataAntigenReactive()) %>%
                   select("Jahr" = year,
                          "KW" = week,
                          "Anzahl Testungen" = sumTests,
                          "Positiv getestet" = sumPositive,
                          "Davon per PCR bestätigt" = sumPCRconfirmed,
                          "Positivquote (%)" = ratioPositive,
                          "Anzahl übermittelnde Labore" = sumLabs), file,
                 sheetName = "SummaryAntigen", append = TRUE)
      write.xlsx(as.data.frame(summaryTidyDataReachCombined()) %>%
                   mutate(sumLabsWorkweek = as.integer(sumLabsWorkweek),
                          sumLabsCapacity = as.integer(sumLabsCapacity),
                          sumDay = as.integer(sumDay),
                          sumWeekWorkdays = as.integer(sumWeekWorkdays),
                          sumWeekCapacityALM = as.integer(sumWeekCapacityALM)) %>%
                   select("Jahr" = year,
                          "KW*, für die die Angabe prognostisch erfolgt ist" = week,
                          "Anzahl übermittelnder Labore (Arbeitstage)" = sumLabsWorkweek,
                          "Anzahl übermittelnder Labore (Reichweite)" = sumLabsCapacity,
                          "Testkapazität pro Tag" = sumDay,
                          "Neu ab KW15: wöchentliche Kapazität anhand von Wochenarbeitstagen" = sumWeekWorkdays,
                          "Reale Kapazität anhand von Testreichweite" = sumWeekCapacity,
                          "Reale Kapazität anhand von Testreichweite und Wochenarbeitstagen (ALM Kriterien)" = sumWeekCapacityALM
                   ), file,
                 sheetName = "SummaryReachCombined", append = TRUE)
      write.xlsx(as.data.frame(testReachWeekFunctionReactive()), file,
                 sheetName = "SummaryReachDays", append = TRUE)
      write.xlsx(as.data.frame(summaryTidyDataBacklogReactive()) %>%
                   select("Jahr" = year,
                          "KW" = week,
                          "Anzahl übermittelnde Labore" = sumLabszero,
                          "Rückstau" = sumBacklog), file,
                 sheetName = "SummaryBacklog", append = TRUE)
      write.xlsx(as.data.frame(summariseInputCombinedReactive()), file,
                 sheetName = "SummaryPCRCombined", append = TRUE)
      write.xlsx(as.data.frame(fileVOXCOSGenPCR_SingleLabsReactive()), file,
                 sheetName = "S Gen", append = TRUE)
      # write.xlsx(as.data.frame(summariseInputDataSequncingReactive()), file,
      #            sheetName = "Sequenzierung ALM", append = TRUE)
      # write.xlsx(as.data.frame(summariseVOXCOseqReactive()), file,
      #            sheetName = "Sequenzierung VOXCO", append = TRUE)
    }
  )
  
  
  output$downloadSequenzierung <- downloadHandler(
    filename = function() {
      paste("Sequenzierung_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(as.data.frame(summary_sequencing()), file,
                 sheetName = "Sequenzierung", append = FALSE)
      write.xlsx(as.data.frame(summary_sequencing_ALM()), file,
                 sheetName = "Sequenzierung_ALM", append =TRUE)
      write.xlsx(as.data.frame(summary_sequencing_VOXCO()), file,
                 sheetName = "Sequenzierung_VOXCO", append = TRUE)
      # write.xlsx(as.data.frame(sequencinq_ALM_VOXCO_final()), file,
      #             sheetName = "Sequenzierung_Einzeln", append = TRUE)
      
    }
  )
  
  output$downloadSingleLabs <- downloadHandler(
    filename = function() {
      paste("SingleLabsResults_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".xlsx", sep="")
    },
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(week(max(fileALMVOXCOPCR_data_input()$year_week)))){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year_week == max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  output$downloadSingleLabs1 <- downloadHandler(
    filename = c("SingleLabsResults_2020_KW5-14.xlsx"),
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(5:14)){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year == 2020) %>%
                     filter(week == i), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  output$downloadSingleLabs2 <- downloadHandler(
    filename = function() {
      paste("SingleLabsResults_",
            2020, "_KW_",
            15, "-",
            23,
            ".xlsx", sep="")
    },
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(15:23)){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year == 2020)%>%
                     filter(week == i), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  output$downloadSingleLabs3 <- downloadHandler(
    filename = function() {
      paste("SingleLabsResults_",
            2020, "_KW_",
            24, "-",
            33,
            ".xlsx", sep="")
    },
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(24:33)){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year == 2020) %>%
                     filter(week == i), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  
  output$downloadSingleLabs4 <- downloadHandler(
    filename = function() {
      paste("SingleLabsResults_",
            2020, "_KW_",
            34, "-",
            42,
            ".xlsx", sep="")
    },
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(34:42)){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year == 2020)%>%
                     filter(week == i), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  output$downloadSingleLabs5 <- downloadHandler(
    filename = function() {
      paste("SingleLabsResults_",
            2020, "_KW_",
            43, "-",
            53,
            ".xlsx", sep="")
    },
    content = function(file) {
      # write.xlsx(fileALMVOXCOPCR_data_input(), file,
      #            sheetName = "SinglePCR", append = FALSE)
      for(i in c(43:53)){
        write.xlsx(fileALMVOXCOPCR_data_input() %>%
                     filter(year == 2020)%>%
                     filter(week == i), file,
                   sheetName = paste("SinglePCR", i, sep = "_"), append = TRUE)
      }
    }
  )
  
  
  # Generate downloadable plots with all PCR and Serology data
  output$downloadFigPCR <- downloadHandler(
    filename = function() {
      paste("PlotsPCR_",
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".png", sep="")
    },
    content = function(file) {
      ggarrange(summaryQCRatioPlotReactive(),
                summaryLabsPlotReactive(),
                summaryTestsPlotReactive(),
                summaryPositivePlotReactive()
                , ncol = 2, nrow = 2) %>%
        ggexport(filename = file,   width = 960, height = 800)
    }
  )
  
  output$downloadFigSerology <- downloadHandler(
    filename = function() {
      paste("PlotsSerology_", 
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".png", sep="")
    },
    content = function(file) {
      ggarrange(summaryQCRatioPlotSerologyReactive(),
                summaryLabsPlotSerologReactive(),
                summaryTestsPlotSerologyReactive(),
                summaryPositivePlotSerologyReactive()
                ,ncol = 2, nrow = 2) %>%
        ggexport(filename = file, width = 960, height = 800)
    }
  )
  
  output$downloadFigLagebericht <- downloadHandler(
    filename = function() {
      paste("PlotLagebericht_", 
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".png", sep="")
    },
    content = function(file) {
      LabsizeTablePlotExport() %>%
        ggsave(filename = file,
               width = 18, height = 8.7)
    }
  )
  
  
  output$downloadFigBacklog <- downloadHandler(
    filename = function() {
      paste("PlotRueckstau_", 
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".png", sep="")
    },
    content = function(file) {
      plotBacklogExport() %>%
        ggsave(filename = file,
               width = 18, height = 8)
    }
  )
  
  
  output$downloadFigBoxplotQC <- downloadHandler(
    filename = function() {
      paste("PlotBoxplotQC_", 
            year(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), "_KW_",
            week(max(sort(unique(fileALMVOXCOPCR_data_input()$year_week)))), ".png", sep="")
    },
    content = function(file) {
      individualQCRatioPlotReactive() %>%
        ggsave(filename = file,
               width = 15, height = 8)
    }
  )
  
  
  
  
  
}

# 03:       Run shiny app ####

shinyApp(ui = ui, server = server)