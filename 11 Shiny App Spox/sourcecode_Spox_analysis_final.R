##
# Source code for analysis of Spox measurements
# Daniel Stern
# RKI
# Version 1
# 2024-08-03
##

##
# Outline
#
# The following code is used to generate a shiny for automated data analysis
# and plotting for the Spox Study. The app should be used to analyse results
# of the multiplex-assay for Mpox-analysis, show QC-Plot and generate 
# downloadable content with finished analysis and QC information for each
# measurement.
#
# The app is structured in the following segments
# 0) Load packages and prepare workspace
# 1) Load models, data and metadata
# 2) Fit data using the DrLumi package and generate quantified results
# 3) Determine serostatus for IgG and IgM compared to the delta results
# 4) Predict data for classification using the different models for LDA
# 5) Predict data for NT using the lasso regression model
# 6) Generate graphical plots for QC
#     a) Plot standards
#     b) Plot controls (Background and PC)
#     c) Output Report for single serum with drop down menu for selection
# 7) Generate downloadable output
#     a) CSV output with all results
#     b) CSV output with short results
#     c) Word report
##

##
# 0) Load packages and prepare workspace ####
# Clean global environment
rm(list = ls(all.names = TRUE))

# Load libraries
library(rio)
library(lubridate)
library(tidyverse)
library(ggthemes)
#library(drLumi)
library(caret)
library(glmnet)
library(docxtractr)
library(officer)

# Dependencies of drLumi
library(reshape)
library(Hmisc)
library(chron)
library(gdata)
library(plyr)
library(stringr)
library(ggplot2)
library(minpack.lm)
library(stats)
library(irr)
library(rootSolve)
library(msm)



##

##
# 1) Load models, data and metadata ####
# Load models and data
source("functions/function_read_multiplex.R", encoding = "UTF-8") ## Function for data import
source("functions/functions_quantify.R", encoding = "UTF-8") ## Function for data quantification

# Source all  drLumi Files to test, if it works without loading the package
files.sources = list.files("functions/drLumi/")
sapply(paste("functions/drLumi/",files.sources, sep = ""), source)


ecdata <- import("data/EC_input.xlsx") %>%  # Load ec-file for quantification
  pivot_longer(c(-sample), names_to = "analyte", values_to = "ec")
load("data/threshold.Rdata") # Load cutoff values for IgG and IgM 
threshold <- 
  threshold %>% 
  filter(assaytype == "Multiplex")

load("data/LDAmodels.Rdata") # Load LDA models 
load("data/NTLassoReg.Rdata") # Loat NT-Lasso Regressin data
load("data/GLMmodels.RData") # Load GLM Models with all data
load("data/threshold_MaxSpec.RData") # Load Cutoff-Values for 100% specificity




ATI_N_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "ATI.N" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()

Combined_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "Combined" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()

# Load input files
#dataFile1 <- c("input/Pox_25_14_Platte 1 240323.xlsx")
#dataFile2 <- c("input/Pox_25_14_Platte 2 240323.xlsx")
#metaWord <- c("input/Pox 25B Assay 2 Platten 24 Proben - Teil 14 für neue Bead-Konzentration.docx")

#dataFile1 <- c("input/Pox 12 Platte 1 Ergebnis 11.10.22.xlsx")
#dataFile2 <- c("input/Pox 12 Platte 2 Ergebnis 11.10.22.xlsx")
#metaWord <- c("input/Pox 12 Planung.docx")

dataFile1 <- c("input/")
dataFile2 <- c("input/SPox_01_Platte_2.xlsx")
metaWord <- c("input/Protokoll SPox_01.docx")


# Import results from both files and word document with metadata
inputData <- join_data_function(read_plate_function(dataFile1, 1), 
                                read_plate_function(dataFile2, 2),
                                read_metadata_function(metaWord,
                                                       dataFile1, dataFile2))


# Add flags for QC samples
dataInputQCHSA <-
  inputData %>%
  filter(analyte == "HSA", dilution == 100) %>% 
  group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
  dplyr::summarise(highBG = if_else(data > 500, T, F),
                   background = data)

dataInputQCIgG <-
  inputData %>%
  filter(analyte == "ahIgG", dilution == 100, 
         isotype == "IgG") %>% 
  group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
  dplyr::summarise(lowIg = if_else(data < 10000, T, F),
                   positiveIg = data)

dataInputQCIgM <-
  inputData %>%
  filter(analyte == "ahIgM", dilution == 100,
         isotype == "IgM") %>% 
  group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
  dplyr::summarise(lowIg = if_else(data < 6000, T, F),
                   positiveIg = data)

dataInputQC <-
  rbind(dataInputQCIgG, dataInputQCIgM)

# Normalize data
dataInput <-
  inputData %>%
  dplyr::left_join(dataInputQC, by = c("experiment", "plate_assay","sampleID_assay", "isotype",
                                       "dilution")) %>%
  dplyr::left_join(dataInputQCHSA, by = c("experiment", "plate_assay","sampleID_assay", "isotype",
                                          "dilution"))
##

##
# 2) Fit data using drLumi
quantInputIgG_plate_1 <-
  dataInput%>% 
  dplyr::rename(dilution_assay = dilution) %>% 
  filter(isotype == "IgG" & plate_assay == 1)

quantInputIgG_plate_2 <-
  dataInput%>% 
  dplyr::rename(dilution_assay = dilution) %>% 
  filter(isotype == "IgG" & plate_assay == 2)

quantInputIgM_plate_1 <-
  dataInput%>% 
  dplyr::rename(dilution_assay = dilution) %>% 
  filter(isotype == "IgM"& plate_assay == 1)

quantInputIgM_plate_2 <-
  dataInput%>% 
  dplyr::rename(dilution_assay = dilution) %>% 
  filter(isotype == "IgM" & plate_assay == 2)

quantifiedIgG_plate_1 <- quantifyMultiplex(ecdata = ecdata,
                                           multiInputList = quantInputIgG_plate_1)

quantifiedIgG_plate_2 <- quantifyMultiplex(ecdata = ecdata,
                                           multiInputList = quantInputIgG_plate_2)

quantifiedIgM_plate_1 <- quantifyMultiplexReFit(ecdata = ecdata,
                                                multiInputList = quantInputIgM_plate_1,
                                                allanalytesIn = quantifiedIgG_plate_1$standards,
                                                1)

quantifiedIgM_plate_2 <- quantifyMultiplexReFit(ecdata = ecdata,
                                                multiInputList = quantInputIgM_plate_2,
                                                allanalytesIn = quantifiedIgG_plate_1$standards,
                                                1)


standard_plate_1 <- quantifiedIgG_plate_1$standards  
standard_plate_2 <- quantifiedIgG_plate_2$standards 

standard_info_plate_1 <- do.call(rbind, lapply(standard_plate_1, '[[', 2)) %>% 
  add_column(plate_assay = 1, 
             date = unique(dataInput$date), 
             experiment = unique(dataInput$experiment), 
             operator = unique(dataInput$operator), 
             bead_lot = unique(dataInput$bead_lot))
standard_info_plate_2 <- do.call(rbind, lapply(standard_plate_1, '[[', 2)) %>% 
  add_column(plate_assay = 2, 
             date = unique(dataInput$date), 
             experiment = unique(dataInput$experiment), 
             operator = unique(dataInput$operator), 
             bead_lot = unique(dataInput$bead_lot))
standard_info <-
  bind_rows(standard_info_plate_1, standard_info_plate_2)


dataInputQuant <- 
  rbind(quantifiedIgG_plate_1$dataOutput, 
        quantifiedIgG_plate_2$dataOutput, 
        quantifiedIgM_plate_1$dataOutput,
        quantifiedIgM_plate_2$dataOutput)  %>% 
  mutate(analyte = factor(analyte, levels = c("A27L",
                                              "A29", 
                                              "L1R", 
                                              "M1",
                                              "D8L",
                                              "E8",
                                              "H3L",
                                              "A33R", 
                                              "A35R",
                                              "B5R", 
                                              "B6", 
                                              "A5L",
                                              "ATI-C", 
                                              "ATI-N",
                                              "VACV", 
                                              "Delta"), ordered = TRUE)) %>% 
  group_by(analyte, experiment, plate_assay) %>% 
  mutate(dataIn = case_when(is.na(quant_conc) & grepl("above", remarks) ~ max(quant_conc, na.rm = TRUE),
                            is.na(quant_conc) & grepl("below", remarks) ~ min(quant_conc, na.rm = TRUE),
                            TRUE ~ quant_conc),
         dataIn = if_else(dataIn < 1, 1, dataIn)) %>% 
  ungroup()

##
# 03) Determine serostatus for IgG and IgM compared to the delta results
dataInputQuantCat <-
  dataInputQuant %>%  
  filter(dilution_assay == 100) %>% 
  dplyr::left_join(threshold, by = c("isotype" = "assay", "analyte" = "antigene")) %>% 
  mutate(serostatus = if_else(dataIn > median, 1, 0),
         serostatus_cat = case_when(dataIn > high ~ "positive",
                                    dataIn > median & dataIn <= high ~ "borderline positive",
                                    dataIn > low & dataIn <= median ~ "borderline negative",
                                    dataIn <= low ~ "negative"),
         serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"),
                                 ordered = TRUE)) 

dataInputQuantCatDelta <-
  dataInputQuantCat %>% 
  filter(dilution_assay == 100 & analyte == "Delta") %>% 
  dplyr::select(experiment, plate_assay, isotype, sampleID_meta, dataIn,
                serostatus, serostatus_cat)

dataInput <-
  dataInputQuantCat %>% 
  dplyr::left_join(dataInputQuantCatDelta, by = c("experiment", "plate_assay", "isotype", "sampleID_meta"), suffix = c("", ".delta")) %>% 
  unique()
#
# End determine serostatus
##



##
# 04) Determine infection group
dataInputLDAAll <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% 
  dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta,
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 

dataInputLDADelta <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% 
  dplyr::filter(serostatus.delta == 1) %>% 
  dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta, 
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta)

dataInputLDAHigh <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% 
  dplyr::filter(log10(dataIn.delta) >= 2.08) %>% 
  dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta, 
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 


dataInputLDAAllwide <-
  dataInputLDAAll %>% 
  dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

dataInputLDADeltawide <-
  dataInputLDADelta %>% 
  dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

dataInputLDAHighwide <-
  dataInputLDAHigh %>% 
  dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")


predictAll <- predict(lDAPanelsAll , dataInputLDAAllwide[4:ncol(dataInputLDAAllwide)])
predictDelta <- predict(lDAPanelsDelta , dataInputLDADeltawide[4:ncol(dataInputLDADeltawide)])
predictHigh <- predict(lDAPanelsHigh , dataInputLDAHighwide[4:ncol(dataInputLDAHighwide)])


dataInputGLMAllwide <-
  dataInputLDAAll %>% 
  dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn),
         analyte = str_replace(analyte, "-", ".")) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

#predictGLM <- predict(aucCombinedIgGPos, dataInputGLMAllwide)

predictionGLMDelta <- 
  data.frame(GLM_value =
               predict(aucCombinedIgGPos, dataInputGLMAllwide)) %>% 
  mutate(Combined_infected = (if_else(GLM_value > Combined_spec, "infected", "uninfected")))

ATI_infected <-
  dataInputGLMAllwide %>% 
  mutate(ATI_infected = (if_else(ATI.N > ATI_N_spec, "infected", "uninfected"))) %>% 
  dplyr::select(sampleID_meta, ATI_infected)#,
#ATI_infected = if_else(is.na(MPXV), ("Seronegative"),  ATI_infected))

##
# Generate prediction for different groups
# MPXV
predictAllpost <- predictAll$posterior[,4] 
predictDeltapost <- predictDelta$posterior[,4] 
predictHighpost <- predictHigh$posterior[,4] 

# Pre
predictAllpostPre <- predictAll$posterior[,1] 
predictDeltapostPre <- predictDelta$posterior[,1] 
predictHighpostPre <- predictHigh$posterior[,1] 

# MVA
predictAllpostMVA <- predictAll$posterior[,2] 
predictDeltapostMVA <- predictDelta$posterior[,2] 
predictHighpostMVA <- predictHigh$posterior[,2] 

# CPXV
predictAllpostCPXV <- predictAll$posterior[,3] 
predictDeltapostCPXV <- predictDelta$posterior[,3] 
predictHighpostCPXV <- predictHigh$posterior[,3] 

LDAAll <-
  dataInputLDAAllwide %>% 
  dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
  add_column(MPXV_all = predictAllpost,
             Pre_all = predictAllpostPre,
             MVA_all = predictAllpostMVA,
             CPXV_all = predictAllpostCPXV)

LDADelta <-
  dataInputLDADeltawide %>% 
  dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
  add_column(MPXV_delta = predictDeltapost,
             Pre_delta = predictDeltapostPre,
             MVA_delta = predictDeltapostMVA,
             CPXV_delta = predictDeltapostCPXV)

LDAHigh <-
  dataInputLDAHighwide %>% 
  dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
  add_column(MPXV_high = predictHighpost,
             Pre_high = predictHighpostPre,
             MVA_high = predictHighpostMVA,
             CPXV_high = predictHighpostCPXV)



unique(dataInput$analyte)

##
# Generate final data output

# Select 

dataOutputIgG <-
  dataInput %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  select(experiment, plate_assay, date, operator, 
         sampleID_meta, isotype, highBG, lowIg,
         analyte, data, dataIn, serostatus, serostatus_cat, serostatus.delta, serostatus_cat.delta) %>% 
  pivot_wider(names_from = analyte, values_from = c(data, dataIn, serostatus, serostatus_cat)) %>% 
  left_join(LDAAll, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
  left_join(LDADelta, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
  left_join(LDAHigh, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
  mutate(MPXV_all = round(MPXV_all, digits = 2),
         MPXV_delta = round(MPXV_delta, digits = 2),
         MPXV_high = round(MPXV_high, digits = 2),
         Pre_all = round(Pre_all, digits = 2),
         Pre_delta = round(Pre_delta, digits = 2),
         Pre_high = round(Pre_high, digits = 2),
         MVA_all = round(MVA_all, digits = 2),
         MVA_delta = round(MVA_delta, digits = 2),
         MVA_high = round(MVA_high, digits = 2),
         CPXV_all = round(CPXV_all, digits = 2),
         CPXV_delta = round(CPXV_delta, digits = 2),
         CPXV_high = round(CPXV_high, digits = 2),
         Final_all = case_when(MPXV_all > 0.95 ~ "MPXV",
                               Pre_all > 0.95 ~ "Pre",
                               MVA_all > 0.95 ~ "MVA",
                               CPXV_all > 0.95 ~ "CPXV",
                               TRUE ~ "Ambiguous"),
         Final_delta = case_when(MPXV_delta > 0.95 ~ "MPXV",
                                 Pre_delta > 0.95 ~ "Pre",
                                 MVA_delta > 0.95 ~ "MVA",
                                 CPXV_delta > 0.95 ~ "CPXV",
                                 is.na(MPXV_delta) ~ "NA",
                                 TRUE ~ "Ambiguous"),
         Final_high = case_when(MPXV_high > 0.95 ~ "MPXV",
                                Pre_high > 0.95 ~ "Pre",
                                MVA_high > 0.95 ~ "MVA",
                                CPXV_high > 0.95 ~ "CPXV",
                                is.na(MPXV_high) ~ "NA",
                                TRUE ~ "Ambiguous"),
         Final_combined = case_when(Final_high != "NA" 
                                    & Final_high == Final_delta ~ Final_high,
                                    Final_high != "NA" & Final_high != "Ambiguous"
                                    & Final_high != Final_delta ~ Final_high,
                                    Final_high != "NA" & Final_high == "Ambiguous"
                                    & Final_high != Final_delta ~ Final_delta,
                                    Final_high == "NA" & Final_delta != "NA" ~ Final_delta,
                                    Final_delta == "NA" ~ "seronegative"),
         Confidence = case_when(Final_high != "NA" 
                                & Final_high == Final_delta ~ "Highest",
                                Final_high != "NA" & Final_high != "Ambiguous"
                                & Final_high != Final_delta ~ "High",
                                Final_high != "NA" & Final_high == "Ambiguous"
                                & Final_high != Final_delta ~ "Medium",
                                Final_high == "NA" & Final_delta != "NA" ~ "Medium",
                                Final_delta == "NA" ~ "seronegative"),
         ratio_A29_A27L = dataIn_A29/dataIn_A27L,
         ratio_E8_D8L = dataIn_E8/dataIn_D8L,
         ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
         ratio_B6_B5R = dataIn_B6/dataIn_B5R) %>% 
  dplyr::left_join(ATI_infected, by = "sampleID_meta") %>%
  unique() %>%  
  add_column(predictionGLMDelta) %>% 
  mutate(#Combined_infected = if_else(is.na(MPXV_delta), ("seronegative"),  Combined_infected),
    ATI_infected = if_else(is.na(MPXV_delta), ("seronegative"),  ATI_infected),
    Final_combined = case_when(is.na(Final_combined) & MPXV_high > 0.75 ~ "MPXV", # New 2023-08-31
                               is.na(Final_combined) & MVA_high > 0.75 ~ "MVA",
                               TRUE ~ Final_combined))




export(dataOutputIgG, "output/Results_report_wide.xlsx")
export(dataInputQuantCat, "output/Input_quantified.xlsx")


##
# 06) Generate plots for graphical analysis

# Plot 1: Plot controls
#plot_controls <-
dataInputQuant %>% 
  select(dilution_assay, isotype, sampleID_meta, dilution_assay, lowIg, positiveIg, highBG, background) %>% 
  filter(dilution_assay == 100) %>% 
  unique() %>% 
  pivot_longer(cols = c(positiveIg, background), names_to = "Control", values_to = "MFI") %>% 
  pivot_longer(cols = c(lowIg, highBG), names_to = "QC_flags", values_to = "flags") %>%
  # filter(!is.na(flags)) %>% 
  ggplot(mapping = aes(x = sampleID_meta, y = MFI, color = flags)) +
  geom_point() +
  theme_bw()+
  scale_color_colorblind() +
  facet_grid(Control~isotype, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) 

# Plot 2: Plot standards
plot_standard_plate_1 <-
  plot(standard_plate_1) +
  theme_bw() +
  theme( strip.background = element_blank(),
         legend.position = "none")

plot_standard_plate_2 <-
  plot(standard_plate_2) +
  theme_bw() +
  theme( strip.background = element_blank(),
         legend.position = "none")

# Plot 3: Plot rawdata over quantified data
dataInputQuant %>% 
  select(dilution_assay, analyte, isotype, sampleID_meta, data, dataIn) %>% 
  filter(dilution_assay == 100) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = log10(data), y = log10(dataIn), color = isotype)) +
  geom_point() +
  theme_bw()+
  scale_color_colorblind() +
  facet_wrap("analyte") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) 


dataInputQuant %>% 
  left_join(dataOutput, by = c("sampleID_meta", "plate_assay", "experiment")) %>% 
  select(dilution_assay, analyte, isotype, sampleID_meta, data, dataIn, plate_assay, Final_combined) %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  # filter(plate_assay == 2) %>% 
  filter(!(analyte %in% c("VACV"))) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = analyte, y = log10(dataIn), color = sampleID_meta, group = sampleID_meta)) +
  geom_point() +
  geom_line() +
  theme_bw()+
  # scale_color_colorblind() +
  facet_wrap("Final_combined") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) 

plot_predictions <- 
  dataInputQuant %>% 
  left_join(dataOutput, by = c("sampleID_meta", "plate_assay", "experiment")) %>% 
  select(dilution_assay, analyte, isotype, sampleID_meta, data, dataIn, plate_assay, Final_combined) %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  # filter(plate_assay == 2) %>% 
  filter(!(analyte %in% c("VACV"))) %>% 
  unique() %>% 
  mutate(Final_combined = factor(Final_combined, levels = c("MVA",
                                                            "CPXV",
                                                            "MPXV",
                                                            "Pre",
                                                            "Ambiguous",
                                                            "seronegative"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Final_combined, y = log10(dataIn), fill = Final_combined)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  #  geom_line() +
  theme_bw()+
  scale_fill_manual(values = colorblind_pal()(8)[2:7]) +
  facet_wrap("analyte") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank()) 

ggsave("output/plot_standard_plate_1.png", plot_standard_plate_1,
       width = 6, height = 5, dpi = 600)

ggsave("output/plot_standard_plate_2.png", plot_standard_plate_2,
       width = 6, height = 5, dpi = 600)

ggsave("output/plot_controls.png", plot_controls,
       width = 14, height = 5, dpi = 600)

ggsave("output/plot_predictions.png", plot_predictions,
       width = 8, height = 7, dpi = 600)
