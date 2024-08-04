##
# Analyse Sera with a Mpox-specific Multiplex Assay
# Version Alpha
# 2023-06-02
# Daniel Stern
# RKI
##

##
# The following scripts analyses results from Bio-Plex 200 measurements
# which have been performed, to determine the serostatus against 
# Mpox antigens. 
# The assay returns the serostatus against IgG and IgM (Positive, Negative),
# quantitative results expressed as bining antibody units (BAU) by 
# quantification compared against VIG, as well as three estimates for the
# likely source of OPXV-specific antibodies: MVA-Immunisation, CPXV-Infection,
# MPXV-Infection with different levels of confidence ("low", "medium", "high"). 
# Additionally, the levels of neutralizing antibodies potentially present are 
# predicted from the mutliplexed response.
#
# For each measurement, two input files are needed
# 1) Excelfile with results exported from the Bio-Plex 200
# 2) Excelfile with metadata describing the input file
#
# The following analysis steps are performed
# 1) Load and merge input files
# 2) Quantify results 
# 3) Determine serostatus for IgG and IgM
# 4) Predict source of infection
# 5) Predict level of neutralizing antibodies
# 6) Generate report and data output
##

##
# Clean environment
rm(list = ls(all.names = TRUE))

##
# Load necessary libraries
library(rio)
library(tidyverse)
library(lubridate)
library(drLumi)
library(MASS)
library(caret)
library(glmnet)

##
# Source necessary functions and data from previous analysis
source("functions/function_read_multiplex.R", encoding = "UTF-8") ## Function for data import
source("functions/functions_quantify.R", encoding = "UTF-8") ## Function for data quantification
ecdata <- import("data/EC_input.xlsx") %>%  # Load ec-file for quantification
  pivot_longer(c(-sample), names_to = "analyte", values_to = "ec")
load("data/threshold.Rdata") # Load cutoff values for IgG and IgM 
load("data/LDAmodels.Rdata") # Load LDA models 
load("data/NTLassoReg.Rdata")
load("../12 ATI-N cutoff/output/GLMmodels.RData")
load("../12 ATI-N cutoff/output/threshold_MaxSpec.RData")

ATI_N_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "ATI.N" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()

Combined_spec <- threshold_MaxSpec %>% 
  filter(panel == "Pos" & antigen == "Combined" & parameter == "Threshold") %>% 
  dplyr::select(values) %>% 
  pull()



##
# 01) Load Multiplex data and metadata ----
inputDataFiles <- list.files("input/")[grepl("Pox_", list.files("input/")) &
                                         !grepl("Meta", list.files("input/"))]
inputMetaFiles <- list.files("input/")[grepl("Pox_", list.files("input/")) &
                                         grepl("Meta", list.files("input/"))]
outlist <- list()
for(i in 1:length(inputDataFiles)){
  dataout <- read_multiplex(paste("input/", inputDataFiles[i], sep = ""),
                            paste("input/", inputMetaFiles[i], sep = ""))
  outlist[[i]] <- dataout
}

input <- bind_rows(outlist) %>% 
  filter(!is.na(data)) ## Removes NAs introduced by Description column in data file

##
# Remove files from import
rm(outlist, i, dataout, inputDataFiles, inputMetaFiles, read_multiplex)
# End import data
##



##
# 02) Quantify results ----
# Select data containing the results for "Standard1"
dataInputQC <-
  input %>%
  filter(assaytype == "Multiplex" & analyte == "HSA", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, isotype, dilution) %>% 
  summarise(highBG = if_else(data > 500, T, F))

dataInputQCIgG <-
  input %>%
  filter(assaytype == "Multiplex" & analyte == "ahIgG", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, isotype, dilution) %>% 
  summarise(lowIgG = if_else(data < 10000, T, F))

# Normalize data
dataInputNorm <-
  input %>%
  left_join(dataInputQC, by = c("experiment", "plate","sampleID_assay", "isotype",
                                "dilution")) %>%
  left_join(dataInputQCIgG, by = c("experiment", "plate","sampleID_assay", "isotype",
                                   "dilution")) %>%
  mutate(exp_plate = paste(experiment, plate, sep = ".")) %>% 
  unique()

rm(input, dataInputQC, dataInputQCIgG)

# Fit IgG data
multiInputList <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgG") %>% 
  group_split(experiment, plate_assay)
names(multiInputList) <- unique(dataInputNorm$exp_plate[dataInputNorm$assaytype == "Multiplex"])

outlist <- list()
names_list <- list()
for(i in 1:length(multiInputList)){
  outlist[[i]] <- quantifyMultiplex(ecdata = ecdata, multiInputList = multiInputList[[i]])
  names_list[[i]] <-names(multiInputList)[i]
  names_list[[i]]$error <-if_else(c("SC_NON_CONVERGENCE") %in% unique(outlist[[i]]$dataOutput$warning), names(multiInputList)[i], NA_character_)
  names_list[[i]]$analytes <- unique(outlist[[i]]$dataOutput$analyte[outlist[[i]]$dataOutput$warning == "SC_NON_CONVERGENCE"])
  names_list[[i]]$batch <- unique(outlist[[i]]$dataOutput$batch[outlist[[i]]$dataOutput$warning == "SC_NON_CONVERGENCE"])
}


outlist[[1]]$standards$A27L
names_list[[1]]$batch

names(outlist) <- names(multiInputList)
names(names_list) <- names(multiInputList)

####
# Generate algorithm to replace standard curves with not fit with closest
# standard curve, that did not generate an error

##
# Read data from standard curve with error
names <- as.data.frame(do.call(rbind, as.vector(lapply(names_list, `[[`, 1))))
error <- as.data.frame(do.call(rbind, lapply(names_list, `[[`, 2))) %>% 
  filter(!is.na(V1))

outlist[error %in% names]


listdata <- lapply(outlist, `[[`, 2)

##
# DEV: Identify results with errors in fit


c("SC_NON_CONVERGENCE") %in% unique(outlist[[13]]$dataOutput$warning)

analytes <- unique(outlist[[25]]$dataOutput$analyte[outlist[[25]]$dataOutput$warning == "SC_NON_CONVERGENCE"])
batch <- unique(outlist[[25]]$dataOutput$batch[outlist[[25]]$dataOutput$warning == "SC_NON_CONVERGENCE"])


liststandards <- lapply(outlist, `[[`, 1)
listdata <- lapply(outlist, `[[`, 2)

# Fit IgM data using IgG standard curves
multiInputListIgM <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgM") %>% 
  group_split(experiment, plate_assay)

outlistIgM <- list()
for(i in 1:length(multiInputListIgM)){
  outlistIgM[[i]] <- quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgM[[i]],
                                            allanalytesIn = outlist[[i]]$standards, i)
}
listdataIgM <- lapply(outlistIgM, `[[`, 2)


igginput <- do.call(rbind,listdata)

dataInputQuant <-
  rbind(do.call(rbind,listdata),
        do.call(rbind,listdataIgM)) %>% 
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
  group_by(analyte, assaytype, experiment, plate_assay) %>% 
  mutate(quant_conc_fix = case_when(is.na(quant_conc) & grepl("above", remarks) ~ max(quant_conc, na.rm = TRUE),
                                    is.na(quant_conc) & grepl("below", remarks) ~ min(quant_conc, na.rm = TRUE),
                                    TRUE ~ quant_conc),
         quant_conc_fix = if_else(quant_conc_fix < 1, 1, quant_conc_fix)) %>% 
  ungroup()

#
# End quantify results
## 



##
# 03) Determine serostatus for IgG and IgM compared to the delta results
dataInputQuantCat <-
  dataInputQuant %>%  
  rename(dataIn = quant_conc_fix) %>% 
  left_join(threshold, by = c("isotype" = "assay", "analyte" = "antigene", "assaytype")) %>% 
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
  filter(dilution_assay == 100 & analyte == "Delta" & assaytype == "Multiplex") %>% 
  dplyr::select(experiment, plate_assay, isotype, sample, dataIn,
                serostatus, serostatus_cat)

dataInput <-
  dataInputQuantCat %>% 
  left_join(dataInputQuantCatDelta, by = c("experiment", "plate_assay", "isotype", "sample"), suffix = c("", ".delta"))
#
# End determine serostatus
##



##
# 04) Determine infection group
dataInputLDAAll <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::select(experiment, plate_assay, isotype, sample, plate_assay,
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 

dataInputLDADelta <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::filter(serostatus.delta == 1) %>% 
  dplyr::select(experiment, plate_assay, isotype, sample, plate_assay,
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta)

dataInputLDAHigh <- dataInput %>% 
  dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  dplyr::filter(log10(dataIn.delta) >= 2.08) %>% 
  dplyr::select(experiment, plate_assay, isotype, sample, plate_assay,
                analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 


dataInputLDAAllwide <-
  dataInputLDAAll %>% 
  dplyr::select(experiment, sample, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

dataInputLDADeltawide <-
  dataInputLDADelta %>% 
  dplyr::select(experiment, sample, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

dataInputLDAHighwide <-
  dataInputLDAHigh %>% 
  dplyr::select(experiment, sample, analyte, dataIn) %>% 
  mutate(dataIn = log10(dataIn)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")


predictAll <- predict(lDAPanelsAll , dataInputLDAAllwide[3:ncol(dataInputLDAAllwide)])
predictDelta <- predict(lDAPanelsDelta , dataInputLDADeltawide[3:ncol(dataInputLDADeltawide)])
predictHigh <- predict(lDAPanelsHigh , dataInputLDAHighwide[3:ncol(dataInputLDAHighwide)])

predictAllpost <- predictAll$posterior[,4] 
predictDeltapost <- predictDelta$posterior[,4] 
predictHighpost <- predictHigh$posterior[,4] 

LDAAll <-
  dataInputLDAAllwide %>% 
  dplyr::select(experiment, sample) %>% 
  add_column(MPXV_all = predictAllpost)

LDADelta <-
  dataInputLDADeltawide %>% 
  dplyr::select(experiment, sample) %>% 
  add_column(MPXV_delta = predictDeltapost)

LDAHigh <-
  dataInputLDAHighwide %>% 
  dplyr::select(experiment, sample) %>% 
  add_column(MPXV_high = predictHighpost)

##
# Predict level of neutralizing antibodies

predNT <-
  predict(best_model, s = best_lambda, newx = as.matrix(dataInputLDAAllwide[3:ncol(dataInputLDAAllwide)])) %>% 
  as.data.frame() %>% 
  mutate(NT_titer = 10^s1)

NTAll <-
  dataInputLDAAllwide %>% 
  dplyr::select(experiment, sample) %>% 
  add_column(predNT)

##
# Generate final data output
##
# Generate final data output
dataOutput <-
  dataInputQuantCatDelta %>%
  pivot_wider(names_from = isotype, values_from = c(dataIn, serostatus, serostatus_cat)) %>%
  left_join(LDAAll, by = c("experiment", "sampleID_meta", "plate_assay")) %>%
  left_join(LDADelta, by = c("experiment", "sampleID_meta", "plate_assay")) %>%
  left_join(LDAHigh, by = c("experiment", "sampleID_meta", "plate_assay")) %>%
  left_join(NTAll, by = c("experiment", "sampleID_meta", "plate_assay")) %>%
  mutate(dataIn_IgG = round(dataIn_IgG, digits = 0),
         dataIn_IgM = round(dataIn_IgM, digits = 0),
         #    MPXV_all = round(MPXV_all, digits = 2),
         MPXV_delta = round(MPXV_delta, digits = 2),
         MPXV_high = round(MPXV_high, digits = 2),
         #     Pre_all = round(Pre_all, digits = 2),
         Pre_delta = round(Pre_delta, digits = 2),
         Pre_high = round(Pre_high, digits = 2),
         #     MVA_all = round(MVA_all, digits = 2),
         MVA_delta = round(MVA_delta, digits = 2),
         MVA_high = round(MVA_high, digits = 2),
         #     CPXV_all = round(CPXV_all, digits = 2),
         CPXV_delta = round(CPXV_delta, digits = 2),
         CPXV_high = round(CPXV_high, digits = 2),
         NT_titer = round(NT_titer, digits = 0),
         #   Final_all = case_when(MPXV_all > 0.75 ~ "MPXV",
         #                         Pre_all > 0.75 ~ "Pre",
         #                         MVA_all > 0.75 ~ "MVA",
         #                         CPXV_all > 0.75 ~ "CPXV",
         #                         TRUE ~ "Ambiguous"),
         Final_delta = case_when(MPXV_delta > 0.75 ~ "MPXV",
                                 Pre_delta > 0.75 ~ "Pre",
                                 MVA_delta > 0.75 ~ "MVA",
                                 CPXV_delta > 0.75 ~ "CPXV",
                                 is.na(MPXV_delta) ~ "NA",
                                 TRUE ~ "Ambiguous"),
         Final_high = case_when(MPXV_high > 0.75 ~ "MPXV",
                                Pre_high > 0.75 ~ "Pre",
                                MVA_high > 0.75 ~ "MVA",
                                CPXV_high > 0.75 ~ "CPXV",
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
                                Final_delta == "NA" ~ "seronegative"))
export(dataOutput, "output/Results_report.xlsx")
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

