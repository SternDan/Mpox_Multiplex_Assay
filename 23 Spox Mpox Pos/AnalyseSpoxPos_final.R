##
# Generate Age-group harmonized input data for 
# Analysis in other R-scripts
# Daniel Stern
# 2024/08/04 
# Final version -> only generation of output needed by other scripts
# Previous script was used to further analyse spox data
# see archive for code
##


# Load libraries
library(rio)
library(tidyverse)
library(ggpubr)
library(stringr)
library(ggthemes)

# Clean environment
rm(list = ls(all.names = TRUE))

# Load data
# Load Metadata for Spox containing PIN and SampleID
dataSPox <- import("input/SPox_Proben.xlsx") %>% 
  rename(sampleID_metadata = `Lab ID`,
         dateSampleIn = `Eingang am`,
         patientID = `Patienten ID (Barcode)`,
         sender = Einsender,
         comments = Kommentar)
# Load prediction of ML-Algorithms
dataPred <- import("input/predictionSpoxComplete.xlsx")

# Load input for ML -> filter only on Spox panel
dataInputSpox <- import("input/dataInputAll.csv") %>% 
  filter(panel %in% c("SPox")) %>% 
  select(sampleID_metadata, isotype, analyte, c(data:MPX_vac_final), panel_detail) %>% 
  unique()
dataAgeGroupsSpox <- import("input/MPOX_clinic_serology_MAIN.xlsx") %>% 
  select(PIN, sampleID_metadata = `LAB ID`, AGE, MPX_VACC_STATUS_1, 
         MPX_VACC_STATUS_2, POX_VACC_STATUS_1) %>% 
  mutate(Age_group_all = factor(AGE, levels = c("18-29",
                                                "30-39",
                                                "40-49",
                                                "50-59",
                                                "60+"), 
                                labels = c("< 30", 
                                           "< 40", 
                                           "< 50", 
                                           "< 60", 
                                           "60+"),
                                ordered = TRUE)) %>% 
  filter(!is.na(sampleID_metadata))

## Load all data for ML
dataInputAll <- import("input/dataInputAll.csv") 

## Filter MVA data from acute panel
dataInputAcuteMVA <-
  dataInputAll %>% 
  filter(panel == "MVA") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "L1R", "M1", 
                                              "D8L", "E8",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), 
                          labels = c("A27L", "A29L",
                                     "L1R", "M1R", 
                                     "D8L", "E8L",
                                     "H3L", "A33R", "A35R",
                                     "B5R", "B6R", "A5L", "ATI-C", 
                                     "ATI-N", "Delta", "VACV"), ordered = TRUE),
         Age_group = case_when(year_birth >= 1994 ~ "< 30",
                               year_birth >= 1984 ~ "< 40",
                               year_birth >= 1974 ~ "< 50",
                               year_birth >= 1964 ~ "< 60"),
         Age_group = factor(Age_group, levels = c("< 30", 
                                                  "< 40", 
                                                  "< 50", 
                                                  "< 60", 
                                                  "60+"),
                            ordered = TRUE))

## Filter MPXV data from acute panel
dataInputAcuteMPXV <-
  dataInputAll %>% 
  filter(panel == "MPXV") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "L1R", "M1", 
                                              "D8L", "E8",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), 
                          labels = c("A27L", "A29L",
                                     "L1R", "M1R", 
                                     "D8L", "E8L",
                                     "H3L", "A33R", "A35R",
                                     "B5R", "B6R", "A5L", "ATI-C", 
                                     "ATI-N", "Delta", "VACV"), ordered = TRUE),
         Age_group = case_when(year_birth >= 1994 ~ "< 30",
                               year_birth >= 1984 ~ "< 40",
                               year_birth >= 1974 ~ "< 50",
                               year_birth >= 1964 ~ "< 60",
                               year_birth < 1964 ~ "60+"),
         Age_group = factor(Age_group, levels = c("< 30", 
                                                  "< 40", 
                                                  "< 50", 
                                                  "< 60", 
                                                  "60+"),
                            ordered = TRUE))


## Combine data for Spox panel
dataCombined <- 
  dataInputSpox %>% 
  left_join(dataPred, by = c("sampleID_metadata")) %>% 
  left_join(dataSPox, by = c("sampleID_metadata")) %>% 
  left_join(dataAgeGroupsSpox, by = c("sampleID_metadata")) %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "L1R", "M1", 
                                              "D8L", "E8",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), 
                          labels = c("A27L", "A29L",
                                     "L1R", "M1R", 
                                     "D8L", "E8L",
                                     "H3L", "A33R", "A35R",
                                     "B5R", "B6R", "A5L", "ATI-C", 
                                     "ATI-N", "Delta", "VACV"), ordered = TRUE))

##
# Load NK data and join with age information
dataInputNK <- import("input/dataInputAll.csv") %>% 
  filter(panel == "Pre_New") %>% 
  mutate(sample = as.numeric(str_remove(sampleID_metadata, "NK-"))) %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "L1R", "M1", 
                                              "D8L", "E8",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), 
                          labels = c("A27L", "A29L",
                                     "L1R", "M1R", 
                                     "D8L", "E8L",
                                     "H3L", "A33R", "A35R",
                                     "B5R", "B6R", "A5L", "ATI-C", 
                                     "ATI-N", "Delta", "VACV"), ordered = TRUE))

dataNKMeta <- import("input/MMR Kontrollseren MPox Panel.xlsx") %>% 
  mutate(sample = as.numeric(str_extract(Pseudonym, "\\d+")),
         Age_cohort = factor(`Altersgruppe (YOB)`, levels = c("1935-1944",
                                                              "1945-1954",
                                                              "1955-1964",
                                                              "1965-1974",
                                                              "1975-1984",
                                                              "1985-1994",
                                                              "1995-2004",
                                                              "2005-2014",
                                                              "2015-2020"),
                             ordered = TRUE),
         Age_group = case_when(Age_cohort == "1935-1944" ~ "60+",
                               Age_cohort == "1945-1954" ~ "60+",
                               Age_cohort == "1955-1964" ~ "60+",
                               Age_cohort == "1965-1974" ~ "< 60",
                               Age_cohort == "1975-1984" ~ "< 50",
                               Age_cohort == "1985-1994" ~ "< 40",
                               Age_cohort == "1995-2004" ~ "< 30",
                               Age_cohort == "2005-2014" ~ "< 30",
                               Age_cohort == "2015-2020" ~ "< 30",
         ),
         Age_group = factor(Age_group, levels = c("< 30", 
                                                  "< 40", 
                                                  "< 50", 
                                                  "< 60", 
                                                  "60+"),
                            ordered = TRUE))



dataInputNKcombined <-
  dataInputNK %>% 
  left_join(dataNKMeta, by = "sample")


###
# Combine different panels into one dataframe for export and subsequent 
# analysis
# 1) dataCombined: Spox Panel
dataCombined_rbind <-
  dataCombined %>% 
  mutate(panel = "SPox") %>% 
  select(panel, panel_detail, isotype, analyte, sampleID_metadata, PIN,
         data, dataIn, serostatus, serostatus_cat, Age_group = Age_group_all,
         POX_VACC_STATUS_ML = POX_VACC_STATUS_1.x,
         POX_VACC_STATUS_SELF = POX_VACC_STATUS_1.y,
         MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
         MPX_vac_final,
         MPX_SYMPTOME, MPX_diagnose_final,
         real, starts_with("pred")) %>% 
  unique()

# 2) dataInputNKcombined: NK panel
dataInputNKcombined_rbind <-
  dataInputNKcombined %>% 
  mutate(PIN = NA,
         POX_VACC_STATUS_ML = NA,
         POX_VACC_STATUS_SELF = NA,
         MPX_VACC_STATUS_1 = NA,
         MPX_VACC_STATUS_2 = NA,
         MPX_vac_final = NA,
         MPX_SYMPTOME = NA,
         MPX_diagnose_final = NA,
         real = NA,
         pred_AllIgGIgM = NA,
         pred_AllIgG = NA,
         pred_EpiIgGIgM = NA,
         pred_EpiIgG = NA,
         pred_RFAllIgGIgM = NA,
         predxgboostAllIgGIgM = NA,
         pred_LDAAllIgGIgM = NA) %>% 
  select(panel, panel_detail, isotype, analyte, sampleID_metadata, PIN,
         data, dataIn, serostatus, serostatus_cat, Age_group,
         POX_VACC_STATUS_ML,
         POX_VACC_STATUS_SELF,
         MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
         MPX_vac_final,
         MPX_SYMPTOME, MPX_diagnose_final,
         real, starts_with("pred")) %>% 
  unique()

# 3: dataInputAcuteMVA
dataInputAcuteMVA_rbind <-
  dataInputAcuteMVA %>% 
  mutate(PIN = NA,
         POX_VACC_STATUS_ML = NA,
         POX_VACC_STATUS_SELF = NA,
         MPX_VACC_STATUS_1 = NA,
         MPX_VACC_STATUS_2 = NA,
         MPX_vac_final = NA,
         MPX_SYMPTOME = NA,
         MPX_diagnose_final = NA,
         real = NA,
         pred_AllIgGIgM = NA,
         pred_AllIgG = NA,
         pred_EpiIgGIgM = NA,
         pred_EpiIgG = NA,
         pred_RFAllIgGIgM = NA,
         predxgboostAllIgGIgM = NA,
         pred_LDAAllIgGIgM = NA) %>% 
  select(panel, panel_detail, isotype, analyte, sampleID_metadata, PIN,
         data, dataIn, serostatus, serostatus_cat, Age_group,
         POX_VACC_STATUS_ML,
         POX_VACC_STATUS_SELF,
         MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
         MPX_vac_final,
         MPX_SYMPTOME, MPX_diagnose_final,
         real, starts_with("pred")) %>% 
  unique()

# 4: dataInputActueMPXV
dataInputAcuteMPXV_rbind <-
  dataInputAcuteMPXV %>% 
  mutate(PIN = NA,
       #  POX_VACC_STATUS_ML = NA,
         POX_VACC_STATUS_SELF = NA,
         MPX_VACC_STATUS_1 = NA,
         MPX_VACC_STATUS_2 = NA,
         MPX_vac_final = NA,
         MPX_SYMPTOME = NA,
        # MPX_diagnose_final = NA,
         real = NA,
         pred_AllIgGIgM = NA,
         pred_AllIgG = NA,
         pred_EpiIgGIgM = NA,
         pred_EpiIgG = NA,
         pred_RFAllIgGIgM = NA,
         predxgboostAllIgGIgM = NA,
         pred_LDAAllIgGIgM = NA) %>% 
  select(panel, panel_detail, isotype, analyte, sampleID_metadata, PIN,
         data, dataIn, serostatus, serostatus_cat, Age_group,
         POX_VACC_STATUS_ML = POX_VACC_STATUS_1,
         POX_VACC_STATUS_SELF,
         MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
         MPX_vac_final,
         MPX_SYMPTOME, MPX_diagnose_final,
         real, starts_with("pred")) %>% 
  unique()


##
# Combine panels for export
dataClustering <- 
  dataCombined_rbind %>% 
  rbind(dataInputNKcombined_rbind, dataInputAcuteMVA_rbind,
        dataInputAcuteMPXV_rbind)

save(dataClustering, file = "output/dataClustering.Rdata")

