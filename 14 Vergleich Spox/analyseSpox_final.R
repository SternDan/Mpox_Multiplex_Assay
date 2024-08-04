##
# Analysis SPox study
# Assay validation
# Compare results from predictions with metadata 
# Daniel Stern RKI
# 04.08.2024
# V3: Ohne L1 und M1
# V4: More data and optimized algorithm
# V5: Final data
# Final: Final check for data sharing
##

# Prepare environment and load packages
rm(list = ls(all.names = TRUE))
library("rio")
library("tidyverse")
library("ggthemes")
library("ggpubr")
library("caret")
library("pROC")
library("freqtables")
library("epitools")

# Source functions
source("functions/determineCutoffFunction.R")
source("functions/rocFunction.R")

# Import data
inputDataFiles <- list.files("input_final/")

outlist <- list()
for(i in 1:length(inputDataFiles)){
  dataout <- import(paste("input_final/", inputDataFiles[i], sep = ""))
  outlist[[i]] <- dataout
}

dataInMeasured <- bind_rows(outlist)

dataInMeta <- import("MPOX_clinic2.xlsx") %>% 
  dplyr::rename(sampleID_meta = `LAB ID`)

dataIn <-
  dataInMeasured %>% 
  left_join(dataInMeta, by= c("sampleID_meta")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)),
         MPX_vac_final = case_when(MPX_VACC_STATUS_1 == "Yes" ~ "Yes",
                                   MPX_vac_jessen == "Yes" ~ "Yes",
                                   MPX_vac_driesener == "Yes" ~ "Yes", 
                                   MPX_VACC_STATUS_1 == "No" ~ "No",
                                   MPX_vac_jessen == "No" ~ "No",
                                   MPX_vac_driesener == "No" ~ "No"),
         MPX_diagnose_final = if_else(MPX_diagnose_final == "Unbekannt", "Yes", MPX_diagnose_final),
         LDA_prediction = if_else(is.na(Final_combined), "No Pred", "Prediction"),
         Final_combined =          case_when(!is.na(Final_high) 
                                             & Final_high != "Ambiguous" ~ Final_high,
                                             !is.na(Final_delta)
                                             & Final_delta != "Ambiguous" ~ Final_delta,
                                             Final_all != "Ambiguous" ~ Final_all,
                                             is.na(Final_delta)  ~ "seronegative")) %>% 
  filter(grepl("SPox", sampleID_meta))

save(dataIn, file = "data/dataInSpox.R")
export(dataIn, "output/2023-09-14_SPox_classified_IgG.xlsx")


# Import datainput before classification IgG and IgM
inputDataFilesQuant <- list.files("input_quant/")

outlistquant <- list()
for(i in 1:length(inputDataFilesQuant)){
  dataout <- import(paste("input_quant/", inputDataFilesQuant[i], sep = ""))
  outlistquant[[i]] <- dataout
}

dataInQuant <- bind_rows(outlistquant)

dataInQuant <-
  dataInQuant %>% 
  left_join(dataInMeta, by= c("sampleID_meta")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)),
         MPX_vac_final = case_when(MPX_VACC_STATUS_1 == "Yes" ~ "Yes",
                                   MPX_vac_jessen == "Yes" ~ "Yes",
                                   MPX_vac_driesener == "Yes" ~ "Yes", 
                                   MPX_VACC_STATUS_1 == "No" ~ "No",
                                   MPX_vac_jessen == "No" ~ "No",
                                   MPX_vac_driesener == "No" ~ "No"),
         MPX_diagnose_final = if_else(MPX_diagnose_final == "Unbekannt", "Yes", MPX_diagnose_final)) %>% 
  filter(grepl("SPox", sampleID_meta))

save(dataInQuant, file = "data/dataInQuantSpox.R")
export(dataInQuant, "output/2023-09-08_SPox_quantified_IgG_IgM.xlsx")


##
# Load NK-Data and combine with metadata
inputDataFilesNC <- list.files("input_nc/")

outlistNC <- list()
for(i in 1:length(inputDataFilesNC)){
  dataout <- import(paste("input_nc/", inputDataFilesNC[i], sep = ""))
  outlistNC[[i]] <- dataout
}

dataInNC <- bind_rows(outlistNC) %>% 
  filter(!grepl("Dummy", sampleID_meta)) %>% 
  filter(sampleID_meta != "DS-2") %>% 
  mutate(sampleID_meta_merge = grep("[0-9][0-9]", sampleID_meta))

dataInMetaNC <- import("MMR Kontrollseren MPox Panel.xlsx") %>% 
  mutate(sampleID_meta = str_remove(Pseudonym, "NRZ-MMR-ctrlMpox-"),
         sampleID_meta_merge = grep("[0-9][0-9]", sampleID_meta),
         age_group = `Altersgruppe (YOB)`,
         childhood_immu = if_else(age_group %in% c("1935-1944",
                                                   "1945-1954",
                                                   "1955-1964",
                                                   "1965-1974"), 1, 0),
         measles_IgG = case_when(`Masern IgG` == "+" ~ 1,
                                 `Masern IgG` == "-" ~ 0,
                                 `Masern IgG` == "?" ~ 99),
         measles_IgM = case_when(`Masern IgM` == "+" ~ 1,
                                 `Masern IgM` == "-" ~ 0,
                                 `Masern IgM` == "?" ~ 99),
         measles_PCR = case_when(`Masern PCR` == "+" ~ 1,
                                 `Masern PCR` == "-" ~ 0,
                                 `Masern PCR` == "?" ~ 99),
         mumps_IgG = case_when(`Mumps IgG` == "+" ~ 1,
                               `Mumps IgG` == "-" ~ 0,
                               `Mumps IgG` == "?" ~ 99),
         mumps_IgM = case_when(`Mumps IgM` == "+" ~ 1,
                               `Mumps IgM` == "-" ~ 0,
                               `Mumps IgM` == "?" ~ 99),
         mumps_PCR = case_when(`Mumps PCR` == "+" ~ 1,
                               `Mumps PCR` == "-" ~ 0,
                               `Mumps PCR` == "?" ~ 99),
         rubella_IgG = case_when(`Röteln IgG` == "+" ~ 1,
                               `Röteln IgG` == "-" ~ 0,
                               `Röteln IgG` == "?" ~ 99),
         rubella_IgM = case_when(`Röteln IgM` == "+" ~ 1,
                               `Röteln IgM` == "-" ~ 0,
                               `Röteln IgM` == "?" ~ 99),
         rubella_PCR = case_when(`Röteln PCR` == "+" ~ 1,
                               `Röteln PCR` == "-" ~ 0,
                               `Röteln PCR` == "?" ~ 99)) %>% 
  select(sampleID_meta_merge:rubella_PCR)

dataInNCCombined <-
  dataInNC %>% 
  left_join(dataInMetaNC, by= c("sampleID_meta_merge")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)),
         LDA_prediction = if_else(is.na(Final_combined), "No Pred", "Prediction"),
         Final_combined =          case_when(!is.na(Final_high) 
                                             & Final_high != "Ambiguous" ~ Final_high,
                                             !is.na(Final_delta)
                                             & Final_delta != "Ambiguous" ~ Final_delta,
                                             Final_all != "Ambiguous" ~ Final_all,
                                             is.na(Final_delta)  ~ "seronegative"))

save(dataInNCCombined, file = "data/dataInNCCombined.R")


# Import datainput NC before classification IgG and IgM
inputDataFilesQuantNC <- list.files("input_nc_quant/")

outlistquantNC <- list()
for(i in 1:length(inputDataFilesQuantNC)){
  dataout <- import(paste("input_nc_quant/", inputDataFilesQuantNC[i], sep = ""))
  outlistquantNC[[i]] <- dataout
}

dataInQuantNC <- bind_rows(outlistquantNC) %>% 
  filter(!grepl("Dummy", sampleID_meta)) %>% 
  filter(sampleID_meta != "DS-2") %>% 
  mutate(sampleID_meta_merge = as.numeric(str_remove(sampleID_meta, "NK-")))
  

dataInQuantNCCombined <-
  dataInQuantNC %>% 
  left_join(dataInMetaNC, by= c("sampleID_meta_merge")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)))

save(dataInQuantNCCombined, file = "data/dataInQuantNCCombined.R")


# Import data of repetition measurements
inputDataFilesRep <- list.files("input_rep/")

outlistRep <- list()
for(i in 1:length(inputDataFilesRep)){
  dataout <- import(paste("input_rep/", inputDataFilesRep[i], sep = ""))
  outlistRep[[i]] <- dataout
}

dataInMeasuredRep <- bind_rows(outlistRep)


dataInRep <-
  dataInMeasuredRep %>% 
  left_join(dataInMeta, by= c("sampleID_meta")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)),
         MPX_vac_final = case_when(MPX_VACC_STATUS_1 == "Yes" ~ "Yes",
                                   MPX_vac_jessen == "Yes" ~ "Yes",
                                   MPX_vac_driesener == "Yes" ~ "Yes", 
                                   MPX_VACC_STATUS_1 == "No" ~ "No",
                                   MPX_vac_jessen == "No" ~ "No",
                                   MPX_vac_driesener == "No" ~ "No"),
         MPX_diagnose_final = if_else(MPX_diagnose_final == "Unbekannt", "Yes", MPX_diagnose_final),
         LDA_prediction = if_else(is.na(Final_combined), "No Pred", "Prediction"),
         Final_combined =          case_when(!is.na(Final_high) 
                                             & Final_high != "Ambiguous" ~ Final_high,
                                             !is.na(Final_delta)
                                             & Final_delta != "Ambiguous" ~ Final_delta,
                                             Final_all != "Ambiguous" ~ Final_all,
                                             is.na(Final_delta)  ~ "seronegative")) %>% 
  filter(grepl("SPox", sampleID_meta))

save(dataInRep, file = "data/dataInSpoxRep.R")


# Import datainput before classification IgG and IgM
inputDataFilesQuantRep <- list.files("input_rep_quant/")

outlistquantRep <- list()
for(i in 1:length(inputDataFilesQuantRep)){
  dataout <- import(paste("input_rep_quant/", inputDataFilesQuantRep[i], sep = ""))
  outlistquantRep[[i]] <- dataout
}

dataInQuantRep <- bind_rows(outlistquantRep)


dataInQuantRep <-
  dataInQuantRep %>% 
  left_join(dataInMeta, by= c("sampleID_meta")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)),
         MPX_vac_final = case_when(MPX_VACC_STATUS_1 == "Yes" ~ "Yes",
                                   MPX_vac_jessen == "Yes" ~ "Yes",
                                   MPX_vac_driesener == "Yes" ~ "Yes", 
                                   MPX_VACC_STATUS_1 == "No" ~ "No",
                                   MPX_vac_jessen == "No" ~ "No",
                                   MPX_vac_driesener == "No" ~ "No"),
         MPX_diagnose_final = if_else(MPX_diagnose_final == "Unbekannt", "Yes", MPX_diagnose_final)) %>% 
  filter(grepl("SPox", sampleID_meta))

save(dataInQuantRep, file = "data/dataInQuantSpoxRep.R")
