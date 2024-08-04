## 
# Generate Dataframe containing both MPox data and SPox data as datainput for ML
# Daniel Stern
# RKI
# 2024-08-04
# Version Final
##


# Load libraries and clean environment
library("rio")
library("tidyverse")
library("lubridate")
library("ggpubr")
library("ggthemes")

rm(list = ls(all.names = TRUE))

# Load data
load("../14 Vergleich Spox/data/dataInQuantSpox.R") # dataInQuant -> Spox Dataframe IgG
load("../14 Vergleich Spox/data/dataInSpox.R") # dataIn -> Spox Dataframe IgG for analysis
load("../14 Vergleich Spox/data/dataInQuantNCCombined.R") # dataInQuantNCCombined -> New NC panel
load("../14 Vergleich Spox/data/dataInNCCombined.R") # dataInNCCombined -> IgG Analysed and classified panel
load("../14 Vergleich Spox/data/dataInQuantSpoxRep.R") # dataInQuantRep -> Repeated measurement of select samples
load("../14 Vergleich Spox/data/dataInSpoxRep.R") # dataInRep -> IgG Analyse repeated measurments
load("../4 Method Comparison/output/dataInputQuantCat.Rdata")
load("../4 Method Comparison/output/dataInputQuantCatNK.Rdata")
load("../4 Method Comparison/output/threshold.Rdata")
load("../2 Import metadata/output/metadata_MVA_time.Rdata")
load("../2 Import metadata/output/metadata_MPXV_patients.Rdata")
cutoffDelta <- log10(threshold$median[threshold$assay == "IgG" & 
                                        threshold$antigene == "Delta" & 
                                        threshold$assaytype == "Multiplex"])

##
# Import data that have been matched with metadata
load("../7 Patient panel merge/output/dataInputMPXVmeta_all.Rdata") # dataInputMPXVmeta -> New: including all data
load("../7 Patient panel merge/output/dataInputMVAmeta.Rdata") # dataInputMVAmeta



metaMPXV <- dataInputMPXVmeta %>% 
  dplyr::select(sampleID_metadata, date_birthday, vaccination_pox) %>% 
  unique() %>% 
  mutate(year_birth = year(date_birthday),
         childhoodImmu = if_else(year_birth < 1975, 1, 0), 
         panel = "MPXV",
         serostatus_factor = as.factor("Post infection")) %>% 
  dplyr::select(-date_birthday)


metaMVA <- dataInputMVAmeta %>% 
  dplyr::select(sampleID_metadata, year_birth, childhoodImmu, serostatus_factor) %>% 
  unique() %>% 
  mutate(vaccination_pox = NA_character_,
         panel = "MVA")

metaCombined <- rbind(metaMPXV, metaMVA) %>% 
  mutate(panel_pre = if_else(serostatus_factor == "Pre", "MVA_Pre", panel))

preSera <- metaCombined$sampleID_metadata[metaCombined$serostatus_factor == "Pre"]

##
# Prepare dataframe
# 1) Combine panel and NK measurement
# 2) Mutate MVA Sera: Pre to 
dataInput <-
  rbind(dataInputQuantCat, dataInputQuantCatNK) %>% 
  filter(!is.na(panel_corr)) %>% 
  mutate(panel_in = case_when(panel_corr == "NK" ~ "Pre",
                              sampleID_metadata %in% preSera ~ "Pre",
                              TRUE ~ panel_corr),
         panel_inf = if_else(panel_in %in% c("CPXV", "MPXV"), "Infected", panel_in))


##
# Select only patient samples and select necessary columns
dataInputRaw <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(grepl("P-", sampleID_metadata) | grepl("S-", sampleID_metadata) |
                  grepl("MMR-", sampleID_metadata)) %>% 
  dplyr::filter(assaytype == "Multiplex")


dataInputRawMPXV <-
  dataInputRaw %>% 
  dplyr::select(sampleID_metadata) %>% 
  unique() %>% 
  filter(grepl("P-22", sampleID_metadata)) %>% 
  pull()

##
# Combined meta data
metaDataContained <-
  unique(metadata_MPXV_patients$sampleID_metadata)[(unique(metadata_MPXV_patients$sampleID_metadata) %in% dataInputRawMPXV)]

##
# Generate core sample names as sampleIDs diverge between measured and meta data
metadata_sampleID_core <-
  metadata_MPXV_patients %>% 
  filter(sample_type %in% c("Serum", "Blut")) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  mutate(sampleID_core = str_remove(sampleID_metadata, "-001-[0-9][0-9]"),
         sampleID_core = str_remove(sampleID_core, "P-22-"),
         sampleID_core = str_remove(sampleID_core, "^0+")) %>% 
  select(sampleID_core) %>% 
  unique() 

##
# Pull core sample names
metadata_sampleID_coreID <-
  metadata_sampleID_core %>% 
  select(sampleID_core) %>% 
  unique() %>% 
  pull()

##
# Generate core sample names from measured samples
measured_sampleID_core <-
  dataInputQuantCat %>% 
  filter(grepl("P-22", sampleID_metadata)) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  mutate(sampleID_core = str_remove(sampleID_metadata, "-001-[0-9][0-9]"),
         sampleID_core = str_remove(sampleID_core, "P-22-"),
         sampleID_core = str_remove(sampleID_core, "^0+"),
         sampleID_core = str_remove(sampleID_core, "-0")) # Remove to find sample 15

##
# Pull core sample names from measured samples
measured_sampleID_coreID <-
  measured_sampleID_core %>% 
  select(sampleID_core) %>% 
  unique() %>% 
  pull()

##
# Samples that have been measured but are not matched by metadata
unmatched_sampleID <-
  measured_sampleID_coreID[!(measured_sampleID_coreID %in% 
                               metadata_sampleID_coreID)]

##
# Samples that have not been measured
unmeasured_sampleID <-
  metadata_sampleID_coreID[!(metadata_sampleID_coreID %in% 
                               measured_sampleID_coreID)]

dataInputMVA <- 
  dataInputMVAmeta %>% 
  filter(dilution_assay == 100 & assaytype == "Multiplex") %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                batch,
                panel = panel_corr,
                isotype,
                analyte,
                sampleID_metadata,
                data,
                dataIn,  serostatus, serostatus_cat,
                serostatus.delta, serostatus_cat.delta,
                year_birth, 
                childhoodImmu,
                time_passed = timeToImmu, 
                serostatus.meta,
                serostatus_factor) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         POX_VACC_STATUS_1 = NA_character_,
         MPX_diagnose_final = "No",
         MPX_SYMPTOME = NA_character_,
         MPX_vac_final = if_else(serostatus.meta == 0, "No", "Yes")) %>% 
  filter(!is.na(serostatus.meta))

dataInputMPXVdays <- 
  dataInputMPXVmeta %>% 
  select(caseID_metadata, sampleID_metadata, date_sampleIn) %>% 
  unique() %>% 
  group_by(caseID_metadata) %>% 
  summarize(days_firstsample = as.numeric(difftime(date_sampleIn, min(date_sampleIn)), units="days"),
            n_samples = length(date_sampleIn),
            sampleID_metadata = sampleID_metadata,
            date_sampleIn = date_sampleIn) %>% 
  mutate(time_class = case_when(days_firstsample <= 7 ~ "<=7d",
                                days_firstsample <= 14 ~ "7-14d",
                                days_firstsample > 14 ~ ">14d"),
         time_class = factor(time_class, levels = c("<=7d",
                                                    "7-14d",
                                                    ">14d"), ordered = TRUE)) %>% 
  ungroup() %>% 
  select(sampleID_metadata, days_firstsample)

dataInputMPXV <- 
  dataInputMPXVmeta %>% 
  left_join(dataInputMPXVdays, by = "sampleID_metadata") %>% 
  filter(dilution_assay == 100 & assaytype == "Multiplex") %>% 
  dplyr::select(experiment,
                date,
                plate_assay,  
                batch,
                panel = panel_corr,
                isotype,
                analyte,
                sampleID_metadata,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta,
                serostatus_cat.delta,
                year_birth = date_birthday, 
                POX_VACC_STATUS_1 = vaccination_pox,
                time_passed = days_firstsample) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         year_birth = year(year_birth),
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         childhoodImmu = if_else(year_birth >= 1975, 0, 1),
         MPX_diagnose_final = "Yes",
         MPX_SYMPTOME = NA_character_,
         MPX_vac_final = if_else(serostatus.meta == 0, "No", "Yes"))

dataInputSPox <-
  dataInQuant %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                batch = bead_lot,
                isotype,
                analyte,
                sampleID_metadata = sampleID_meta,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta,
                serostatus_cat.delta,
                MPX_diagnose_final,
                MPX_SYMPTOME,
                MPX_vac_final,
                POX_VACC_STATUS_1) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         panel = "SPox",
         plate_assay = as.character(plate_assay),
         year_birth = NA,
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         childhoodImmu = NA,
         time_passed = NA)


dataInputNKdelta <-
  dataInputQuantCatNK %>% 
  filter(analyte == "Delta") %>% 
  dplyr::select(sampleID_metadata, isotype, serostatus, serostatus_cat) %>% 
  filter(!is.na(serostatus_cat))

dataInputNK <-
  dataInputQuantCatNK %>% 
  left_join(dataInputNKdelta, by = c("sampleID_metadata",
                                     "isotype"), suffix = c("", ".delta")) %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                batch,
                isotype,
                analyte,
                sampleID_metadata,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta,
                serostatus_cat.delta
  ) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         panel = "Pre",
         MPX_diagnose_final = NA,
         MPX_SYMPTOME = NA,
         MPX_vac_final = NA,
         POX_VACC_STATUS_1 = NA,
         plate_assay = as.character(plate_assay),
         year_birth = NA,
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         childhoodImmu = NA,
         time_passed = NA)


# CPXV Input
dataInputCPXVdelta <-
  dataInputQuantCat %>% 
  filter(panel_corr == "CPXV") %>% 
  filter(dilution_assay == 100 & assaytype == "Multiplex") %>% 
  filter(analyte == "Delta") %>% 
  dplyr::select(sampleID_metadata, isotype, serostatus, serostatus_cat) %>% 
  filter(!is.na(serostatus_cat))

dataInputCPXV <-
  dataInputQuantCat %>% 
  filter(panel_corr == "CPXV") %>% 
  filter(dilution_assay == 100 & assaytype == "Multiplex") %>% 
  left_join(dataInputCPXVdelta, by = c("sampleID_metadata",
                                       "isotype"), suffix = c("", ".delta")) %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                batch,
                isotype,
                analyte,
                sampleID_metadata,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta,
                serostatus_cat.delta
  ) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         panel = "CPXV",
         MPX_diagnose_final = NA,
         MPX_SYMPTOME = NA,
         MPX_vac_final = NA,
         POX_VACC_STATUS_1 = NA,
         plate_assay = as.character(plate_assay),
         year_birth = NA,
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         childhoodImmu = NA,
         time_passed = NA)

# New NK panel input
dataInputNKNewdelta <-
  # dataInputQuantCatNK %>% 
  dataInQuantNCCombined %>% 
  filter(analyte == "Delta") %>% 
  dplyr::select(sampleID_metadata = sampleID_meta, isotype, serostatus, serostatus_cat) %>% 
  filter(!is.na(serostatus_cat))

dataInputNKNew <-
  dataInQuantNCCombined %>% 
  left_join(dataInputNKNewdelta, by = c("sampleID_meta" = "sampleID_metadata",
                                        "isotype"), suffix = c("", ".delta"))  %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                isotype,
                analyte,
                sampleID_metadata = sampleID_meta,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta = serostatus.delta.delta,
                serostatus_cat.delta = serostatus_cat.delta.delta,
                childhoodImmu = childhood_immu) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         panel = "Pre_New",
         batch = "Batch_NK",
         MPX_diagnose_final = NA,
         MPX_SYMPTOME = NA,
         MPX_vac_final = NA,
         POX_VACC_STATUS_1 = NA,
         plate_assay = as.character(plate_assay),
         year_birth = NA,
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         time_passed = NA)

# Repeated measurement panel
dataInputSPoxRep <-
  dataInQuantRep %>% 
  dplyr::select(experiment,
                date,
                plate_assay,
                batch = bead_lot,
                isotype,
                analyte,
                sampleID_metadata = sampleID_meta,
                data,
                dataIn,
                serostatus,
                serostatus_cat,
                serostatus.delta,
                serostatus_cat.delta,
                MPX_diagnose_final,
                MPX_SYMPTOME,
                MPX_vac_final,
                POX_VACC_STATUS_1) %>% 
  mutate(data = log10(data),
         dataIn = log10(dataIn),
         panel = "SPox_Rep",
         plate_assay = as.character(plate_assay),
         year_birth = NA,
         serostatus.meta = NA_real_,
         serostatus_factor = NA,
         childhoodImmu = NA,
         time_passed = NA)






dataInput <-
  rbind(dataInputMPXV, dataInputMVA, dataInputSPox, dataInputNK, dataInputCPXV,
        dataInputNKNew, dataInputSPoxRep) %>% 
  mutate(panel_detail = case_when(panel == "MVA" & serostatus.meta == 0 ~ "Pre",
                                  panel == "SPox" & MPX_diagnose_final == "Yes" ~ "MPXV",
                                  panel == "SPox" & MPX_diagnose_final == "No" &
                                    MPX_vac_final == "Yes" ~ "MVA",
                                  panel == "SPox" & MPX_diagnose_final == "No" &
                                    MPX_vac_final == "No" ~ "Pre",
                                  panel == "Pre" ~ "Pre",
                                  panel == "MVA" & serostatus.meta != 0 ~ "MVA",
                                  panel == "CPXV" ~ "CPXV",
                                  panel == "MPXV" ~ "MPXV",
                                  TRUE ~ panel))

##
# Select only patient samples and select necessary columns
dataInputRawFin <-
  dataInput 

dataInputRawFin %>% 
  dplyr::select(sampleID_metadata, panel) %>% 
  unique() %>% 
  dplyr::select(panel) %>% 
  ftable()


dataInputMulticlassed <-
  dataInput %>% 
  mutate(POX_VACC = case_when((childhoodImmu == 0 | POX_VACC_STATUS_1 == "Nein" | POX_VACC_STATUS_1 == "No") ~ "No",
                              (childhoodImmu == 1 | POX_VACC_STATUS_1 == "Ja" | POX_VACC_STATUS_1 == "Yes") ~ "Yes"),
         MPX_vac_final_comb = case_when(panel == "MPXV" ~ "No",
                                   panel == "CPXV" ~ "No",
                                   panel == "Pre_New" ~ "No",
                                   TRUE ~ MPX_vac_final),
         MPX_diagnose_final = case_when(panel == "MVA" ~ "No", 
                                        panel == "Pre" ~ "No",
                                        panel == "Pre_New" ~ "No",
                                        panel == "CPXV" ~ "No",
                                        TRUE ~ MPX_diagnose_final)) %>% 
  select(panel, panel_detail, isotype, analyte, sampleID_metadata, data, dataIn, serostatus,
         serostatus_cat, serostatus.delta, serostatus_cat.delta, POX_VACC, MPX_vac_final_comb,
         MPX_diagnose_final)

dataInputMultiClassComplete <- dataInputMulticlassed[complete.cases(dataInputMulticlassed),]

# Export results for subsequent ML input
export(dataInput, file = "output/dataInputAll.csv")
export(dataInputMultiClassComplete, file = "output/dataInputComplete.csv")
export(dataIn, file = "output/dataPredictionSPox.csv")
