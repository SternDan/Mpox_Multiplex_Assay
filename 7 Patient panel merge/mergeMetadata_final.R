##
# Analyse patient panel
# Daniel Stern
# 2024-08-03
# Robert Koch Institute
# Version Final
# Data wrangling for metadata files
##

##
# Clean environment
rm(list = ls(all.names = TRUE))

##
# Load packages
library(rio)
library(tidyverse)


##
# Load data
load("../4 Method Comparison/output/dataInputQuantCat.Rdata")
load("../2 Import metadata/output/metadata_MVA_time.Rdata")
load("../2 Import metadata/output/metadata_MPXV_patients.Rdata")


##
# Add serostatus based on Delta signals at 1:100 dilution
dataInputQuantCatDelta <-
  dataInputQuantCat %>% 
  filter(dilution_assay == 100 & analyte == "Delta" & assaytype == "Multiplex") %>% 
  dplyr::select(experiment, plate_assay, isotype, sampleID_metadata, dataIn,
                serostatus, serostatus_cat)

dataInputQuantCatDeltaMerge <-
  dataInputQuantCat %>% 
  left_join(dataInputQuantCatDelta, by = c("experiment", "plate_assay", "isotype", "sampleID_metadata"),
            suffix = c("", ".delta"))

##
# Determine PCR positivity in metadata samples
# If PCR is positive in at least on sample type -> sample is considered 
# positive. 
# Determine PCR positivity in case ID
# If at least one sample is positive -> Whole case is considered positive.
PCR_sampleStatus <-
  metadata_MPXV_patients %>% 
  select(caseID_metadata, sampleID_metadata, date_sampleIn, sample_type, results_PCR_OPV, results_PCR_MPXV) %>% 
  mutate(PCR_OPV_num = case_when(results_PCR_OPV == "positiv" ~ 1,
                                 results_PCR_OPV == "negativ" ~ 0,
                                 results_PCR_OPV == TRUE ~ NA_real_),
         PCR_MPXV_num = case_when(results_PCR_MPXV == "positiv" ~ 1,
                                  results_PCR_MPXV == "negativ" ~ 0,
                                  results_PCR_MPXV == TRUE ~ NA_real_)) %>% 
  group_by(caseID_metadata, date_sampleIn) %>% 
  summarise(PCR_OPV_sampleStatus = if_else(sum(PCR_OPV_num, na.rm = TRUE) >= 1, 1, 0),
            PCR_MPXV_sampleStatus = if_else(sum(PCR_MPXV_num, na.rm = TRUE) >= 1, 1, 0))

PCR_caseStatus <-
  PCR_sampleStatus %>% 
  group_by(caseID_metadata) %>% 
  summarise(PCR_OPV_caseStatus= if_else(sum(PCR_OPV_sampleStatus, na.rm = TRUE) >= 1, 1, 0),
            PCR_MPXV_caseStatus = if_else(sum(PCR_MPXV_sampleStatus, na.rm = TRUE) >= 1, 1, 0))

metadata_MPXV <-
  metadata_MPXV_patients %>% 
  left_join(PCR_sampleStatus, by = c("caseID_metadata", "date_sampleIn")) %>% 
  left_join(PCR_caseStatus, by = c("caseID_metadata"))


## 
# Test which samples have been tested and if tested samples can be found 
# in metadata table

##
# Pull unique sampleIDs from metadata
metadata_sampleID <-
  metadata_MPXV %>% 
 # filter(sample_type %in% c("Serum", "Blut")) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()

##
# Pull unique sampleIDs from measured samples
measured_sampleID <-
  dataInputQuantCat %>% 
  filter(grepl("P-22", sampleID_metadata)) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()


##
# Generate core sample names as sampleIDs diverge between measured and meta data
metadata_sampleID_core <-
  metadata_MPXV %>% 
  filter(sample_type %in% c("Serum", "Blut")) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  mutate(sampleID_core = str_remove(sampleID_metadata, "-001-[0-9][0-9]"),
         sampleID_core = str_remove(sampleID_core, "P-22-"),
         sampleID_core = str_remove(sampleID_core, "^0+"))

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


##
# Join coreIDs to both the metadata and the measured data for subsequent merge
dataInputMPXV <-
  dataInputQuantCatDeltaMerge %>% 
  filter(panel_corr == "MPXV") %>% 
  left_join(measured_sampleID_core, by = c("sampleID_metadata")) %>% 
  filter(!is.na(sampleID_core)) %>% 
  mutate(last_sampleID = sub('.*(?=.{2}$)', '', sampleID_metadata, perl=T))

metadata_MPXV_core <-
  metadata_MPXV %>% 
  left_join(metadata_sampleID_core, by = c("sampleID_metadata")) %>% 
  mutate(last_sampleID = sub('.*(?=.{2}$)', '', sampleID_metadata, perl=T))

dataInputMPXVmeta <-
  dataInputMPXV %>% 
  left_join(metadata_MPXV_core, by = c("sampleID_core", "last_sampleID"), suffix = c("", ".meta"))

unique(dataInputMPXVmeta$sampleID_metadata)


##
# Clean environment
rm(dataInputMPXV, metadata_MPXV, metadata_MPXV_patients,
   metadata_sampleID_core, PCR_caseStatus, PCR_sampleStatus,
   measured_sampleID, metadata_sampleID, unmatched_sampleID, unmeasured_sampleID)

##
# Join metadata with data for MVA panel
sampleID_meta <-
  metadata_MVA_time %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()

sampleID_measured <-
  dataInputQuantCatDeltaMerge %>% 
  filter(panel_corr == "MVA") %>% 
  filter(grepl("S-", sampleID_metadata)) %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()

##
# Samples that have been measured but are not in the metadata file
sampleID_measured[!(sampleID_measured %in% sampleID_meta)]
# Only preimmun sera without subsequent immunisation


##
# Sample in the metadata that have not been measured
sampleID_meta[!(sampleID_meta %in% sampleID_measured)]
# Samples have been excluded and should be repeated

##
dataInputMVAmeta <-
  dataInputQuantCatDeltaMerge %>% 
  filter(panel_corr == "MVA") %>% 
  filter(grepl("S-", sampleID_metadata)) %>% 
  left_join(metadata_MVA_time, by = c("sampleID_metadata"), suffix = c("", ".meta"))
length(unique(dataInputMVAmeta$sampleID_metadata))

save(dataInputMPXVmeta, file = "output/dataInputMPXVmeta_all.Rdata")
save(dataInputMVAmeta, file = "output/dataInputMVAmeta.Rdata")
