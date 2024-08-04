####
# Generate metadataInput for Mpox Analysis
# Daniel Stern RKI
# 20.02.2023
# generate_metadataInput.R
####

## Description ----
# The following script imports metadata describing sample properties. 

## Prepare environment and load packages ----
rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(lubridate)

## Load Metadata for CPXV MVA panel: 1 - Load IFA data ----
# Data is generated from IFA experiments on either CPXV or VACV slides
load("input/Metadata_CPXV_MVA/dataInputIFA.RData")

metadata_IFA_CPXV_MVA <-
  dataInputIFA %>% 
  filter(!is.na(titer)) %>% 
  filter(!grepl("Humanes", sample)) %>% 
  filter(!(grepl("DS", sample) | grepl("SeCoV", sample))) %>% 
  select(sampleID_metadata = sample, isotype, titer) %>% 
  mutate(panel = "CPXV_MVA")


## Load Metadata for CPXV MVA panel: 2 - Load MVA vaccination data ----
metadata_MVA <- import("input/Metadata_CPXV_MVA/2022_10_21_Imvanex_Studie_anonym.xlsx") %>%
  filter(ImpfungDurch != "Nein") %>% 
  pivot_longer(cols = contains("Imvanex"), names_to = "Impfung",
               values_to = "DatumImpfung")   %>% 
  filter(!is.na(DatumImpfung)) %>% 
  mutate(FallID = str_remove(MaterialID, ".{3}$"),
         Impfung = as.numeric(str_extract(Impfung, "[0-9]")),
         DatumImpfung = if_else(grepl("^44", DatumImpfung) | grepl("^43", DatumImpfung),
                                as_date(as.numeric(DatumImpfung), origin = ("1899-12-30")),
                                as_date(DatumImpfung, format = ("%d.%m.%Y"))),
         Abnahmedatum = as_date(Abnahmedatum, format = ("%d.%m.%Y")),
         childhoodImmu = if_else(Geburtsjahr < 1975, 1, 0),
         timeID = as.numeric(sub('.*(?=.{2}$)', '', MaterialID, perl=T)),
         timeToImmu = Abnahmedatum - DatumImpfung,
         serostatus = case_when(timeToImmu == 0 & Impfung == 1 ~ 0,
                                timeToImmu > 0 & Impfung == 1 ~ 1,
                                timeToImmu == 0 & Impfung == 2 ~ 1,
                                timeToImmu > 0 & Impfung == 2 ~ 2,
                                timeToImmu == 0 & Impfung == 3 ~ 2,
                                timeToImmu > 0 & Impfung == 3 ~ 3)) %>% 
  filter(!(MaterialID == "S-20-004-030-01" & Impfung == 1))


timeSera <-
  metadata_MVA %>% 
  select(FallID, timeID, Abnahmedatum) %>% 
  pivot_wider(names_from = timeID, values_from = Abnahmedatum) %>% 
  mutate(time_1 = `1`- `1`,
         time_2 = `2`-`1`,
         time_3 = `3` - `1`,
         time_4 = `4` - `1`,
         time_5 = `5` - `1`) %>% 
  select(FallID, starts_with("time")) %>% 
  pivot_longer(cols = starts_with("time"), values_to = "days", names_to = "timeID",
               names_prefix = "time_") %>% 
  mutate(timeID = as.numeric(timeID)) %>% 
  filter(!is.na(days))

metadata_MVA_time <-
  metadata_MVA %>% 
  left_join(timeSera, by = c("FallID", "timeID")) %>% 
  mutate(serostatus_fine = if_else(serostatus == 2 & days > 500, 
                                   "2 late", as.character(serostatus)),
         serostatus_factor = factor(serostatus_fine,
                                    levels = c("0", "1", "2", "2 late", "3"),
                                    labels = c("Pre",
                                               "Prime", 
                                               "First Boost", 
                                               "Pre Second Boost", 
                                               "Second Boost"), ordered = TRUE)) %>% 
  select(caseID_metadata = FallID, 
         sampleID_metadata = MaterialID,
         timeID,
         date_serum = Abnahmedatum,
         date_vaccination = DatumImpfung,
         year_birth = Geburtsjahr,
         childhoodImmu,
         timeToImmu,
         days_firstsample = days,
         serostatus,
         serostatus_factor)

## Delete non needed variables
rm(dataInputIFA, metadata_MVA, timeSera, timeToImmu, timeVacc)


## Load metadata for MPXV panel: 3 - Load patient meta data -
metadata_MPXV_PCR <- import("input/Metadata_MPXV/104 Blutproben_Virämietestung.xlsx") %>% 
  filter(!is.na(MaterialID)) %>% 
  select(-`...7`, -`...8`) %>% 
  mutate(date_posPCR = if_else(is.na(`Entnahme pos. PCR`),
                               as.POSIXct(str_remove(Kommentar, "PE am "), 
                                          format = "%d.%m.%Y"),
                               as.POSIXct(`Entnahme pos. PCR`)),
         date_serum = if_else(is.na(`Entnahme Blut`),
                              as.POSIXct(str_remove(Kommentar, "PE am "), 
                                         format = "%d.%m.%Y"),
                              as.POSIXct(`Entnahme Blut`))) %>% 
  select(caseID_metadata = FallID, 
         sampleID_metadata = MaterialID,
         date_posPCR, date_serum, timeToPCR = Differenz)


## Large file containing all metadata ----
metadata_MPXV_patients <- import("input/Metadata_MPXV/MPXVSerenMultiplexEtablierung.xlsx",
                                 sheet = 2) %>% 
  mutate(caseID_metadata = FallID,
         sampleID_metadata = MaterialID,
         date_sampleIn = as.POSIXct(Einsendedatum, format = "%d.%m.%Y"),
         date_birthday = as.POSIXct(Geburtsdatum, format = "%d.%m.%Y"),
         vaccination_pox = `Impfung (Pocken)`,
         sample_type = Probenmaterial,
         results_PCR_OPV = `Erstuntersuchung: PCR Orthopockenviren`,
         CT_PCR_OPV = as.numeric(str_replace(`PCR Orthopockenviren: Erstuntersuchung: Wert`, ",", ".")),
         CT_PCR_c_myc = as.numeric(str_replace(`PCR c-myc: Erstuntersuchung: Wert`, ",", ".")),
         CT_PCR_KoMa = as.numeric(str_replace(`PCR KoMa: Erstuntersuchung: Wert`, ",", ".")),
         results_PCR_MPXV = `Erstuntersuchung: PCR MPXV G2R generisch`,
         results_PCR_MPXV_WA = `Erstuntersuchung: PCR MPXV G2R WA`,
         CT_PCR_MPXV = as.numeric(str_replace(`PCR MPXV G2R generisch: Erstuntersuchung: Wert`, ",", ".")),
         CT_PCR_MPXV_WA = as.numeric(str_replace(`PCR MPXV G2R WA: Erstuntersuchung: Wert`, ",", ".")),
         titer_IgG = if_else(`IFT Orthopockenviren IgG: Erstuntersuchung: Wert` == "-", 0,
                             as.numeric(`IFT Orthopockenviren IgG: Erstuntersuchung: Wert`)),
         titer_IgM = if_else(`IFT Orthopockenviren IgM: Erstuntersuchung: Wert` == "-", 0,
                             as.numeric(`IFT Orthopockenviren IgM: Erstuntersuchung: Wert`))) %>% 
  select(caseID_metadata:titer_IgM)


## Smaller metadatafile containing only IFA data ----
metadata_IFA_MPXV_patients <-
  metadata_MPXV_patients %>% 
  filter(!is.na(titer_IgG)) %>% 
  select(sampleID_metadata, titer_IgG, titer_IgM) %>% 
  pivot_longer(cols = c(starts_with("titer")), names_to = "isotype",
               names_prefix = "titer_", values_to = "titer") %>% 
  mutate(panel = "MPXV")

## Combine IFA data from different panels ----
metadata_IFA <-
  metadata_IFA_CPXV_MVA %>% 
  rbind(metadata_IFA_MPXV_patients)

## Metadata with only the patient sample
metadata_MPXV_patients_only <-
  metadata_MPXV_patients %>% 
  select(caseID_metadata, date_birthday, vaccination_pox) %>% 
  unique()

## Save metadata
save(metadata_IFA, file = "output/metadata_IFA.Rdata")
save(metadata_MVA_time, file = "output/metadata_MVA_time.Rdata")
save(metadata_MPXV_patients_only, file = "output/metadata_MPXV_patients_only.Rdata")
save(metadata_MPXV_patients, file = "output/metadata_MPXV_patients.Rdata")
