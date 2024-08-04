####
# Generate dataInput for Mpox Analysis
# Daniel Stern RKI
# 22.05.2023
# generate_dataInput.R
####

## Description ----
# The following script imports experimental data and metadata describing 
# the experimental samples and merges them into one dataframe for subsequent
# analysis. 
# To this aim, Excel files containing experiment data from either ELISA or
# Bio-Plex Luminex assays are imported. Metadata containing basic information
# on each sample are also loaded and added to the experimental data.
# The final data frame contains one transformation: the calculation of the 
# difference between VACV and Hep-2 lysates termed delta. Otherwise, no further
# analysis or transformations are performed.
# The aim of this file is to provide one large dataset for all subsequent analysis
# to faciliate data sharing and publication on GitHub without sharing 
# sensitive data.
# A second script will handle all metadata about the samples and bringing it in 
# to a form so that it can be joined easily for subsequent analysis where
# necessary

##
# Changes made to original files before data import
# P-22-01608-001-08 -> Changed to P-22-01610-001-08 name error

## Prepare environment and load packages ----
rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(lubridate)
library(purrr)

## Load functions for data import ----
source("function_read_multiplex.R", encoding = "UTF-8")
source("function_importELISA.R", encoding = "UTF-8")

## Load Multiplex data for validation CPXV MVA panel ----
inputDataFiles <- list.files("input/CPXV_MVA/")[grepl("Pox ", list.files("input/CPXV_MVA/")) &
                                                  !grepl("Meta", list.files("input/CPXV_MVA/"))]
inputMetaFiles <- list.files("input/CPXV_MVA/")[grepl("Pox ", list.files("input/CPXV_MVA/")) &
                                                  grepl("Meta", list.files("input/CPXV_MVA/"))]
outlist <- list()
for(i in 1:length(inputDataFiles)){
  dataout <- read_multiplex(paste("input/CPXV_MVA/", inputDataFiles[i], sep = ""),
                            paste("input/CPXV_MVA/", inputMetaFiles[i], sep = ""))
  outlist[[i]] <- dataout
}

input_CPXV_MVA <- bind_rows(outlist)


## Load and bind sampleID_metadata
sampleID_meta_CPXV_MVA <- import("input/CPXV_MVA/SampleID_assay_metadata_code.xlsx")

input_CPXV_MVA_meta <-
  input_CPXV_MVA %>% 
  left_join(sampleID_meta_CPXV_MVA, by = c("experiment", "sampleID_assay"))  %>% 
  mutate(panel = "CPXV_MVA")

## Delete non needed variables
rm(dataout, input_CPXV_MVA, outlist, sampleID_meta_CPXV_MVA,
   i, inputDataFiles, inputMetaFiles)


## Load Multiplex data for MPXV panel ----
inputDataFiles <- list.files("input/MPXV/")[grepl("Pox_", list.files("input/MPXV/")) &
                                              !grepl("Meta", list.files("input/MPXV/"))]
inputMetaFiles <- list.files("input/MPXV/")[grepl("Pox_", list.files("input/MPXV/")) &
                                              grepl("Meta", list.files("input/MPXV/"))]
outlist <- list()
for(i in 1:length(inputDataFiles)){
  dataout <- read_multiplex(paste("input/MPXV/", inputDataFiles[i], sep = ""),
                            paste("input/MPXV/", inputMetaFiles[i], sep = ""))
  outlist[[i]] <- dataout
}

input_MPXV <- bind_rows(outlist)


## Load and bind sampleID_metadata
path <- c("input/MPXV/Proben-Listen_Sample_Meta.xlsx")

# Read all excel worksheets using the purrr package
rawDataInputList <-
  path %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  purrr:::map(readxl::read_excel, path = path)


namesList <- names(rawDataInputList)

metadataAssayRaw <-
  do.call(rbind, rawDataInputList) 

metadataAssay <-
  metadataAssayRaw %>% 
  cbind(name = row.names(metadataAssayRaw)) %>%
  mutate(name = str_remove(name, "Teil "),
         name = str_remove(name, "^ "),
         sample = as.integer(Nummer),
         plate_assay = as.integer(str_remove(Platte, "Platte "))) %>% 
  separate(name, into = c("assay", "date_meta"), sep = " ") %>% 
  mutate(date = as.POSIXct(date_meta, format = "%d.%m.%Y", tz = "UTC"),
         # date_character = as.character(date),
         sample = if_else(plate_assay == 2, as.integer(sample - 12), sample),
         assay = as.numeric(assay)) %>% 
  select(-Nummer, -date_meta) %>% 
  filter(!is.na(sample)) %>% 
  mutate(sampleID_assay = paste("S", sample, sep = ""),
         experiment = paste("Pox_25_", assay, sep = ""),
         plate = str_replace(Platte, "Platte ", "plate_")) %>% 
  select(experiment, plate, sampleID_assay, sampleID_metadata = Probe)

row.names(metadataAssay) <- NULL

input_MPXV_meta <-
  input_MPXV %>% 
  left_join(metadataAssay, by = c("experiment", "sampleID_assay", "plate")) %>% 
  mutate(sampleID_metadata = if_else(sampleID_assay == "NK", 
                                     "SeCoV-1036", sampleID_metadata))  %>% 
  mutate(panel = "MPXV") %>% 
  filter(!is.na(data)) # Workaround

## Delete non needed variables
rm(dataout, input_MPXV, outlist, metadataAssay, metadataAssayRaw,
   i, inputDataFiles, inputMetaFiles, namesList, path, rawDataInputList)


## Load Multiplex data for NK panel ----
inputDataFiles <- list.files("input/NK/")[grepl("Pox_", list.files("input/NK/")) &
                                              !grepl("Meta", list.files("input/NK/"))]
inputMetaFiles <- list.files("input/NK/")[grepl("Pox_", list.files("input/NK/")) &
                                              grepl("Meta", list.files("input/NK/"))]
outlist <- list()
for(i in 1:length(inputDataFiles)){
  dataout <- read_multiplex(paste("input/NK/", inputDataFiles[i], sep = ""),
                            paste("input/NK/", inputMetaFiles[i], sep = ""))
  outlist[[i]] <- dataout
}

input_NK <- bind_rows(outlist)

## Load Metadata for NK panel ----
metadataAssayNK <- import("input/NK/Proben-Listen_Sample_Meta_NK.xlsx") %>% 
  mutate(sampleID_assay = paste("S", Nummer, sep = "")) %>% 
  select(sampleID_assay, sampleID_metadata = Probe)

input_NK_meta <-
  input_NK %>% 
  left_join(metadataAssayNK, by = "sampleID_assay") %>% 
  mutate(panel = "NK")


rm(dataout, outlist, input_NK, metadataAssayNK,
   i, inputDataFiles, inputMetaFiles)


## Load ELISA data for validation CPXV MVA panel ----
## Load sampleID
sampleIDELISA <- import("input/ELISA_CPXV_MVA/Probenliste.xlsx") 

outlist <- list()
outlist[[1]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221014_IgG_Platte 1.xlsx",
                            "221014",
                            "plate_1",
                            "IgG", 4)
outlist[[2]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221014_IgG_Platte 2.xlsx",
                            "221014",
                            "plate_2",
                            "IgG", 4)
outlist[[3]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221014_IgG_Platte 3.xlsx",
                            "221014",
                            "plate_3",
                            "IgG", 4)
outlist[[4]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221014_IgG_Platte 4.xlsx",
                            "221014",
                            "plate_4",
                            "IgG", 4)
outlist[[5]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221021_IgM_Platte 1.xlsx",
                            "221021",
                            "plate_1",
                            "IgM", 4)
outlist[[6]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221021_IgM_Platte 2.xlsx",
                            "221021",
                            "plate_2",
                            "IgM", 4)
outlist[[7]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221021_IgM_Platte 3.xlsx",
                            "221021",
                            "plate_3",
                            "IgM", 4)
outlist[[8]] <- importELISA("Pox_10",
                            "input/ELISA_CPXV_MVA/221021_IgM_Platte 4.xlsx",
                            "221021",
                            "plate_4",
                            "IgM", 4)
input_ELISA_CPXV_MVA <- bind_rows(outlist)


## Load ELISA metadata
inputELISAantigen <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                            range = "A2:M10") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "analyte") %>%
  mutate(well = paste(`...1`, col, sep = "")) %>%
  dplyr::select(well, analyte)

inputELISAdilutions <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                              range = "A14:M22") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "dilution") %>%
  mutate(well = paste(`...1`, col, sep = "")) %>%
  dplyr::select(well, dilution)

outlist <- list()
outlist[[1]] <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                       range = "A26:M34") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_1") %>%
  dplyr::select(well, plate, sampleID_assay)

outlist[[2]] <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                       range = "A38:M46") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_2") %>%
  dplyr::select(well, plate, sampleID_assay)

outlist[[3]] <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                       range = "A50:M58") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_3") %>%
  dplyr::select(well, plate, sampleID_assay)

outlist[[4]] <- import("input/ELISA_CPXV_MVA/Plattenlayout_ELISA.xlsx",
                       range = "A62:M70") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_4") %>%
  dplyr::select(well, plate, sampleID_assay)

input_ELISA_layout <- bind_rows(outlist)

## Bind sampleID_metadata
input_ELISA_CPXV_MVA_meta <-
  input_ELISA_CPXV_MVA %>% 
  left_join(inputELISAantigen, by = c("well")) %>% 
  left_join(inputELISAdilutions, b = c("well")) %>% 
  left_join(input_ELISA_layout, by = c("well", "plate")) %>% 
  left_join(sampleIDELISA, by = c("sampleID_assay")) %>% 
  mutate(assaytype = "ELISA",
         batch = "ELISA") %>% 
  select(experiment, filename, date, plate, batch, well, assaytype,
         sampleID_assay, dilution, isotype, analyte, data,
         sampleID_metadata) %>% 
  mutate(panel = "CPXV_MVA")

inputDelta <-
  input_ELISA_CPXV_MVA_meta %>% 
  select(- well) %>% 
  pivot_wider(names_from = analyte, values_from = data) %>%
  mutate(Delta = VACV - Hep2,
         Delta = Delta + abs(min(Delta))) %>%
  pivot_longer(cols = c(VACV, Hep2, Delta), names_to = "analyte", values_to = "data") %>% 
  filter(analyte == "Delta") %>% 
  mutate(well = "delta")

input_ELISA_CPXV_MVA_meta <-
  input_ELISA_CPXV_MVA_meta %>% 
  rbind(inputDelta)

## Delete non needed variables
rm(input_ELISA_CPXV_MVA, input_ELISA_layout,
   outlist, sampleIDELISA, inputELISAantigen, 
   inputELISAdilutions, inputDelta)


#### dev 
## Load ELISA data for validation MPXV panel ----
## Load sampleID
sampleIDELISA <- import("input/ELISA_MPXV/Probenliste.xlsx") 

outlist <- list()
outlist[[1]] <- importELISA("Pox_28",
                            "input/ELISA_MPXV/230209_IgG_Platte 1_WH.xlsx",
                            "230209",
                            "plate_1",
                            "IgG",5)
outlist[[2]] <- importELISA("Pox_28",
                            "input/ELISA_MPXV/230209_IgG_Platte 2_WH.xlsx",
                            "230209",
                            "plate_2",
                            "IgG",5)
outlist[[3]] <- importELISA("Pox_28",
                            "input/ELISA_MPXV/230216_IgM_Platte 1.xlsx",
                            "221014",
                            "plate_1",
                            "IgM", 5)
outlist[[4]] <- importELISA("Pox_28",
                            "input/ELISA_MPXV/230216_IgM_Platte 1.xlsx",
                            "230216",
                            "plate_2",
                            "IgM", 5)

input_ELISA_MPXV <- bind_rows(outlist)


## Load ELISA metadata
inputELISAantigen <- import("input/ELISA_MPXV/Plattenlayout_ELISA.xlsx",
                            range = "A2:M10") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "analyte") %>%
  mutate(well = paste(`...1`, col, sep = "")) %>%
  dplyr::select(well, analyte)

inputELISAdilutions <- import("input/ELISA_MPXV/Plattenlayout_ELISA.xlsx",
                              range = "A14:M22") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "dilution") %>%
  mutate(well = paste(`...1`, col, sep = "")) %>%
  dplyr::select(well, dilution)

outlist <- list()
outlist[[1]] <- import("input/ELISA_MPXV/Plattenlayout_ELISA.xlsx",
                       range = "A26:M34") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_1") %>%
  dplyr::select(well, plate, sampleID_assay)

outlist[[2]] <- import("input/ELISA_MPXV/Plattenlayout_ELISA.xlsx",
                       range = "A38:M46") %>%
  pivot_longer(cols = -`...1`, names_to = "col", values_to = "sampleID_assay") %>%
  mutate(well = paste(`...1`, col, sep = ""), 
         plate = "plate_2") %>%
  dplyr::select(well, plate, sampleID_assay)

input_ELISA_layout <- bind_rows(outlist)

## Bind sampleID_metadata
input_ELISA_MPXV_meta <-
  input_ELISA_MPXV %>% 
  left_join(inputELISAantigen, by = c("well")) %>% 
  left_join(inputELISAdilutions, b = c("well")) %>% 
  left_join(input_ELISA_layout, by = c("well", "plate")) %>% 
  left_join(sampleIDELISA, by = c("sampleID_assay")) %>% 
  mutate(assaytype = "ELISA",
         batch = "ELISA") %>% 
  select(experiment, filename, date, plate, batch, well, assaytype,
         sampleID_assay, dilution, isotype, analyte, data,
         sampleID_metadata) %>% 
  mutate(panel = "MPXV") 

inputDelta <-
  input_ELISA_MPXV_meta %>% 
  select(- well) %>% 
  pivot_wider(names_from = analyte, values_from = data) %>%
  mutate(Delta = VACV - Hep2,
         Delta = Delta + abs(min(Delta))) %>%
  pivot_longer(cols = c(VACV, Hep2, Delta), names_to = "analyte", values_to = "data") %>% 
  filter(analyte == "Delta") %>% 
  mutate(well = "delta")

input_ELISA_MPXV_meta <-
  input_ELISA_MPXV_meta %>% 
  rbind(inputDelta)


## Delete non needed variables
rm(input_ELISA_MPXV, input_ELISA_layout,
   outlist, sampleIDELISA, inputELISAantigen, 
   inputELISAdilutions, inputDelta)

dataInput <-
  rbind(input_CPXV_MVA_meta,
        input_ELISA_CPXV_MVA_meta,
        input_MPXV_meta,
        input_ELISA_MPXV_meta)

dataInputNK <- input_NK_meta
save(dataInput, file = "output/dataInput.Rdata")
save(dataInputNK, file = "output/dataInputNK.Rdata")





