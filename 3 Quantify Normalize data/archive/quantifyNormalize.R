##
# Quantify and normalize Multiplex and ELISA data
# Date 06.03.2023
# Daniel Stern
# Robert Koch Institute
# Version 1.1
##

##
# Clean environment
rm(list = ls(all.names = TRUE))

##
# Load packages
library(rio)
library(tidyverse)
library(drLumi)

##
# Source function used to quantify the results using the drLumi package
source("functions_quantify.R", encoding = "UTF-8")

##
# Load data
load("../1 Import data/output/dataInput.Rdata")
load("../1 Import data/output/dataInputNK.Rdata")

##
# Load ec-file
ecdata <- import("input/EC_input.xlsx") %>% 
  pivot_longer(c(-sample), names_to = "analyte", values_to = "ec")

##
# Normalize data

##
# Define samples to be excluded due to errors 
# To Do: Include plate and date to be able to include new measurements
excludeSamples <- c("P-22-01381-001-18", "P-22-013-0-001-03", 
                    "P-22-013-0-001-09", "P-22-01386-001-08",
                    "S-20-004-006-01", "S-20-004-006-02", "S-20-004-006-03")

# Select data containing the results for "Standard1"
standard <- dataInput %>%
  filter(sampleID_assay == "Standard1") %>%
  dplyr::select(-filename, -date, -well, -assaytype,
                -sampleID_assay, -dilution, -sampleID_metadata, - panel)

# 
dataInputQC <-
  dataInput %>%
  filter(assaytype == "Multiplex" & analyte == "HSA", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, sampleID_metadata, isotype, dilution) %>% 
  summarise(highBG = if_else(data > 500, T, F))

dataInputQCIgG <-
  dataInput %>%
  filter(assaytype == "Multiplex" & analyte == "ahIgG", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, sampleID_metadata, isotype, dilution) %>% 
  summarise(lowIgG = if_else(data < 10000, T, F))


# Normalize data
dataInputNorm <-
  dataInput %>%
  left_join(dataInputQC, by = c("experiment", "plate","sampleID_assay", "sampleID_metadata", "isotype",
                                "dilution")) %>%
  left_join(dataInputQCIgG, by = c("experiment", "plate","sampleID_assay", "sampleID_metadata", "isotype",
                                   "dilution")) %>%
  left_join(standard, by = c("experiment", "plate", "isotype",
                             "analyte"),
            suffix = c("", ".std")) %>%
  mutate(data.norm = log10(data/data.std)+3) %>% 
  mutate(exp_plate = paste(experiment, plate, sep = ".")) %>% 
  unique()

rm(dataInput, dataInputQC, dataInputQCIgG, standard)

multiInputList <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgG") %>% 
  group_split(experiment, plate_assay)

names(multiInputList) <- unique(dataInputNorm$exp_plate[dataInputNorm$assaytype == "Multiplex"])

outlist <- list()
for(i in 1:length(multiInputList)){
  outlist[[i]] <- quantifyMultiplex(ecdata = ecdata, multiInputList = multiInputList)
}

names(outlist) <- names(multiInputList)

# Error with fit
# Pox_15.plate_2 -> use standard from Pox_15.plate_1
# Pox_25_10.plate_2 -> standard from Pox_25_10.plate_1
# Pox_25_3.plate_2 -> Exclude


#summary(outlist[[6]]$standards)
#summary(outlist[[13]]$standards)
#summary(outlist[[25]]$standards)

#names(outlist[6])
#names(outlist[13])
#names(outlist[25])

## Re-Fit using standard curves from plate 1
outlistRe <- list()
outlistRe[[6]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputList,
                         allanalytesIn = outlist[[5]]$standards, i = 6)

outlistRe[[13]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputList,
                         allanalytesIn = outlist[[12]]$standards, i = 13)
outlistFinal <-
  outlist

outlistFinal[[6]] <- outlistRe[[6]]
outlistFinal[[13]] <- outlistRe[[13]]
names(outlistFinal) <- names(multiInputList)
outlistFinal <- outlistFinal[-25]

#
liststandards <- lapply(outlistFinal, `[[`, 1)
listdata <- lapply(outlistFinal, `[[`, 2)


## Fit IgM with IgG Standards
multiInputListIgM <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgM") %>% 
  group_split(experiment, plate_assay)


outlistIgM <- list()
for(i in 1:length(multiInputListIgM)){
  outlistIgM[[i]] <- quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgM,
                                            allanalytesIn = outlist[[i]]$standards, i)
  # outlist[i] <- multiInputList[[i]]$experiment
}

## Re-Fit using standard curves from plate 1
outlistReIgM <- list()
outlistReIgM[[6]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgM,
                         allanalytesIn = outlist[[5]]$standards, i = 6)

outlistReIgM[[13]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgM,
                         allanalytesIn = outlist[[12]]$standards, i = 13)
outlistFinalIgM <-
  outlistIgM

outlistFinalIgM[[6]] <- outlistReIgM[[6]]
outlistFinalIgM[[13]] <- outlistReIgM[[13]]
names(outlistFinalIgM) <- names(multiInputList)
outlistFinalIgM <- outlistFinalIgM[-25]

listdataIgM <- lapply(outlistFinalIgM, `[[`, 2)

## Fit IgA with IgG Standards
multiInputListIgA <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgA") %>% 
  group_split(experiment, plate_assay)


outlistIgA <- list()
for(i in 1:length(multiInputListIgA)){
  outlistIgA[[i]] <- quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgA,
                                            allanalytesIn = outlist[[i]]$standards, i)
  # outlist[i] <- multiInputList[[i]]$experiment
}

## Re-Fit using standard curves from plate 1
outlistReIgA <- list()
outlistReIgA[[6]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgA,
                         allanalytesIn = outlist[[5]]$standards, i = 6)

outlistReIgA[[13]] <-
  quantifyMultiplexReFit(ecdata = ecdata, multiInputList = multiInputListIgA,
                         allanalytesIn = outlist[[12]]$standards, i = 13)
outlistFinalIgA <-
  outlistIgA

outlistFinalIgA[[6]] <- outlistReIgA[[6]]
outlistFinalIgA[[13]] <- outlistReIgA[[13]]
names(outlistFinalIgA) <- names(multiInputList)
outlistFinalIgA <- outlistFinalIgA[-25]

listdataIgA <- lapply(outlistFinalIgA, `[[`, 2)



## Analyse ELISA
dataInputNorm %>% 
  filter(assaytype == "ELISA") %>%
  filter(!grepl("Standard", sampleID_assay)) %>% 
  select(dilution) %>% 
  unique()




##
# Select on needed data for ELISA analysis and mutate OD to 10^OD
# as drLumi perfomrs log-transformation


##
# Select ec-data only for delta
ecdataELISA <-
  ecdata %>% 
  filter(analyte == "Delta")

## Analyse ELISA IgG
multiInputListIgGELISA <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "ELISA" & isotype == "IgG") %>% 
  filter(analyte == "Delta") %>% 
  mutate(data = 10^data) %>% 
  group_split(experiment, plate_assay)


outlistELISAIgG <- list()
for(i in 1:length(multiInputListIgGELISA)){
  outlistELISAIgG[[i]] <- quantifyELISA(ecdata = ecdataELISA, multiInputList = multiInputListIgGELISA)
}

liststandardsELISAIgG <- lapply(outlistELISAIgG, `[[`, 1)
listdataELISAIgG <- lapply(outlistELISAIgG, `[[`, 2)

## Analyse ELISA IgM
multiInputListIgMELISA <-
  dataInputNorm %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "ELISA" & isotype == "IgM") %>% 
  filter(analyte == "Delta") %>% 
  mutate(data = 10^data) %>% 
  group_split(experiment, plate_assay)


outlistELISAIgM <- list()
for(i in 1:length(multiInputListIgMELISA)){
  outlistELISAIgM[[i]] <- quantifyMultiplexReFitELISA(ecdata = ecdata, multiInputList = multiInputListIgMELISA,
                                                 allanalytesIn = outlistELISAIgG[[i]]$standards, i)
}

liststandardsELISAIgM <- lapply(outlistELISAIgM, `[[`, 1)
listdataELISAIgM <- lapply(outlistELISAIgM, `[[`, 2)

dataInputQuant <-
  rbind(do.call(rbind,listdata),
        do.call(rbind,listdataIgA),
        do.call(rbind,listdataIgM),
        do.call(rbind,listdataELISAIgG),
        do.call(rbind,listdataELISAIgM)) %>% 
  mutate(data = if_else(assaytype == "ELISA", log10(data), data)) %>% 
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
  mutate(panel_corr = case_when(panel == "CPXV_MVA" & grepl("S-", sampleID_metadata) ~ "MVA",
                                panel == "CPXV_MVA" & grepl("P-", sampleID_metadata) ~ "CPXV",
                                panel == "MPXV" ~ "MPXV",)) %>%
  ungroup() %>% 
  filter(!(sampleID_metadata %in% excludeSamples))

save(dataInputQuant, file = "output/dataInputQuant.Rdata")


##
# Quantify NK-Panel
# Select data containing the results for "Standard1"
standardNK <- dataInputNK %>%
  filter(sampleID_assay == "Standard1") %>%
  dplyr::select(-filename, -date, -well, -assaytype,
                -sampleID_assay, -dilution, -sampleID_metadata, - panel)

# 
dataInputQCNK <-
  dataInputNK %>%
  filter(assaytype == "Multiplex" & analyte == "HSA", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, sampleID_metadata, isotype, dilution) %>% 
  summarise(highBG = if_else(data > 500, T, F))

dataInputQCIgGNK <-
  dataInputNK %>%
  filter(assaytype == "Multiplex" & analyte == "ahIgG", dilution == 100) %>% 
  group_by(experiment, plate, sampleID_assay, sampleID_metadata, isotype, dilution) %>% 
  summarise(lowIgG = if_else(data < 10000, T, F))


# Normalize data
dataInputNormNK <-
  dataInputNK %>% 
  filter(sampleID_assay != "NK") %>% 
  left_join(dataInputQCNK, by = c("experiment", "plate","sampleID_assay", "sampleID_metadata", "isotype",
                                "dilution")) %>%
  left_join(dataInputQCIgGNK, by = c("experiment", "plate","sampleID_assay", "sampleID_metadata", "isotype",
                                   "dilution")) %>%
  left_join(standardNK, by = c("experiment", "plate", "isotype",
                             "analyte"),
            suffix = c("", ".std")) %>%
  mutate(data.norm = log10(data/data.std)+3) %>% 
  mutate(exp_plate = paste(experiment, plate, sep = ".")) %>% 
  unique()

rm(dataInputNK, dataInputQCNK, dataInputQCIgGNK, standardNK)

# Quantify IgG for NK
multiInputListNK <-
  dataInputNormNK %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgG") %>% 
  group_split(experiment, plate_assay)

names(multiInputListNK) <- unique(dataInputNormNK$exp_plate[dataInputNormNK$assaytype == "Multiplex"])

outlistNK <- list()
for(i in 1:length(multiInputListNK)){
  outlistNK[[i]] <- quantifyMultiplexNK(ecdata = ecdata, multiInputList = multiInputListNK)
}

names(outlistNK) <- names(multiInputListNK)

liststandardsNK <- lapply(outlistNK, `[[`, 1)
listdataNK <- lapply(outlistNK, `[[`, 2)

# Quantify IgM for NK
multiInputListNKIgM <-
  dataInputNormNK %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgM") %>% 
  group_split(experiment, plate_assay)

names(multiInputListNKIgM) <- unique(dataInputNormNK$exp_plate[dataInputNormNK$assaytype == "Multiplex"])

outlistNKIgM <- list()
for(i in 1:length(multiInputListNKIgM)){
  outlistNKIgM[[i]] <- quantifyMultiplexReFitNK(ecdata = ecdata, multiInputList = multiInputListNKIgM,
                                            allanalytesIn = outlistNK[[i]]$standards, i)
}

names(outlistNKIgM) <- names(multiInputListNKIgM)

liststandardsNKIgM <- lapply(outlistNKIgM, `[[`, 1)
listdataNKIgM <- lapply(outlistNKIgM, `[[`, 2)

# Quantify IgA for NK
multiInputListNKIgA <-
  dataInputNormNK %>% 
  rename(plate_assay = plate, dilution_assay = dilution) %>% 
  filter(assaytype == "Multiplex" & isotype == "IgA") %>% 
  group_split(experiment, plate_assay)

names(multiInputListNKIgA) <- unique(dataInputNormNK$exp_plate[dataInputNormNK$assaytype == "Multiplex"])

outlistNKIgA <- list()
for(i in 1:length(multiInputListNKIgA)){
  outlistNKIgA[[i]] <- quantifyMultiplexReFitNK(ecdata = ecdata, multiInputList = multiInputListNKIgA,
                                                allanalytesIn = outlistNK[[i]]$standards, i)
}

names(outlistNKIgA) <- names(multiInputListNKIgA)

liststandardsNKIgA <- lapply(outlistNKIgA, `[[`, 1)
listdataNKIgA <- lapply(outlistNKIgA, `[[`, 2)

dataInputQuantNK <-
  rbind(do.call(rbind,listdataNK),
        do.call(rbind,listdataNKIgA),
        do.call(rbind,listdataNKIgM)) %>% 
  mutate(data = if_else(assaytype == "ELISA", log10(data), data)) %>% 
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
  mutate(quant_conc_fix = case_when(is.na(quant_conc) & grepl("above", remarks) ~ quant_single_dil, #max(quant_conc, na.rm = TRUE), # Deviation from other fitting due to singel diution
                                    is.na(quant_conc) & grepl("below", remarks) ~ min(quant_conc, na.rm = TRUE),
                                    TRUE ~ quant_conc),
         quant_conc_fix = if_else(quant_conc_fix < 1, 1, quant_conc_fix)) %>% 
  mutate(panel_corr = case_when(panel == "CPXV_MVA" & grepl("S-", sampleID_metadata) ~ "MVA",
                                panel == "CPXV_MVA" & grepl("P-", sampleID_metadata) ~ "CPXV",
                                panel == "MPXV" ~ "MPXV",
                                panel == "NK" ~ "NK")) %>%
  ungroup() %>% 
  filter(!(sampleID_metadata %in% excludeSamples))

save(dataInputQuantNK, file = "output/dataInputQuantNK.Rdata")


export(dataInputQuant, "output/dataInputQuant.xlsx")

