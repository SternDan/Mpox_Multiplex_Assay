## 
# Generate Dataframe containing both MPox data and SPox data as datainput for ML
# Daniel Stern
# RKI
# 2023-09-14
# Version 2.1
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
#load("../9 Correlation NT/input/lassoData.Rdata")

load("../4 Method Comparison/output/dataInputQuantCat.Rdata")
load("../4 Method Comparison/output/dataInputQuantCatNK.Rdata")
load("../4 Method Comparison/output/threshold.Rdata")
cutoffDelta <- log10(threshold$median[threshold$assay == "IgG" & 
                                        threshold$antigene == "Delta" & 
                                        threshold$assaytype == "Multiplex"])

##
# Import data that have been matched with metadata
load("../7 Patient panel merge/output/dataInputMPXVmeta.Rdata") # dataInputMPXVmeta
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

dataInputRaw %>% 
  dplyr::select(sampleID_metadata, panel_in) %>% 
  unique() %>% 
  dplyr::select(panel_in) %>% 
  ftable()

dataInputRawMPXV <-
  dataInputRaw %>% 
  dplyr::select(sampleID_metadata) %>% 
  unique() %>% 
  filter(grepl("P-22", sampleID_metadata)) %>% 
  pull()

metaDataContained <-
unique(metadata_MPXV_patients$sampleID_metadata)[(unique(metadata_MPXV_patients$sampleID_metadata) %in% dataInputRawMPXV)]
dataInputRawMPXV[!(dataInputRawMPXV %in% metaDataContained)]

load("../2 Import metadata/output/metadata_MVA_time.Rdata")
load("../2 Import metadata/output/metadata_MPXV_patients.Rdata")

#load("../6 LDA/output/dataInputIgG.Rdata") # 

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
  dataInputMVAexport %>% 
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

dataInputMPXV <- 
  dataInputMPXVexport %>% 
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

dataInput <-
  rbind(dataInputMPXV, dataInputMVA, dataInputSPox, dataInputNK, dataInputCPXV) %>% 
  mutate(panel_detail = case_when(panel == "MVA" & serostatus.meta == 0 ~ "Pre",
                                panel == "SPox" & MPX_diagnose_final == "Yes" ~ "MPXV",
                                panel == "SPox" & MPX_diagnose_final == "No" &
                                  MPX_vac_final == "Yes" ~ "MVA",
                                panel == "SPox" & MPX_diagnose_final == "No" &
                                  MPX_vac_final == "No" ~ "Pre",
                                panel == "Pre" ~ "Pre",
                                panel == "MVA" & serostatus.meta != 0 ~ "MVA",
                                panel == "CPXV" ~ "CPXV",
                                panel == "MPXV" ~ "MPXV")
         )


length(unique(dataInput$sampleID_metadata))
unique(dataInput$panel_detail)

dataInputIgGMPOX <-
  dataInputIgG %>% 
  filter(panel_in == "MPXV")

unique(dataInputMPXV$sampleID_metadata) %in% unique(dataInputIgGMPOX$sampleID_metadata)


export(dataInput, file = "output/dataInputAll.csv")
export(dataIn, file = "output/dataPredictionSPox.csv")

names(dataInput)
summary(dataInput)
unique(dataInput$panel)

na_input <-
dataInput %>% 
  filter(is.na(panel_detail))

dataInputComplete <-
  dataInput %>% 
  filter(!is.na(panel_detail)) %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  filter(analyte != "VACV")

# Plot comparison for MPox Positive Sera used for validation and for testing
plotComparisonMpoxSpoxPositive <-
dataInput %>% 
  filter(panel %in% c("MPXV", "SPox")) %>%
  filter(MPX_diagnose_final == "Yes") %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  ggplot(mapping = aes(x = panel, y = dataIn, fill = panel)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75
  ), alpha = 0.1) +
  facet_grid(isotype ~ analyte) +
  theme_bw()+
  stat_compare_means(comparisons = list(c("MPXV", "SPox")),
                     label = "p.signif") +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("Multiplex results") 

plotComparisonMpoxSpoxNegative <-
dataInput %>% 
  filter(panel %in% c("MPXV", "SPox")) %>%
  filter((MPX_diagnose_final == "No" & panel == "SPox") |
           panel == "MPXV") %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  ggplot(mapping = aes(x = panel, y = dataIn, fill = panel)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75
  ), alpha = 0.1) +
  facet_grid(isotype ~ analyte) +
  theme_bw()+
  stat_compare_means(comparisons = list(c("MPXV", "SPox")),
                     label = "p.signif") +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("Multiplex results") 


ggsave("output/plotComparisonMpoxSpoxPositive.png", plotComparisonMpoxSpoxPositive,
       width = 10, height = 7, dpi = 600)

ggsave("output/plotComparisonMpoxSpoxNegative.png", plotComparisonMpoxSpoxNegative,
       width = 10, height = 7, dpi = 600)


igmPositive <- 
  dataInput %>% 
  filter(isotype == "IgM" & analyte %in% "Delta" & panel == "SPox") %>% 
  filter(serostatus == 1)

unique(igmPositive$sampleID_metadata)

igmSample <- 
  dataInput %>% 
  filter(sampleID_metadata %in% c("SPox-0133", "SPox-0141", "SPox-0428", "SPox-0532","SPox-0542", "SPox-0626", "SPox-0693",
                              "SPox-0731", "SPox-0721", "SPox-0796", "SPox-0782", "SPox-0847"))

# Plot duplicate measuremens
dataInput %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = data, y = dataIn, color = batch)) +
  geom_point() +
  facet_wrap("analyte")
