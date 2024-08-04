## 
# Generate Dataframe containing both MPox data and SPox data as datainput for ML
# Daniel Stern
# RKI
# 2023-09-01
# Version 1
##


# Load libraries and clean environment
library("rio")
library("tidyverse")
library("lubridate")
library("ggpubr")
library("ggthemes")

rm(list = ls(all.names = TRUE))

# Load data
load("../14 Vergleich Spox/data/dataInQuantSpox.R") # dataIn -> Spox Dataframe IgG
load("../14 Vergleich Spox/data/dataInSpox.R") # dataIn -> Spox Dataframe IgG
load("../9 Correlation NT/input/lassoData.Rdata")
load("../4 Method Comparison/output/dataInputQuantCat.Rdata")
load("../4 Method Comparison/output/dataInputQuantCatNK.Rdata")
load("../6 LDA/output/dataInputIgG.Rdata")

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
         MPX_vac_final = if_else(serostatus.meta == 0, "No", "Yes"),
         serostatus.meta = if_else(is.na(serostatus.meta), 0, serostatus.meta),
         serostatus_factor = if_else(is.na(serostatus_factor), ("Pre"), as.character(serostatus_factor)))


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
                                  panel == "MVA" & is.na(serostatus.meta) ~ "Pre",
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


dataInput %>% 
  select(sampleID_metadata, panel) %>% 
  unique() %>% 
  select(panel) %>% 
  ftable()




export(dataInput, file = "output/dataInputAll.csv")
export(dataIn, file = "output/dataPredictionSPox.csv")



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
