##
# Compare duplicate measurements
# Daniel Stern
# RKI
# Version 1.0
# Last modified: 2023-09-15
##

rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggthemes)
library(caret)

##
# Load data

dataPrediction <- import("../15 Generate Dataframe ML/output/dataPredictionSPox.csv")
dataInputAll <- import("../15 Generate Dataframe ML/output/dataInputAll.csv")

quantifiedRep <- import("input/Quantifizierter Input SPox Doppelmessungen.xlsx")
predictedRep <- import("input/Zusammenfassung Ergebnisse SPox_Doppelmessungen.xlsx")

dataInMeta <- import("../14 Vergleich Spox/MPOX_clinic2.xlsx") %>% 
  dplyr::rename(sampleID_meta = `LAB ID`)

predictedRepFinal <-
  predictedRep %>% 
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


dataPredFiltered <- 
  dataPrediction %>% 
  filter(sampleID_meta %in% unique(predictedRepFinal$sampleID_meta)) %>% 
  select(sampleID_meta, Final_combined) %>% 
  mutate(Final_combined = if_else(Final_combined == "", NA_character_, Final_combined))


conf <-
predictedRepFinal %>% 
  filter(POX_VACC_STATUS_1 != "Yes") %>% 
  left_join(dataPredFiltered, by = "sampleID_meta", suffix = c("", ".org")) %>% 
  select(sampleID_meta,Final_combined, Final_combined.org) %>% 
  mutate(difference = case_when(Final_combined.org == "Pre" &
                                  Final_combined == "MPXV" ~ "Different",
                                Final_combined.org == "Pre" &
                                  Final_combined == "seronegative" ~ "Different",
                                Final_combined.org == "Pre" &
                                  is.na(Final_combined) ~ "Different",
                                Final_combined.org == "seronegative" &
                                  is.na(Final_combined) ~ "Different",
                                Final_combined.org == "MPXV" &
                                  Final_combined == "MPXV" ~ "Same",
                                Final_combined.org == "MVA" &
                                  Final_combined == "MVA" ~ "Same",
                                Final_combined.org == "Pre" &
                                  Final_combined == "Pre" ~ "Same",
                                Final_combined.org == "seronegative" &
                                  Final_combined == "seronegative" ~ "Same",
                                is.na(Final_combined.org) &
                                  is.na(Final_combined) ~ "Same"))


caret::confusionMatrix(data = as.factor(conf$Final_combined),
                       reference = as.factor(conf$Final_combined.org))
levels(reference)
levels(data)
missclass <-
predictedRepFinal %>% 
  left_join(dataPredFiltered, by = "sampleID_meta", suffix = c("", ".org")) %>% 
  filter(Final_combined == "MPXV" & Final_combined.org == "Pre")



predictedRepFinal %>% 
 # filter(POX_VACC_STATUS_1 != "Yes") %>% 
#  filter(MPX_vac_final != "Yes") %>% 
  left_join(dataPredFiltered, by = "sampleID_meta", suffix = c("", ".org")) %>% 
  select(Final_combined, MPX_diagnose_final) %>% 
  ftable()


##
# Compare data
dataInputAllFiltered <-
  dataInputAll %>% 
  filter(sampleID_metadata %in% unique(predictedRepFinal$sampleID_meta))

quantifiedRepSelect <-
  quantifiedRep %>% 
  dplyr::select(sampleID_meta, analyte, isotype, data, dataIn)



dataInputAllFiltered %>% 
  filter(isotype == "IgG") %>% 
  left_join(quantifiedRepSelect, by = c("sampleID_metadata" = "sampleID_meta",
                                        "analyte", "isotype"), suffix = c("", ".rep")) %>% 
  left_join(conf, by = c("sampleID_metadata" = "sampleID_meta")) %>% 
  ggplot(mapping = aes(x = dataIn, y = log10(dataIn.rep))) +
  geom_smooth(method = "lm") +
  scale_color_colorblind() +
  geom_point(aes(color = difference, shape = MPX_diagnose_final), size = 2, alpha = 0.5) +
  facet_wrap("analyte")

