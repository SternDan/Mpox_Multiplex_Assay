##
# Analyse and discern samples with positive Mpox diagnosis 
# by the ML-based algorithm
# Information by Uli
#  ich hätte noch einmal eine Bitte bezüglich der serologischen Testergebnisse.
# 
#  Ich hab die PINs für folgende sample in eine Excel-Datei kopiert:
#   
#  1) Proben, für die keine mpox-Diagnose und keine MVA-Impfung angegeben 
# worden war, die aber Orthopocken-Antikörper positiv waren. In zwei weiteren
# Spalten ist entweder das konkrete Alter in Jahren oder die Altersgruppe
# eingetragen (1= 18-29J; 2= 30-39J.; 3= 40-49J; 4= 50-59J.; 5 = 60+ J.)
# 
# Für Altersgruppe 1 ist eine Pockenimpfung in der Kindheit sehr
# unwahrscheinlich, für Altersgruppe 2 eher unwahrscheinlich.
# Könntest Du Dir die serologischen Ergebnisse für diese Proben nochmal ansehen
# und angeben, ob z.B. IgM positiv war, und ob die Antikörpertiter eher für eine
# lange zurückliegende Impfung oder für ein frischeres Geschehen sprechen.
# 
# >>> samples coded in variable dataDeltaPos 


# 2 a) Proben, für die eine mpox-Diagnose im Fragebogen angegeben wurde, die
# aber nicht durch die Algorithmen eindeutig bestätigt werden konnte. Gibt es
# hier ähnliche Muster wie für auffällige Proben aus 1) und aus 2b)?
#  2 b) Proben, für die eine mpox-Symptomatik angegeben wurde, aber keine
# mpox-Diagnose. Die meisten dieser Teilnehmer hatten eine MVA-Impfung erhalten
# und wurden in den Algorithmen auch als MVA geimpft bewertet. Für fünf Proben
# wurde keine MVA-Impfung angegeben, drei davon waren Orthopocken-Antikörper
# positiv, davon zwei eigentlich zu jung für eine Pockenimpfung in der Kindheit.
# Gibt es für einen Teil dieser Proben Hinweise auf ein kürzer zurückliegendes 
# Infektionsereignis? 
#
# >>> samples coded in variable dataMpoxFN (False Negatives)
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



##
# Plot information 
plotMPXVAcuteAgeGroup <-
  dataInputAcuteMPXV %>% 
  select(sampleID_metadata, year_birth, POX_VACC_STATUS_1) %>% 
  filter(!is.na(year_birth)) %>% 
  filter(POX_VACC_STATUS_1 %in% c("Ja", "Nein")) %>% 
  unique() %>% 
  mutate(
    Age_group = case_when(year_birth >= 1994 ~ "< 30",
                          year_birth >= 1984 ~ "< 40",
                          year_birth >= 1974 ~ "< 50",
                          year_birth >= 1964 ~ "< 60",
                          year_birth < 1964 ~ "60+"),
    Age_group = factor(Age_group, levels = (c("< 30", 
                                              "< 40", 
                                              "< 50", 
                                              "< 60", 
                                              "60+")),
                       ordered = TRUE)) %>% 
  group_by(Age_group) %>% 
  summarise(Prop_Vaccinated = mean(POX_VACC_STATUS_1 == "Ja"),
            N = n(),
            Std_Error = sqrt(Prop_Vaccinated*(1-Prop_Vaccinated)/N)) %>% 
  mutate( Lower = Prop_Vaccinated - 1.96*Std_Error, 
          Upper = Prop_Vaccinated + 1.96*Std_Error) %>% 
  ggplot(aes(x = Age_group, y = Prop_Vaccinated)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Proportion Self-reported Pox Vaccinated",
                     limits = c(0,1)) +
  theme_pubr()

plotSpoxAgeGroups <-
  dataAgeGroupsSpox %>% 
  select(sampleID_metadata, Age_group_all, POX_VACC_STATUS_1) %>% 
  filter(!is.na(Age_group_all)) %>% 
  filter(POX_VACC_STATUS_1 %in% c("Yes", "No")) %>% 
  unique() %>% 
  group_by(Age_group_all) %>% 
  summarise(Prop_Vaccinated = mean(POX_VACC_STATUS_1 == "Yes"),
            N = n(),
            Std_Error = sqrt(Prop_Vaccinated*(1-Prop_Vaccinated)/N)) %>% 
  mutate( Lower = Prop_Vaccinated - 1.96*Std_Error, 
          Upper = Prop_Vaccinated + 1.96*Std_Error) %>% 
  ggplot(aes(x = Age_group_all, y = Prop_Vaccinated)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Proportion Self-reported Pox Vaccinated",
                     limits = c(0,1)) +
  theme_pubr()

plotSpoxMVAAgeGroups <-
  dataAgeGroupsSpox %>% 
  select(sampleID_metadata, Age_group_all, MPX_VACC_STATUS_1) %>% 
  filter(!is.na(Age_group_all)) %>% 
  filter(MPX_VACC_STATUS_1 %in% c("Yes", "No")) %>% 
  unique() %>% 
  group_by(Age_group_all) %>% 
  summarise(Prop_Vaccinated = mean(MPX_VACC_STATUS_1 == "Yes"),
            N = n(),
            Std_Error = sqrt(Prop_Vaccinated*(1-Prop_Vaccinated)/N)) %>% 
  mutate( Lower = Prop_Vaccinated - 1.96*Std_Error, 
          Upper = Prop_Vaccinated + 1.96*Std_Error) %>% 
  ggplot(aes(x = Age_group_all, y = Prop_Vaccinated)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Proportion Self-reported MVA Vaccinated",
                     limits = c(0,1)) +
  theme_pubr()

plotAgeGroupsVacc <-
  ggarrange(plotMPXVAcuteAgeGroup, plotSpoxAgeGroups, plotSpoxMVAAgeGroups, align = "hv", labels = 
              c("a", "b", "c"), ncol = 3)

ggsave(file = "output/plotAgeVacc.png", plotAgeCroupsVacc, width = 10, height = 4,
       dpi = 600)











# Load data to check
dataDeltaPos <- import("input/Mpox-samples_serocheck.xlsx", range = "A2:C107", sheet = "Pre_") %>% 
  mutate(Age_group = case_when(Age < 30 ~ 1,
                               Age < 40 ~ 2,
                               Age < 50 ~ 3, 
                               Age < 60 ~ 4, 
                               (!is.na(Age)) ~ 5,
                               TRUE ~ Age_group),
         Age_group = factor(Age_group, levels = c(1,2,3,4,5),
                            labels = c("< 30", 
                                       "< 40", 
                                       "< 50", 
                                       "< 60", 
                                       "60+"),
                            ordered = TRUE))

dataMpoxFN <- import("input/Mpox-samples_serocheck.xlsx", range = "A2:A21", sheet = "MPXV") %>% 
  rename(PIN = `PIN MPXV unconfirmed`)
dataMpoxSymptFN <- import("input/Mpox-samples_serocheck.xlsx", range = "A24:E44", sheet = "MPXV") %>% 
  mutate(Age_group = factor(AGE, levels = c(1,2,3,4,5),
                            labels = c("< 30", 
                                       "< 40", 
                                       "< 50", 
                                       "< 60", 
                                       "60+"),
                            ordered = TRUE))

##
# Analyse Delta Positive Data -> How many are clear positive?
dataPosCombined <- 
  dataCombined %>% 
  filter(patientID %in% c(unique(dataDeltaPos$PIN))) %>% 
  left_join(dataDeltaPos, by = c("patientID" = "PIN"))

plotSpoxPoxIgM <-
  dataPosCombined %>% 
  mutate(sumNA = (is.na(pred_AllIgGIgM) + is.na(pred_AllIgG) + is.na(pred_EpiIgG) +is.na(pred_EpiIgGIgM))) %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotSpoxPoxIgG <-
  dataPosCombined %>% 
  mutate(sumNA = (is.na(pred_AllIgGIgM) + is.na(pred_AllIgG) + is.na(pred_EpiIgG) +is.na(pred_EpiIgGIgM))) %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotPrePos <- 
  ggarrange(plotSpoxPoxIgG, plotSpoxPoxIgM, align = "hv",
            common.legend = TRUE, ncol = 2)
ggsave("output/plotPrePos.png", plotPrePos, width = 12, height = 8, dpi = 600) 


##
# Plot age groups on all datasets
# Plot on Pre
plotSpoxPreIgM <-
  dataCombined %>% 
  filter(panel_detail == "Pre") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        #   legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotSpoxPreIgG <-
  dataCombined %>% 
  filter(panel_detail == "Pre") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotPre <- 
  ggarrange(plotSpoxPreIgG, plotSpoxPreIgM, align = "hv",
            common.legend = TRUE, ncol = 2)
ggsave("output/plotPre.png", plotPre, width = 12, height = 8, dpi = 600) 

# Plot on MVA
plotSpoxMVAIgM <-
  dataCombined %>% 
  filter(panel_detail == "MVA") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        #   legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotSpoxMVAIgG <-
  dataCombined %>% 
  filter(panel_detail == "MVA") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotMVA <- 
  ggarrange(plotSpoxMVAIgG, plotSpoxMVAIgM, align = "hv",
            common.legend = TRUE, ncol = 2)
ggsave("output/plotMVA.png", plotMVA, width = 12, height = 8, dpi = 600) 


# Plot on Mpox
plotSpoxMPXVIgM <-
  dataCombined %>% 
  filter(panel_detail == "MPXV") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        #   legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotSpoxMPXVIgG <-
  dataCombined %>% 
  filter(panel_detail == "MPXV") %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group_all) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group_all, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotMPXV <- 
  ggarrange(plotSpoxMPXVIgG, plotSpoxMPXVIgM, align = "hv",
            common.legend = TRUE, ncol = 2)
ggsave("output/plotMPXV.png", plotMPXV, width = 12, height = 8, dpi = 600) 


##
# Plot NK Data
plotNKIgM <-
  dataInputNKcombined %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        #   legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotNKIgG <-
  dataInputNKcombined %>% 
  select(sampleID_metadata, isotype, analyte, serostatus_cat, Age_group) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, fill = serostatus_cat))+
  geom_bar(position = "fill") +
  facet_wrap("analyte") +
  theme_pubr()+
  scale_fill_manual(name = "Serostatus", values = c("white", "grey90", colorblind_pal()(8)[c(3,4)]))+ 
  theme(strip.background = element_blank(),
        #   legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotNK <- 
  ggarrange(plotNKIgG, plotNKIgM, align = "hv",
            common.legend = TRUE, ncol = 2)
ggsave("output/plotNK.png", plotNK, width = 12, height = 8, dpi = 600) 

##
# Arrange IgG plots only
plotIgG <-
  ggarrange(plotSpoxPreIgG,  plotSpoxPoxIgG, plotSpoxMVAIgG, plotNKIgG, plotSpoxMPXVIgG,
            ncol = 2, nrow = 3,
            align = "hv", common.legend = TRUE,
            labels = c("Pre", "Pre Pos", "MVA", "NK", "MPXV"))

ggsave("output/plotIgG.png", plotIgG, width = 12, height = 20, dpi = 600) 

##
# Plot Quantified data -> only for selected antigens
# Select the following antigens based on the bar plot
# - Good antigens: D8L, E8L, 
# - Good antigens EEV: A33R, A35R
# - Good antigens EEV 2: B5R, B6R
# - High spec Antigen: ATI-N



##
# Plot different panels with dataIn
plotDataInPanels <- 
  dataPanel %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte %in% c("D8L", "E8L",
                        "A33R", "A35R",
                        "B5R", "B6R", 
                        "ATI-N")) %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, y = dataIn, fill = analyte))+
  geom_boxplot() +
  facet_grid(analyte ~ panel_detail, scales = "free_y") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotDataInPanels.png", plotDataInPanels, width = 8, height = 10,
       dpi = 600)


##
# Bearbeitung der Frage 2 a: Proben, die MPox positive sind, aber nicht eindeutig
# bestätigt werden konnten. Gibt es ein ähnliches Muster wie bei Frage 1 oder 2b?
dataPanelMPXVFN <- 
  dataPanel %>% 
  filter(PIN %in% dataMpoxFN$PIN) %>% 
  mutate(panel_detail = "MPXV False Neg")

# Bearbeitung der Frage 2 b: 
dataPanelMPXVFP <- 
  dataPanel %>% 
  filter(PIN %in% dataMpoxSymptFN$PIN) %>% 
  mutate(panel_detail = "MPXV False Pos")


plotDataInPanelsGood <- 
  dataPanel %>% 
  rbind(dataPanelMPXVFN, dataPanelMPXVFP) %>% 
  filter(panel_detail %in% c("NC", "Pre", "Pre Pos", "MVA", "MPXV False Pos", "MPXV False Neg", "MPXV")) %>% 
  filter(analyte %in% c("D8L", "E8L",
                        "A33R", "A35R",
                        "B5R", "B6R", 
                        "ATI-N")) %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE),
         panel_detail = factor(panel_detail, levels = c("NC",
                                                        "Pre", 
                                                        "Pre Pos", 
                                                        "MVA", 
                                                        "MPXV False Pos",
                                                        "MPXV False Neg", 
                                                        "MPXV"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, y = dataIn, fill = analyte))+
  geom_boxplot() +
  facet_grid(analyte ~ panel_detail, scales = "free_y") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotDataInPanelsGood.png", plotDataInPanelsGood, width = 10, height = 10,
       dpi = 600)


plotDataInPanelsBad <- 
  dataPanel %>% 
  rbind(dataPanelMPXVFN, dataPanelMPXVFP) %>% 
  filter(panel_detail %in% c("NC", "Pre", "Pre Pos", "MVA", "MPXV False Pos", "MPXV False Neg", "MPXV")) %>% 
  filter(!(analyte %in% c("D8L", "E8L",
                          "A33R", "A35R",
                          "B5R", "B6R", 
                          "ATI-N"))) %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgG") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE),
         panel_detail = factor(panel_detail, levels = c("NC",
                                                        "Pre", 
                                                        "Pre Pos", 
                                                        "MVA", 
                                                        "MPXV False Pos",
                                                        "MPXV False Neg", 
                                                        "MPXV"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, y = dataIn, fill = analyte))+
  geom_boxplot() +
  facet_grid(analyte ~ panel_detail, scales = "free_y") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(1:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotDataInPanelsBad.png", plotDataInPanelsBad, width = 10, height = 10,
       dpi = 600)

## Plot IgM
plotDataInPanelsGoodIgM <- 
  dataPanel %>% 
  rbind(dataPanelMPXVFN, dataPanelMPXVFP) %>% 
  filter(panel_detail %in% c("NC", "Pre", "Pre Pos", "MVA", "MPXV False Pos", "MPXV False Neg", "MPXV")) %>% 
  filter(analyte %in% c("D8L", "E8L",
                        "A33R", "A35R",
                        "B5R", "B6R", 
                        "ATI-N")) %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE),
         panel_detail = factor(panel_detail, levels = c("NC",
                                                        "Pre", 
                                                        "Pre Pos", 
                                                        "MVA", 
                                                        "MPXV False Pos",
                                                        "MPXV False Neg", 
                                                        "MPXV"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, y = dataIn, fill = analyte))+
  geom_boxplot() +
  facet_grid(analyte ~ panel_detail, scales = "free_y") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotDataInPanelsGoodIgM.png", plotDataInPanelsGoodIgM, width = 10, height = 10,
       dpi = 600)


plotDataInPanelsBadIgM <- 
  dataPanel %>% 
  rbind(dataPanelMPXVFN, dataPanelMPXVFP) %>% 
  filter(panel_detail %in% c("NC", "Pre", "Pre Pos", "MVA", "MPXV False Pos", "MPXV False Neg", "MPXV")) %>% 
  filter(!(analyte %in% c("D8L", "E8L",
                          "A33R", "A35R",
                          "B5R", "B6R", 
                          "ATI-N"))) %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgM") %>% 
  mutate(serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"), 
                                 ordered = TRUE),
         panel_detail = factor(panel_detail, levels = c("NC",
                                                        "Pre", 
                                                        "Pre Pos", 
                                                        "MVA", 
                                                        "MPXV False Pos",
                                                        "MPXV False Neg", 
                                                        "MPXV"), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group, y = dataIn, fill = analyte))+
  geom_boxplot() +
  facet_grid(analyte ~ panel_detail, scales = "free_y") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(1:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotDataInPanelsBadIgM.png", plotDataInPanelsBadIgM, width = 10, height = 10,
       dpi = 600)

##
# Analyse Mpox diagnosed samples
# 281848 – 64 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 281854 – 42 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 281859 – 33 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 281866 – 58 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 291874 – 77 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 291878 -- 37 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 291886 – 56 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 291909 – 62 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 362153 – 54 Jahre alt, eine Mpox-Diagnose wurde dokumentiert.
# 362154 -- 60 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 362168 – 56 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 362186 – 55 Jahre alt, eine Mpox-Diagnose wurde nicht dokumentiert.
# 141302   ist jetzt aktuell 65 Jahre alt und hat anamnestisch keine MPX Infektion gehabt
# 141304   aktuell 32 Jahre, keine MPX Diagnose in der Anamnese (könnte aber passen, hatte unklare klinische Symptome) 
# 382245  aktuell 52 Jahre , keine MPX Diagnose (klinisch auch nicht auszuschließen)
# 382262  ist aktuell 33 Jahre, cave : hier ist wohl was schief gelaufen bei uns (der ist gar kein MSM und es gibt natürlich auch keine MPX Diag Anamnese und es gibt auch kein dokumentiertes Einverständnis, da muss irgendwas verwechselt worden sein, vielleicht besser komplett rausnehmen aus der Studie und Auswertung
# 382264 aktuell 35 Jahre, keine MPX Anamnese
# 392278 aktuell 28 Jahre,  der hatte gesichert  MPX Diagnose  !!!
# 392288 , aktuell 34 Jahre, der hatte auch gesichert MPX Diagn ( PCR pos.)  !
# 392297, aktuell 23 Jahr, keine MPX Dx Anamnese

# Generate dataframe with information on samples
dataSPoxInfo <- data.frame(PIN = 
                             c(281848, 281854, 281859, 281866, 291874, 291878, 291886,
                               291909, 362153, 362154, 362168, 362186, 141302, 141304,
                               382245, 382262, 382264, 392278, 392288, 392297),
                           AGE = c(64, 42, 33, 58, 77, 37, 56, 62, 54, 60, 56, 55,
                                   65, 32, 52, 33, 35, 28, 34, 23),
                           MPOX_DIAGNOSED = c("No", "No", "No", "No", "No", "No", "No", "No",
                                              "Yes", "No", "No", "No",
                                              "No", "Possible", "Possible", "False", "No", "Yes", "Yes", "No"))


plotSpoxUnknown <- 
  dataCombined %>% 
  left_join(dataSPoxInfo, by = "PIN") %>% 
  filter(PIN %in% dataSPoxInfo$PIN | panel_detail %in% c("Pre", "MPXV", "MVA")) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "SPox", "MPXV"),
                               ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = panel_detail, y = dataIn))+
  geom_boxplot(aes(fill = panel_detail)) +
  geom_point(aes(color = as.factor(MPOX_DIAGNOSED))) +
  facet_wrap("analyte") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotSpoxUnknown.png", plotSpoxUnknown, width = 10, height = 10,
       dpi = 600)



plotSpoxUnknownIgM <- 
  dataCombined %>% 
  left_join(dataSPoxInfo, by = "PIN") %>% 
  filter(PIN %in% dataSPoxInfo$PIN | panel_detail %in% c("Pre", "MPXV", "MVA")) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgM") %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "SPox", "MPXV"),
                               ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = panel_detail, y = dataIn))+
  geom_boxplot(aes(fill = panel_detail)) +
  geom_point(aes(color = as.factor(MPOX_DIAGNOSED))) +
  facet_wrap("analyte") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("output/plotSpoxUnknownIgM.png", plotSpoxUnknownIgM, width = 10, height = 10,
       dpi = 600)



plotSpoxUnknownYoung <- 
  dataCombined %>% 
  left_join(dataSPoxInfo, by = "PIN") %>% 
  filter(MPOX_DIAGNOSED != "False" | is.na(AGE.y)) %>% 
  filter(AGE.y < 50 | is.na(AGE.y)) %>% 
  filter(PIN %in% dataSPoxInfo$PIN | panel_detail %in% c("Pre", "MPXV", "MVA")) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "SPox", "MPXV"),
                               ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = panel_detail, y = dataIn))+
  geom_boxplot(aes(fill = panel_detail)) +
  geom_point(aes(color = as.factor(MPOX_DIAGNOSED))) +
  facet_wrap("analyte") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotSpoxUnknownOld <- 
  dataCombined %>% 
  left_join(dataSPoxInfo, by = "PIN") %>% 
  filter(MPOX_DIAGNOSED != "False"| is.na(AGE.y)) %>% 
  filter(AGE.y >= 50 | is.na(AGE.y)) %>% 
  filter(PIN %in% dataSPoxInfo$PIN | panel_detail %in% c("Pre", "MPXV", "MVA")) %>% 
  filter(analyte != "VACV") %>% 
  filter(isotype == "IgG") %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "SPox", "MPXV"),
                               ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = panel_detail, y = dataIn))+
  geom_boxplot(aes(fill = panel_detail)) +
  geom_point(aes(color = as.factor(MPOX_DIAGNOSED))) +
  facet_wrap("analyte") +
  theme_bw()+
  scale_fill_manual(name = "Analyte", values =  colorblind_pal()(8)[c(2:8)])+ 
  theme(strip.background = element_blank(),
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



ggsave("output/plotSpoxUnknown.png", plotSpoxUnknown, width = 10, height = 10,
       dpi = 600)


ggsave("output/plotSpoxUnknownYoung.png", plotSpoxUnknownYoung, width = 10, height = 10,
       dpi = 600)
ggsave("output/plotSpoxUnknownOld.png", plotSpoxUnknownOld, width = 10, height = 10,
       dpi = 600)
