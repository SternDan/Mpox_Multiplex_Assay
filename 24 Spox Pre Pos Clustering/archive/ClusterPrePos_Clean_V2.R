####
# Analyse OPV positive samples from SPox study
# Daniel Stern RKI Uli Markus RKI
# 2024/02/06
####

####
# Aim of the analysis
# Within the Spox study, several sera, that are classified as "Pre",
# meaning that neither MVA vaccination nor Mpox infections have occured,
# show detectable levels of anti-OPXV antibodies. 
# As those sera have been trained by the ML-algorithms as "Pre" or they were
# inconclusive, in this analysis we want to test, if those samples would cluster
# more closely with Mpox infected samples, or with MVA immunized samples. 
# As an additional class of samples, NC samples will also be included
# In a first attempt, we will try to cluster the samples stratified by age
# In a second attemt, we will try to test if OPXV positive samples are also
# clearly positive for the other antigens, tested, only focusing
# In a third attempt, we will try to use ATI-N as a marker for infection
# in MVA immunized persons. To this aim, we will need to load the MVA 
# immunized data in the panel of younger persons.
####

library(rio)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(ggalt)
library(factoextra)
library(cluster)
library(pROC)
library(data.table)
library(caret)

rm(list = ls(all.names = TRUE))

# Source functions
source("determineCutoffFunction.R", encoding = "UTF-8")
source("rocFunction.R", encoding = "UTF-8")

## Load datainput generated for clustering
load("input/dataClustering.Rdata")

## Unify the format for the Age_groups varibale
dataClustering_pre <- 
  dataClustering %>% 
  mutate(Age_group = case_when(Age_group == "18-29" ~ "< 30",
                               Age_group == "30-39" ~ "< 40",
                               Age_group == "40-49" ~ "< 50",
                               Age_group == "50-59" ~ "< 60",
                               TRUE ~ Age_group),
         Age_group = factor(Age_group, levels = c("< 30", 
                                                  "< 40", 
                                                  "< 50",
                                                  "< 60", 
                                                  "60+"),
                            ordered = TRUE))




predRFwide <-
  dataClustering %>% 
  select(sampleID_metadata, pred_RFAllIgGIgM) %>% 
  unique() %>% 
  mutate(nm1= str_c('predRF_', rowid(sampleID_metadata))) %>% 
  pivot_wider(names_from = nm1, values_from = pred_RFAllIgGIgM)

predXGboostwide <-
  dataClustering %>% 
  select(sampleID_metadata, predxgboostAllIgGIgM) %>% 
  unique() %>% 
  mutate(nm1= str_c('predxgboost_', rowid(sampleID_metadata))) %>% 
  pivot_wider(names_from = nm1, values_from = predxgboostAllIgGIgM)

predLDAwide <-
  dataClustering %>% 
  select(sampleID_metadata, pred_LDAAllIgGIgM) %>% 
  unique() %>% 
  mutate(nm1= str_c('preLDA_', rowid(sampleID_metadata))) %>% 
  pivot_wider(names_from = nm1, values_from = pred_LDAAllIgGIgM)

dataClustering_post <- 
  dataClustering_pre %>% 
  select(-pred_RFAllIgGIgM, - predxgboostAllIgGIgM, 
         -pred_LDAAllIgGIgM) %>% 
  unique() %>% 
  left_join(predRFwide, by = "sampleID_metadata") %>% 
  left_join(predXGboostwide, by = "sampleID_metadata") %>% 
  left_join(predLDAwide, by = "sampleID_metadata") 




####
## Prepare wide dataframes for further analysis
## SPox dataframe
dataInWideSpox <-
  dataClustering_post %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel %in% "SPox") %>% 
  filter(panel_detail != "SPox") %>% 
  select(-PIN, - data) %>% 
  pivot_wider(names_from = analyte, values_from = c(dataIn, serostatus, serostatus_cat)) %>% 
  mutate(panel_clustering = case_when(panel_detail == "Pre" & serostatus_Delta == 1 ~ "Pre positive",
                                      panel_detail == "Pre" & serostatus_Delta == 0 ~ "Pre negative",
                                      TRUE ~ panel_detail))

## NC dataframe
dataInWideNC <-
  dataClustering_post %>% 
  filter(!is.na(Age_group)) %>% 
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel %in% "Pre_New") %>% 
  filter(panel_detail != "SPox") %>% 
  select(-PIN, - data) %>% 
  pivot_wider(names_from = analyte, values_from = c(dataIn, serostatus, serostatus_cat)) %>% 
  mutate(panel_clustering = case_when(panel_detail == "Pre" & serostatus_Delta == 1 ~ "Pre positive",
                                      panel_detail == "Pre" & serostatus_Delta == 0 ~ "Pre negative",
                                      TRUE ~ panel_detail))


## MVA acute dataframe
dataInWideMVA <-
  dataClustering_post %>% 
  filter(!is.na(Age_group)) %>% 
  filter(!sampleID_metadata %in% c("S-20-004-007-02", "S-20-004-007-01")) %>% # Remove duplicates
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel %in% "MVA") %>% 
  filter(panel_detail != "SPox") %>% 
  select(-PIN, - data) %>% 
  pivot_wider(names_from = analyte, values_from = c(dataIn, serostatus, serostatus_cat)) %>% 
  mutate(panel_clustering = case_when(panel_detail == "Pre" & serostatus_Delta == 1 ~ "Pre positive",
                                      panel_detail == "Pre" & serostatus_Delta == 0 ~ "Pre negative",
                                      TRUE ~ panel_detail))

## MPXV acture dataframe
dataInWideMPXV <-
  dataClustering_post %>% 
  filter(!is.na(Age_group)) %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08")) %>% # Remove duplicates
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel %in% "MPXV") %>% 
  filter(panel_detail != "SPox") %>% 
  select(-PIN, - data) %>% 
  pivot_wider(names_from = analyte, values_from = c(dataIn, serostatus, serostatus_cat)) %>% 
  mutate(panel_clustering = case_when(panel_detail == "Pre" & serostatus_Delta == 1 ~ "Pre positive",
                                      panel_detail == "Pre" & serostatus_Delta == 0 ~ "Pre negative",
                                      TRUE ~ panel_detail))


## Bind all wide dataframes together
dataInWide <-
  rbind(dataInWideSpox, dataInWideNC, dataInWideMPXV, dataInWideMVA) %>% 
  rowwise() %>% 
  mutate(sum_seropositive = sum(serostatus_D8L, serostatus_L1R,
                                serostatus_B5R, serostatus_A33R,
                                serostatus_M1R, `serostatus_ATI-N`,
                                serostatus_B6R, serostatus_E8L, 
                                serostatus_A5L, serostatus_A29L,
                                serostatus_H3L, serostatus_A35R,
                                `serostatus_ATI-C`, serostatus_A27L,
                                serostatus_Delta),
         ratioAG_seropositive = sum_seropositive/15,
         clear_positive = sum(serostatus_cat_D8L == "positive", serostatus_cat_L1R == "positive",
                              serostatus_cat_B5R == "positive", serostatus_cat_A33R == "positive",
                              serostatus_cat_M1R == "positive", `serostatus_cat_ATI-N` == "positive",
                              serostatus_cat_B6R == "positive", serostatus_cat_E8L == "positive", 
                              serostatus_cat_A5L == "positive", serostatus_cat_A29L == "positive",
                              serostatus_cat_H3L == "positive", serostatus_cat_A35R == "positive",
                              `serostatus_cat_ATI-C` == "positive", serostatus_cat_A27L == "positive",
                              serostatus_cat_Delta == "positive"),
         ratioAG_clear_seropositive = clear_positive/15)


dataSeropositive <-
  dataInWide %>% 
  select(sampleID_metadata, ratioAG_seropositive, ratioAG_clear_seropositive) %>% 
  unique()


#####
## Dev: Find duplicates to be removed before the dataframes are widened
# duplicatesMVA <-
#   dataClustering %>% 
#   filter(!is.na(Age_group)) %>% 
#   filter(isotype == "IgG") %>% 
#   filter(analyte != "VACV") %>% 
#   filter(panel %in% "MVA") %>% 
#   filter(panel_detail != "SPox") %>% 
#   select(-PIN, - data)  %>%
#   dplyr::group_by(panel, panel_detail, isotype, sampleID_metadata, Age_group,
#                   POX_VACC_STATUS_ML, POX_VACC_STATUS_SELF, MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
#                   MPX_vac_final, MPX_SYMPTOME, MPX_diagnose_final, real, pred_AllIgGIgM, pred_AllIgG,
#                   pred_EpiIgGIgM, pred_EpiIgG, pred_RFAllIgGIgM, predxgboostAllIgGIgM, pred_LDAAllIgGIgM,
#                   analyte) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L) 
# 
# duplicatesMPXV <-
#   dataClustering %>% 
#   filter(!is.na(Age_group)) %>% 
#   filter(isotype == "IgG") %>% 
#   filter(analyte != "VACV") %>% 
#   filter(panel %in% "MPXV") %>% 
#   filter(panel_detail != "SPox") %>% 
#   select(-PIN, - data)  %>%
#   dplyr::group_by(panel, panel_detail, isotype, sampleID_metadata, Age_group,
#                   POX_VACC_STATUS_ML, POX_VACC_STATUS_SELF, MPX_VACC_STATUS_1, MPX_VACC_STATUS_2,
#                   MPX_vac_final, MPX_SYMPTOME, MPX_diagnose_final, real, pred_AllIgGIgM, pred_AllIgG,
#                   pred_EpiIgGIgM, pred_EpiIgG, pred_RFAllIgGIgM, predxgboostAllIgGIgM, pred_LDAAllIgGIgM,
#                   analyte) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L) 
# END DEV
####

####
# Generate plots on Clustering data to test, if the signal distribution is 
# different between the different panels

# Plot Results for ATI-N and Delta in younger patients < 40 years
plotATI_N_young <-
  dataClustering_post %>%
  left_join(dataSeropositive, by = "sampleID_metadata") %>% 
  filter(panel_detail != "SPox") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre","Pre SPox",  "MVA",  "MVA SPox", "MPXV", "MPXV SPox"),
                            labels = c("Pre acute", "Pre epi", "MVA acute", "MVA epi", "MPXV acute", "MPXV epi"),
                            ordered = TRUE)) %>% 
  filter(analyte %in% c("ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = dataIn)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  stat_compare_means(comparisons = list(c("Pre acute", "Pre epi"),
                                        c("MVA acute", "MVA epi"),
                                        c("MPXV acute", "MPXV epi")),label = "p.signif",
                     method = "t.test") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,6)) +
  scale_x_discrete(name = "") +
  scale_color_viridis_c(name = "Ratio seropositivity") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Plot Results for ATI-N and Delta in younger patients > 50 years
plotATI_N_old <-
  dataClustering_post %>% 
  left_join(dataSeropositive, by = "sampleID_metadata") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(!(Age_group %in% c("< 30", "< 40", "< 50"))) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre","Pre SPox",  "MVA",  "MVA SPox", "MPXV", "MPXV SPox"),
                            labels = c("Pre acute", "Pre epi", "MVA acute", "MVA epi", "MPXV acute", "MPXV epi"),
                            ordered = TRUE)) %>% 
  filter(analyte %in% c("ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = dataIn)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  stat_compare_means(comparisons = list(c("Pre acute", "Pre epi"),
                                        c("MVA acute", "MVA epi"),
                                        c("MPXV acute", "MPXV epi")),label = "p.signif",
                     method = "t.test") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,6)) +
  scale_x_discrete(name = "") +
  scale_color_viridis_c(name = "Ratio seropositivity") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotCombined_ATI_N <- 
  ggarrange(plotATI_N_young, plotATI_N_old, ncol = 2, align = "hv",
            common.legend = TRUE, labels = c("g", "h"))

#ggsave(filename ="output/plotCombined_ATI_N.png", plotCombined_ATI_N , width = 7, height = 5, 
#       dpi = 600)



# Plot Results for relevant antigens in younger patients < 40 years
plotYoungAntigens <- 
  dataClustering_post %>% 
  left_join(dataSeropositive, by = "sampleID_metadata") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  # filter(POX_VACC_STATUS_SELF != "Yes" | is.na(POX_VACC_STATUS_SELF)) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox"),
                            labels = c("Pre acute", "Pre epi", "MVA acute", "MVA epi", "MPXV acute", "MPXV epi"),
                            ordered = TRUE)) %>% 
  # filter(panel_ATI %in% c("Pre acute", "Pre epi", "MPXV acute", "MPXV epi")) %>% 
  filter(analyte %in% c("D8L", "E8L", "A33R", "A35R", "B5R", "B6R", "ATI-N", "Delta")) %>% 
  filter(data > 0) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = data)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(1,5.5)) +
  scale_x_discrete(name = "") +
  scale_color_viridis_c(name = "Ratio seropositivity") +
  stat_compare_means(comparisons = list(c("Pre acute", "Pre epi"), c("MVA acute", "MVA epi"), c("MPXV acute", "MPXV epi")),label = "p.signif",
                     method = "t.test") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Plot Results for relevant antigens in older patients > 40 years
plotOldAntigens <-
  dataClustering_post %>% 
  left_join(dataSeropositive, by = "sampleID_metadata") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(!(Age_group %in% c("< 30", "< 40", "< 50"))) %>% 
  # filter(POX_VACC_STATUS_SELF != "Yes" | is.na(POX_VACC_STATUS_SELF)) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox"),
                            labels = c("Pre acute", "Pre epi", "MVA acute", "MVA epi", "MPXV acute", "MPXV epi"),
                            ordered = TRUE)) %>% 
  # filter(panel_ATI %in% c("Pre acute", "Pre epi", "MPXV acute", "MPXV epi")) %>% 
  filter(analyte %in% c("D8L", "E8L", "A33R", "A35R", "B5R", "B6R", "ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = data)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(1,5.5)) +
  scale_x_discrete(name = "") +
  scale_color_viridis_c(name = "Ratio seropositivity") +
  stat_compare_means(comparisons = list(c("Pre acute", "Pre epi"), c("MVA acute", "MVA epi"), c("MPXV acute", "MPXV epi")),label = "p.signif",
                     method = "t.test") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotCombined <- 
  ggarrange(plotYoungAntigens, plotOldAntigens, ncol = 2, align = "hv",
            common.legend = TRUE, labels = c("e", "f"))

ggsave(filename ="output/plotAntigens.png", plotCombined, width = 10, height = 8, 
       dpi = 600)




##
# Plot results of ratios for homologue pairs of antigens, 
# stratified for D8L serostatus
plotRatiosYoung <-
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  mutate(ratio_A29L_A27L = dataIn_A29L/dataIn_A27L,
         ratio_M1R_L1R = dataIn_M1R/dataIn_L1R,
         ratio_E8L_D8L = dataIn_E8L/dataIn_D8L,
         ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
         ratio_B6R_B5R = dataIn_B6R/dataIn_B5R) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox")) %>% 
  filter(isotype == "IgG") %>% 
  select(serostatus_D8L, panel_ATI, sampleID_metadata, starts_with("ratio")) %>% 
  pivot_longer(starts_with("ratio_"), names_to = "analyte", names_prefix = "ratio_",
               values_to = "ratio") %>% 
  mutate(analyte = factor(analyte, levels = c("A29L_A27L", "M1R_L1R", "E8L_D8L",
                                              "A35R_A33R", "B6R_B5R"), 
                          labels = c("A29L/A27L", "M1R/L1R", "E8L/D8L",
                                     "A35R/A33R", "B6R/B5R"))) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = panel_ATI, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), 
                                        c("MVA", "MVA SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotRatiosOld <-
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  filter(!(Age_group %in% c("< 30", "< 40", "< 50"))) %>% 
  mutate(ratio_A29L_A27L = dataIn_A29L/dataIn_A27L,
         ratio_M1R_L1R = dataIn_M1R/dataIn_L1R,
         ratio_E8L_D8L = dataIn_E8L/dataIn_D8L,
         ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
         ratio_B6R_B5R = dataIn_B6R/dataIn_B5R) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox")) %>% 
  filter(isotype == "IgG") %>% 
  select(serostatus_D8L, panel_ATI, sampleID_metadata, starts_with("ratio")) %>% 
  pivot_longer(starts_with("ratio_"), names_to = "analyte", names_prefix = "ratio_",
               values_to = "ratio") %>% 
  mutate(analyte = factor(analyte, levels = c("A29L_A27L", "M1R_L1R", "E8L_D8L",
                                              "A35R_A33R", "B6R_B5R"), 
                          labels = c("A29L/A27L", "M1R/L1R", "E8L/D8L",
                                     "A35R/A33R", "B6R/B5R"))) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = panel_ATI, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), 
                                        c("MVA", "MVA SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotCombinedRatios <- 
  ggarrange(plotRatiosYoung, plotRatiosOld, ncol = 2, align = "hv",
            common.legend = TRUE)

ggsave(filename ="output/plotCombinedRatios.png", plotCombinedRatios, width = 10, height = 8, 
       dpi = 600)


####
# Perform ROC Analysis based on the wide dataframe
# Generate two dataframes: first dataframe for the prediction on the pre-panel
# Second dataframe: for the prediction on the MVA-panel based on ATI-N alonw
dataInWideAntigens <- 
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "MPXV SPox"))

dataInWideATI_N <- 
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "MVA", "MPXV SPox")) %>% 
  mutate(panel_ATI_N = if_else(panel_ATI == "MPXV SPox", panel_ATI, "Pre/MVA"))

##
# Perform ROC-Analysis for the following antigens:
# Delta, D8L, E8L, A33R, A35R, B5R, B6R, ratioAG_seropositive
roc_Delta <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_Delta)

roc_D8L <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_D8L)

roc_E8L <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_E8L)

roc_A33R <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_A33R)

roc_A35R <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_A35R)

roc_B5R <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_B5R)

roc_B6R <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$dataIn_B6R)

roc_RatioAG <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$ratioAG_seropositive)

roc_RatioAG_clear <-
  roc(factor(dataInWideAntigens$panel_ATI, levels = c("Pre", "MPXV SPox"),
             ordered = TRUE), dataInWideAntigens$ratioAG_clear_seropositive)


##
# Perform ROC-Analysis for ATI-N for the following antigen
# ATI-N
roc_ATI_N <-
  roc(factor(as.character(dataInWideATI_N$panel_ATI_N), levels = c("Pre/MVA", "MPXV SPox"),
             ordered = TRUE),dataInWideATI_N$`dataIn_ATI-N`)




##
# Determine cutoff values for the MVA-panel based on the following ratios:
# E8L/D8L, A35R/A33R, B6R/B5R
# The Ratios A29L/A27L as well as M1R/L1R were excluded as they were not
# very conclusive.
# To determine the ROC with regard to the population-based cutoffs, the MVA
# panel is set as the control panel, the MPXV SPox is set as the case panel

# Generate dataframe with ratios
dataInWideRatios <- 
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  mutate(ratio_A29L_A27L = dataIn_A29L/dataIn_A27L,
         ratio_M1R_L1R = dataIn_M1R/dataIn_L1R,
         ratio_E8L_D8L = dataIn_E8L/dataIn_D8L,
         ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
         ratio_B6R_B5R = dataIn_B6R/dataIn_B5R) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("MVA", "MPXV SPox")) %>% 
  filter(isotype == "IgG") %>% 
  mutate(panel_ATI = factor(panel_ATI, levels = c("MVA", "MPXV SPox"),
                            ordered = TRUE))




##
# Perform ROC-Analysis for the following antigens:
# E8L/D8L, A35R/A33R, B6R/B5R
roc_E8L_D8L <-
  roc(dataInWideRatios$panel_ATI, dataInWideRatios$ratio_E8L_D8L)

roc_A35R_A33R <-
  roc(dataInWideRatios$panel_ATI, dataInWideRatios$ratio_A35R_A33R)

roc_B6R_B5R <-
  roc(dataInWideRatios$panel_ATI, dataInWideRatios$ratio_B6R_B5R)




##
# Read Assay parameters and determine cutoff values
# Calculate 95% confidence interval
k = 1
ci_list <- list()
for(k in c(1:length(ls(all.names = TRUE, pattern = ("^roc_"))))){
  ci_list[[k]] <- ci.coords(get(ls(all.names = TRUE, pattern = ("^roc_"))[k]), x="best", input = "threshold", ret=rets,  best.method = c("closest.topleft"),
                            best.policy = "random")
  assign(paste("ci",  ls(all.names = TRUE, pattern = ("^roc_"))[k], sep = "_"), ci_list[[k]])
  k = k+1
}

# Export results of roc-analysis
k = 1
param_list <- list()
for(k in c(1:length(ls(all.names = TRUE, pattern = ("^ci_roc"))))){
  param_list[[k]] <- readRocParam(ls(all.names = TRUE, pattern = ("^ci_roc"))[k])
  k = k+1
}
# Save parameters in a dataframe
paramDataFrame <- rocParamFunc(param_list)
paramDataFrame_transpose_sep <- transformRocParamFunc(param_list, "IFA IgG")
sens <- prepDataFrameFunc(paramDataFrame_transpose_sep, "sensitivity", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "sensitivity")
spec <- prepDataFrameFunc(paramDataFrame_transpose_sep, "specificity", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "specificity")
acc <- prepDataFrameFunc(paramDataFrame_transpose_sep, "accuracy", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "accuracy")

# Determine threshold 
threshold <-
  paramDataFrame_transpose_sep %>%
  dplyr::select(antigene, CI, threshold) %>%
  mutate(antigene = str_remove(antigene, "ci_roc_")) %>% 
  mutate(antigene = if_else(antigene == "ATI_N", "ATI-N", antigene))

threshold_median <-
  threshold %>% 
  filter(CI == "median")

# For highest specificites -> Lower for more excluded values from training
# min(roc_ATI_N$thresholds[roc_ATI_N$specificities == 1])

##
# Predict on the SPox-Panel Pre and the SPox-Panel MVA
predictSPox <- 
  dataClustering_post %>% 
  filter(isotype == "IgG") %>% 
  mutate(panel_cutoff = case_when(panel == "MPXV" ~ "MPXV",
                                  panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                                  panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                                  panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                                  panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                                  panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_cutoff = factor(panel_cutoff, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                               ordered = TRUE)) %>% 
  left_join(threshold_median, by = c("analyte" = "antigene")) %>%
  filter(!is.na(threshold)) %>% 
  rowwise() %>% 
  mutate(SPox_Pre_prediction = if_else(dataIn > threshold, "positive", "negative"),
         SPox_Pre_prediction = case_when(analyte == "ATI-N" & dataIn > 
                                           min(roc_ATI_N$thresholds[roc_ATI_N$specificities == 1]) ~ "positive",
                                         analyte == "ATI-N" & dataIn <= 
                                           min(roc_ATI_N$thresholds[roc_ATI_N$specificities == 1]) ~ "negative",
                                         
                                         analyte == "A33R" & dataIn > 
                                           min(roc_A33R$thresholds[roc_A33R$specificities == 1]) ~ "positive",
                                         analyte == "A33R" & dataIn <= 
                                           min(roc_A33R$thresholds[roc_A33R$specificities == 1]) ~ "negative",
                                         
                                         analyte == "Delta" & dataIn > 
                                           min(roc_Delta$thresholds[roc_Delta$specificities == 1]) ~ "positive",
                                         analyte == "Delta" & dataIn <= 
                                           min(roc_Delta$thresholds[roc_Delta$specificities == 1]) ~ "negative",
                                         
                                         analyte == "A35R" & dataIn > 
                                           min(roc_A35R$thresholds[roc_A35R$specificities == 1]) ~ "positive",
                                         analyte == "A35R" & dataIn <= 
                                           min(roc_A35R$thresholds[roc_A35R$specificities == 1]) ~ "negative",
                                         
                                         analyte == "B5R" & dataIn > 
                                           min(roc_B5R$thresholds[roc_B5R$specificities == 1]) ~ "positive",
                                         analyte == "B5R" & dataIn <= 
                                           min(roc_B5R$thresholds[roc_B5R$specificities == 1]) ~ "negative",
                                         
                                         analyte == "B6R" & dataIn > 
                                           min(roc_B6R$thresholds[roc_B6R$specificities == 1]) ~ "positive",
                                         analyte == "B6R" & dataIn <= 
                                           min(roc_B6R$thresholds[roc_B6R$specificities == 1]) ~ "negative",
                                         
                                         analyte == "D8L" & dataIn > 
                                           min(roc_D8L$thresholds[roc_D8L$specificities == 1]) ~ "positive",
                                         analyte == "D8L" & dataIn <= 
                                           min(roc_D8L$thresholds[roc_D8L$specificities == 1]) ~ "negative",
                                         TRUE ~ SPox_Pre_prediction))



predictSPox_Ratios <- 
  dataInWide %>% 
  #  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40", "< 50"))) %>% 
  mutate(ratio_A29L_A27L = dataIn_A29L/dataIn_A27L,
         ratio_M1R_L1R = dataIn_M1R/dataIn_L1R,
         ratio_E8L_D8L = dataIn_E8L/dataIn_D8L,
         ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
         ratio_B6R_B5R = dataIn_B6R/dataIn_B5R) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  select(panel, panel_detail, sampleID_metadata, Age_group, panel_ATI, starts_with("ratio_")) %>% 
  unique() %>% 
  pivot_longer(cols = starts_with("ratio_"), names_prefix = "ratio_",
               names_to = "antigene", values_to = "ratio") %>% 
  left_join(threshold_median, by = c("antigene")) %>% 
  filter(!is.na(threshold)) %>% 
  rowwise() %>% 
  mutate(SPox_Ratio_prediction = if_else(ratio > threshold, "positive", "negative"),
         SPox_Ratio_prediction = case_when(antigene == "E8L_D8L" & ratio > 
                                             min(roc_E8L_D8L$thresholds[roc_E8L_D8L$specificities == 1]) ~ "positive",
                                           antigene == "E8L_D8L" & ratio <= 
                                             min(roc_E8L_D8L$thresholds[roc_E8L_D8L$specificities == 1]) ~ "negative",
                                           
                                           antigene == "A35R_A33R" & ratio > 
                                             min(roc_A35R_A33R$thresholds[roc_A35R_A33R$specificities == 1]) ~ "positive",
                                           antigene == "A35R_A33R" & ratio <= 
                                             min(roc_A35R_A33R$thresholds[roc_A35R_A33R$specificities == 1]) ~ "negative",
                                           
                                           antigene == "B6R_B5R" & ratio > 
                                             min(roc_B6R_B5R$thresholds[roc_B6R_B5R$specificities == 1]) ~ "positive",
                                           antigene == "B6R_B5R" & ratio <= 
                                             min(roc_B6R_B5R$thresholds[roc_B6R_B5R$specificities == 1]) ~ "negative",
                                           TRUE ~ SPox_Ratio_prediction)) #%>% 
#select(-threshold) %>% 
#pivot_wider(names_from = antigene, values_from = c(ratio, SPox_Ratio_prediction))




##
# Plot ROC curves
scaleFUN <- function(x) sprintf("%.1f", x)

plotROCfunction <- function(inputList, subIn){
  ggroc(inputList) +
    geom_abline(slope = 1, intercept = 1, linetype = 2, colour = "red")+
    labs(colour = "")+
    scale_x_reverse(labels = scaleFUN) +
    scale_y_continuous(labels = scaleFUN) +
    scale_color_colorblind(name = "", drop = F) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    # labs(subtitle = subIn) +
    theme_bw() +
    theme(plot.subtitle=element_text(hjust=0.5),
          legend.position = "top") 
} 

plotrocAntigens <- 
  plotROCfunction(list("Delta" = roc_Delta,
                       "D8L" = roc_D8L,
                       "E8L" = roc_E8L,
                       "A33R" = roc_A33R,
                       "A35R" = roc_A35R,
                       "B5R" = roc_B5R,
                       "B6R" = roc_B6R,
                       "ATI-N" = roc_ATI_N), "IgG")

ggsave(filename = "output/plotROCAntigen.png", plotrocAntigens,
       width = 6, height = 6, dpi = 600)


plotrocRatios <- 
  plotROCfunction(list("E8L/D8L" = roc_E8L_D8L,
                       "A35R/A33R" = roc_A35R_A33R,
                       "B6R/B5R" = roc_B6R_B5R), "IgG")

ggsave(filename = "output/plotrocRatios.png", plotrocRatios,
       width = 6, height = 6, dpi = 600)

predict_Pre_Young_Spox <-
  predictSPox %>% 
  # filter(panel_cutoff == "Pre SPox") %>% 
  filter(analyte %in% c("Delta", "D8L", "E8L", "A33R", "A35R", "B5R", "B6R", "ATI-N")) 



##
# Generate plot with ratios of positive/negative in the different age groups
# and the different panels for the different analytes
plotAntigensYoungPre <-
  predict_Pre_Young_Spox %>% 
  filter(panel_detail != "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_cutoff %in% c("Pre", "Pre SPox")) %>% 
  filter(!is.na(Age_group)) %>% 
  group_by(panel_cutoff, analyte) %>% 
  summarise(Prop_Positive = mean(SPox_Pre_prediction == "positive"),
            N = n(),
            Std_Error = sqrt(Prop_Positive*(1-Prop_Positive)/N)) %>% 
  mutate( Lower = Prop_Positive - 1.96*Std_Error, 
          Upper = Prop_Positive + 1.96*Std_Error) %>% 
  ungroup() %>% 
  ggplot(aes(x = panel_cutoff, y = Prop_Positive)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Panel") +
  scale_y_continuous(name = "Proportion Seropositive") +
  facet_wrap("analyte", scales = "free_y") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotAntigensYoungMVA <-
  predict_Pre_Young_Spox %>% 
  filter(panel_detail != "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_cutoff %in% c("MVA", "MVA SPox")) %>% 
  filter(!is.na(Age_group)) %>% 
  group_by(panel_cutoff, analyte) %>% 
  summarise(Prop_Positive = mean(SPox_Pre_prediction == "positive"),
            N = n(),
            Std_Error = sqrt(Prop_Positive*(1-Prop_Positive)/N)) %>% 
  mutate( Lower = Prop_Positive - 1.96*Std_Error, 
          Upper = Prop_Positive + 1.96*Std_Error) %>% 
  ungroup() %>% 
  ggplot(aes(x = panel_cutoff, y = Prop_Positive)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Panel") +
  scale_y_continuous(name = "Proportion Seropositive") +
  facet_wrap("analyte", scales = "free_y") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



plotRatiosMVA <- 
  predictSPox_Ratios %>% 
  filter(panel_detail != "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox")) %>% 
  filter(!is.na(Age_group)) %>% 
  group_by(panel_ATI, antigene) %>% 
  summarise(Prop_Positive = mean(SPox_Ratio_prediction == "positive"),
            N = n(),
            Std_Error = sqrt(Prop_Positive*(1-Prop_Positive)/N)) %>% 
  mutate( Lower = Prop_Positive - 1.96*Std_Error, 
          Upper = Prop_Positive + 1.96*Std_Error) %>% 
  ungroup() %>% 
  ggplot(aes(x = panel_ATI, y = Prop_Positive)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  scale_x_discrete(name = "Panel") +
  scale_y_continuous(name = "Proportion Seropositive") +
  facet_wrap("antigene", scales = "free_y") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))




plotPredictYoungCombined <-
  ggarrange(plotAntigensYoungPre, plotAntigensYoungMVA,plotRatiosMVA, align = "hv",
            ncol = 3, labels = "auto")

ggsave(filename = "output/plotPredictSeropos.png", width = 12, height = 6,
       dpi = 600)


plotSPoxSelfReporteVACC <-
  predict_Pre_Young_Spox %>% 
  filter(panel_detail != "SPox") %>% 
  filter(POX_VACC_STATUS_SELF %in% c("Yes", "No")) %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_cutoff %in% c("Pre", "Pre SPox", "MVA SPox")) %>% 
  filter(!is.na(Age_group)) %>% 
  group_by(panel_cutoff, analyte, POX_VACC_STATUS_SELF) %>% 
  summarise(Prop_Positive = mean(SPox_Pre_prediction == "positive"),
            N = n(),
            Std_Error = sqrt(Prop_Positive*(1-Prop_Positive)/N)) %>% 
  mutate( Lower = Prop_Positive - 1.96*Std_Error, 
          Upper = Prop_Positive + 1.96*Std_Error) %>% 
  ungroup() %>% 
  ggplot(aes(x = panel_cutoff, y = Prop_Positive, fill = POX_VACC_STATUS_SELF)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = "dodge") +
  scale_x_discrete(name = "Panel") +
  scale_fill_manual(name = "Self reported Smallpox VACC", values = colorblind_pal()(8)[2:8]) +
  scale_y_continuous(name = "Proportion Seropositive") +
  facet_wrap("analyte") +
  theme_pubr() +
  #ä scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(filename = "output/plotPredictSeroposPOXVACC.png", 
       plotSPoxSelfReporteVACC, width = 7, height = 7,
       dpi = 600)

## Select and filter ratios
predict_ratio_all <- 
  predictSPox_Ratios %>% 
  filter(panel == "SPox") %>% 
  unique() %>%  
  select(-threshold) %>% 
  pivot_wider(names_from = antigene, values_from = c(ratio, SPox_Ratio_prediction)) %>% 
  select(sampleID_metadata, starts_with("SPox"))


## Filter SPOx panel for predictions on young positives
predict_all <-
  predict_Pre_Young_Spox %>% 
  select(sampleID_metadata, panel, panel_detail, panel_cutoff, analyte, Age_group,
         POX_VACC_STATUS_SELF, MPX_VACC_STATUS_1, MPX_VACC_STATUS_2, MPX_vac_final,
         MPX_SYMPTOME, MPX_diagnose_final, real:preLDA_2, SPox_Pre_prediction) %>%
  filter(panel == "SPox") %>% 
  unique() %>% 
  pivot_wider(names_from = analyte, values_from = SPox_Pre_prediction) %>% 
  mutate(real = factor(real),
         preLDA_1 = factor(preLDA_1)) %>% 
  left_join(predict_ratio_all, by = "sampleID_metadata")





## Filter only young positives
predict_young_pre_positive <-
  predict_all %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  filter(panel_detail == "Pre") %>% 
  rowwise() %>% 
  mutate(sum_positive = sum(D8L == "positive", 
                            E8L == "positive",
                            A33R == "positive",
                            A35R == "positive",
                            B5R == "positive",
                            B6R == "positive",
                            Delta == "positive",
                            `ATI-N` == "positive"),
         sum_Ratio_positive = sum(SPox_Ratio_prediction_E8L_D8L == "positive",
                                  SPox_Ratio_prediction_A35R_A33R == "positive",
                                  SPox_Ratio_prediction_B6R_B5R == "positive") ,
         fraq_pred = sum(!is.na(pred_AllIgGIgM),
                         !is.na(pred_AllIgG),
                         !is.na(pred_EpiIgGIgM),
                         !is.na(pred_EpiIgG)))


# Plot for prediction on young positive
plotPredictYoungPos <- 
  predict_young_pre_positive %>% 
  ggplot(mapping = aes(x = sum_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotPredictYoungPosRatio <- 
  predict_young_pre_positive %>% 
  ggplot(mapping = aes(x = sum_Ratio_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotPredictYoungPosPosRatio <- 
  predict_young_pre_positive %>% 
  ggplot(mapping = aes(x = sum_positive, y = sum_Ratio_positive, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



# Predict for MPXV positive
predict_young_mpxv_positive <-
  predict_all %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  filter(panel_detail == "MPXV") %>% 
  rowwise() %>% 
  mutate(sum_positive = sum(D8L == "positive", 
                            E8L == "positive",
                            A33R == "positive",
                            A35R == "positive",
                            B5R == "positive",
                            B6R == "positive",
                            Delta == "positive",
                            `ATI-N` == "positive"),
         sum_Ratio_positive = sum(SPox_Ratio_prediction_E8L_D8L == "positive",
                                  SPox_Ratio_prediction_A35R_A33R == "positive",
                                  SPox_Ratio_prediction_B6R_B5R == "positive") ,
         fraq_pred = sum(!is.na(pred_AllIgGIgM),
                         !is.na(pred_AllIgG),
                         !is.na(pred_EpiIgGIgM),
                         !is.na(pred_EpiIgG)))

plotPredictYounMPXVgPos <- 
  predict_young_mpxv_positive %>% 
  ggplot(mapping = aes(x = sum_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  # facet_wrap("MPX_vac_final") +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotPredictYounMPXVgPosRatio <- 
  predict_young_mpxv_positive %>% 
  ggplot(mapping = aes(x = sum_Ratio_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  # facet_wrap("MPX_vac_final") +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotPredictYounMPXVgPosPosRatio <- 
  predict_young_mpxv_positive %>% 
  ggplot(mapping = aes(x = sum_positive, y = sum_Ratio_positive, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  # facet_wrap("MPX_vac_final") +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


predict_young_MVA_positive <-
  predict_all %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  filter(panel_detail == "MVA") %>% 
  rowwise() %>% 
  mutate(sum_positive = sum(D8L == "positive", 
                            E8L == "positive",
                            A33R == "positive",
                            A35R == "positive",
                            B5R == "positive",
                            B6R == "positive",
                            Delta == "positive",
                            `ATI-N` == "positive"),
         sum_positive_ATI = sum(`ATI-N` == "positive") ,
         sum_Ratio_positive = sum(SPox_Ratio_prediction_E8L_D8L == "positive",
                                  SPox_Ratio_prediction_A35R_A33R == "positive",
                                  SPox_Ratio_prediction_B6R_B5R == "positive") ,
         fraq_pred = sum(!is.na(pred_AllIgGIgM),
                         !is.na(pred_AllIgG),
                         !is.na(pred_EpiIgGIgM),
                         !is.na(pred_EpiIgG)))



plotPredictYoungMVAPosATI <- 
  predict_young_MVA_positive%>% 
  ggplot(mapping = aes(x = as.factor(sum_positive_ATI), y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotPredictYoungMVAPosAll <- 
  predict_young_MVA_positive%>% 
  ggplot(mapping = aes(x = sum_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotPredictYoungMVAPosAllRatio <- 
  predict_young_MVA_positive%>% 
  ggplot(mapping = aes(x = sum_Ratio_positive, y = fraq_pred, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotPredictYoungMVAPosAllPosRatio <- 
  predict_young_MVA_positive%>% 
  ggplot(mapping = aes(x = sum_positive, y = sum_Ratio_positive, color = preLDA_1)) +
  geom_jitter(width = 0.2) +
  #  theme_pubr() +
  scale_color_manual(name = "ML LDA Prediction", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotPredictCombined <-
  ggarrange(plotPredictYoungPos, plotPredictYounMPXVgPos,
            plotPredictYoungMVAPosAll, plotPredictYoungMVAPosATI, 
            common.legend = TRUE, align = "hv", ncol = 2, nrow = 2, 
            labels = "auto")

plotPredictCombinedRatio <-
  ggarrange(plotPredictYoungPosRatio, plotPredictYounMPXVgPosRatio,
            plotPredictYoungMVAPosAllRatio, 
            common.legend = TRUE, align = "hv", ncol = 2, nrow = 2, 
            labels = "auto")

ggsave(file = "output/plotPredictCombined.png", 
       plotPredictCombined, width = 8, height = 8, 
       dpi = 600)

young_pre_export <-
  predict_young_pre_positive %>% 
  filter(sum_positive > 0 | sum_Ratio_positive > 1) %>% 
  mutate(sum_positive_ATI = NA)

young_MVA_export <-
  predict_young_MVA_positive %>% 
  filter((sum_positive_ATI == 1 | sum_Ratio_positive > 1) & sum_positive > 0)


young_positive_export <-
  rbind(young_pre_export, young_MVA_export) 

export(young_positive_export, file = "output/Pre_MVA_Positive_SPox_2.xlsx")
sampleIDsExcluds <- data.frame(excludeIDs = sort(unique(young_positive_export$sampleID_metadata)))

export(sampleIDsExcluds, file = "output/Samples_Pre_MVA_Positive_SPox_2.csv")

include <-
  unique(predict_all$sampleID_metadata)[!unique(predict_all$sampleID_metadata) %in% 
                                          sort(unique(young_positive_export$sampleID_metadata))]

##
# Test if exclusion improves prediction

predConfAll <-
  predict_all

predConfExcl <-
  predict_all %>% 
  filter(sampleID_metadata %in% include)

confusionMatrix(predConfAll$real, predConfAll$preLDA_1)
confusionMatrix(predConfExcl$real, predConfExcl$preLDA_1)
confusionMatrix(predConfAll$real, as.factor(predConfAll$predxgboost_1))
confusionMatrix(predConfExcl$real, as.factor(predConfExcl$predxgboost_1))
confusionMatrix(predConfAll$real, as.factor(predConfAll$predRF_1))
confusionMatrix(predConfExcl$real, as.factor(predConfExcl$predRF_1))


