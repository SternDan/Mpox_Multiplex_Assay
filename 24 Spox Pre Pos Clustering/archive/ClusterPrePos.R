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

rm(list = ls(all.names = TRUE))

## Load datainput generated for clustering
load("input/dataClustering.Rdata")

## Unify the format for the Age_groups varibale
dataClustering <- 
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



####
## Prepare wide dataframes for further analysis
## SPox dataframe
dataInWideSpox <-
  dataClustering %>% 
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
  dataClustering %>% 
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
  dataClustering %>% 
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
  dataClustering %>% 
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
  dataClustering %>%
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
                            ordered = TRUE)) %>% 
  filter(analyte %in% c("ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = dataIn)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"),
                                        c("MVA", "MVA SPox"),
                                        c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_pubr() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Plot Results for ATI-N and Delta in younger patients > 50 years
plotATI_N_old <-
  dataClustering %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
                            ordered = TRUE)) %>% 
  filter(analyte %in% c("ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = dataIn)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"),
                                        c("MVA", "MVA SPox"),
                                        c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_pubr() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotCombined_ATI_N <- 
  ggarrange(plotATI_N_young, plotATI_N_old, ncol = 2, align = "hv",
            common.legend = TRUE)

ggsave(filename ="output/plotCombined_ATI_N.png", plotCombined_ATI_N , width = 7, height = 5, 
       dpi = 600)



# Plot Results for relevant antigens in younger patients < 40 years
plotYoungAntigens <- 
  dataClustering %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MPXV", "MPXV SPox")) %>% 
  filter(analyte %in% c("D8L", "E8L", "A33R", "A35R", "B5R", "B6R")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = data)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Plot Results for relevant antigens in older patients > 40 years
plotOldAntigens <-
  dataClustering %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MPXV", "MPXV SPox")) %>% 
  filter(analyte %in% c("D8L", "E8L", "A33R", "A35R", "B5R", "B6R")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = data)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  geom_jitter(aes(color = ratioAG_seropositive), width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  scale_color_viridis_c() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotCombined <- 
  ggarrange(plotYoungAntigens, plotOldAntigens, ncol = 2, align = "hv",
            common.legend = TRUE)

ggsave(filename ="output/plotAntigens.png", plotCombined, width = 10, height = 8, 
       dpi = 600)




##
# Plot results of ratios for homologue pairs of antigens, 
# stratified for D8L serostatus
plotRatiosYoung <-
  dataInWide %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
# Perform Probit Analysis to separate sera in different groups according to 
# the following antigens and/or ratios stratified by age < 40 or >= 50
# 1) ATI-N:   Groups  -> Pre/MVA combined
#                     -> MPXV SPox
# 2) Delta, D8L, E8L, A33R, A35R, B5R, B6R -> only for age < 40 years
#             Groups  -> Pre 
#                     -> MPXV SPox
# 3) Ratios: A29L/A27L, M1R/L1R, E8L/D8L, A35R/A33R, B6R/B5R  
#             Groups  -> Pre/MVA combined
#                     -> MPXV SPox

#####
# Probit Analysis:
# # https://community.rstudio.com/t/plotting-probit-regression-with-ggplot2/7843/2

probX = function(p, model) {
  data.frame(prob=p, 
             xval = (qnorm(p) - coef(model)[1])/coef(model)[2])
}
scaleFUN <- function(x) sprintf("%.1f", x)

# Generate datainput with all antigens
dataInputProbitAntigens <-   
  dataClustering %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
  left_join(dataSeropositive, by = "sampleID_metadata") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(panel_detail != "SPox") %>% 
  #  filter(Age_group %in% c("< 30", "< 40")) %>% 
  # filter(POX_VACC_STATUS_SELF != "Yes" | is.na(POX_VACC_STATUS_SELF)) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(analyte %in% c("ATI-N", "Delta", "D8L", "E8L", "A33R", "A35R", "B5R", "B6R")) %>% 
  filter(isotype == "IgG")

# Generate dataframe for probit analysis on ATI-N
dataIn_Probit_ATI_N <-
  dataInputProbitAntigens %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_ATI %in% c("Pre", "MVA", "MPXV SPox")) %>% 
  filter(analyte %in% c("ATI-N")) %>% 
  mutate(timeProbit = if_else(panel_ATI == "MPXV SPox", 1, 0))

logr_vm_N <- glm(timeProbit ~ dataIn, data=dataIn_Probit_ATI_N, family=binomial(link="probit"))

d_N = probX(c(0.25,0.5,0.95), logr_vm_N)

# Save GLM
save(logr_vm_N, file = "output/Probit_ATI-N_Antigen_glm.RData")

# Berechnung der Werte auf der Kurve mit Konfidenzintervallen
# https://www.geo.fu-berlin.de/en/v/soga/Basics-of-statistics/Logistic-Regression/Logistic-Regression-in-R---An-Example/index.html
probs = predict(logr_vm_N, 
                newdata = data.frame(dataIn = dataIn_Probit_ATI_N$dataIn), 
                type = "response", 
                se.fit = TRUE)

# Speichern in einem großen gemeinsamen Dataframe
outDataframe <-
  dataIn_Probit_ATI_N %>%
  add_column(as_tibble(probs$fit)) %>%
  rename(Fit = value) %>%
  add_column(as_tibble(probs$se)) %>%
  rename(SE = value) %>%
  mutate(upper = Fit + 1.96*SE,
         lower = Fit - 1.96*SE) 

plotProbit_ATI_N_Norm <- 
  outDataframe %>%
  ggplot(mapping = aes(x=dataIn, y = timeProbit)) +
  geom_point(aes(shape = as.factor(panel_ATI)), size = 2) +
  stat_smooth(method="glm", method.args=list(family=binomial(link="probit")) , se = T, fullrange = T,
              color = colorblind_pal()(8)[4]) +
  theme_pubr() +
  scale_shape_manual(name = "ATI-N IgG", values = c(1, 19,2)) +
  geom_segment(data=d_N, aes(x=xval, xend=xval, y=0, yend=prob), colour=colorblind_pal()(8)[3]) +
  geom_segment(data=d_N, aes(x=0, xend=xval, y=prob, yend=prob), colour=colorblind_pal()(8)[3]) +
  geom_text(data=d_N, aes(label=round((xval), 2), x=xval, y=-0.03), size=3, colour=colorblind_pal()(8)[3]) +
  scale_x_continuous(name = "DataIn", expand = c(0,0)) +
  scale_y_continuous(name = "Panels", breaks = c(0,1), labels = c("Pre/MVA", "MPXV SPox"), expand = c(0.05,0))+
  annotation_logticks(sides = "b")+
  scale_color_manual(name = "ATI-N IgG", values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank())

ggsave("output/plotProbit_ATI_N_Norm.png", plotProbit_ATI_N_Norm, width = 6, height = 4, dpi = 600)

# ROC Auswertung mit den gefitteten Daten durchführen
rocIgM_IgG_ratio <-
  roc(as.factor(outDataframe$timeProbit), outDataframe$Fit)



plot(rocATI_N)

# Plot ROC Curves
plotROC <-
  ggroc(list("ATI-N" = rocIgM_IgG_ratio)) +
  geom_abline(slope = 1, intercept = 1, linetype = 2, colour = "red")+
  labs(colour = "")+
  scale_x_reverse(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  scale_color_manual(name = "Method", values = colorblind_pal()(8)[2:8]) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  theme_bw()

ggsave("output/plotROC_ATI_N.png", plotROC, width = 5.5, height = 4, dpi = 600)



#####
# DEV
# Develop function for the probit-Analysis
# Generate dataframe for probit analysis on ATI-N

# ageFilter <- c("< 30", "< 40")
# panelIn <- c("Pre", "MVA", "MPXV SPox")
# analyteIn <- c("ATI-N")
# saveModel <- c("output/Probit_ATI-N_Antigen_glm.RData")
# savePlot <- c("output/plotProbit_ATI_N_Norm.png")
# saveROC <- c("output/plotROC_ATI_N.png")

analyseProbit <- function(inputData, ageFilter, panelDataIn, 
                          analyteIn,
                          saveModel, savePlot,
                          saveROC, labelIn){
  dataIn_Probit_ATI_N <-
    inputData %>% 
    filter(Age_group %in% ageFilter) %>% 
    filter(panel_ATI %in% panelDataIn) %>% 
    filter(analyte %in% analyteIn) %>% 
    mutate(timeProbit = if_else(panel_ATI == "MPXV SPox", 1, 0))
  
  logr_vm_N <- glm(timeProbit ~ dataIn, data=dataIn_Probit_ATI_N, family=binomial(link="probit"))
  d_N = probX(c(0.25,0.5,0.95), logr_vm_N)
  
  # Save GLM
  save(logr_vm_N, file = saveModel)
  probs = predict(logr_vm_N, 
                  newdata = data.frame(dataIn = dataIn_Probit_ATI_N$dataIn), 
                  type = "response", 
                  se.fit = TRUE)
  
  # Speichern in einem großen gemeinsamen Dataframe
  outDataframe <-
    dataIn_Probit_ATI_N %>%
    add_column(as_tibble(probs$fit)) %>%
    rename(Fit = value) %>%
    add_column(as_tibble(probs$se)) %>%
    rename(SE = value) %>%
    mutate(upper = Fit + 1.96*SE,
           lower = Fit - 1.96*SE) 
  
  plotProbit_ATI_N_Norm <- 
    outDataframe %>%
    ggplot(mapping = aes(x=dataIn, y = timeProbit)) +
    geom_point(aes(shape = as.factor(panel_ATI)), size = 2) +
    stat_smooth(method="glm", method.args=list(family=binomial(link="probit")) , se = T, fullrange = T,
                color = colorblind_pal()(8)[4]) +
    theme_pubr() +
    scale_shape_manual(name = analyteIn, values = c(1, 19,2)) +
    geom_segment(data=d_N, aes(x=xval, xend=xval, y=0, yend=prob), colour=colorblind_pal()(8)[3]) +
    geom_segment(data=d_N, aes(x=0, xend=xval, y=prob, yend=prob), colour=colorblind_pal()(8)[3]) +
    geom_text(data=d_N, aes(label=round((xval), 2), x=xval, y=-0.03), size=3, colour=colorblind_pal()(8)[3]) +
    scale_x_continuous(name = "DataIn", expand = c(0,0)) +
    scale_y_continuous(name = "Panels", breaks = c(0,1), labels = labelIn, expand = c(0.05,0))+
    annotation_logticks(sides = "b")+
    scale_color_manual(name = analyteIn, values = colorblind_pal()(8)[2:8]) +
    theme(strip.background = element_blank())
  
  ggsave(savePlot, plotProbit_ATI_N_Norm, width = 6, height = 5, dpi = 600)
  
  # ROC Auswertung mit den gefitteten Daten durchführen
  rocIgM_IgG_ratio <-
    roc(as.factor(outDataframe$timeProbit), outDataframe$Fit)
  
  # Plot ROC Curves
  plotROC <-
    ggroc(list(rocIgM_IgG_ratio)) +
    geom_abline(slope = 1, intercept = 1, linetype = 2, colour = "red")+
    labs(colour = "")+
    scale_x_reverse(labels = scaleFUN) +
    scale_y_continuous(labels = scaleFUN) +
    scale_color_manual(name = "Method", values = colorblind_pal()(8)[2:8]) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    theme_bw()
  
  ggsave(saveROC, plotROC, width = 5.5, height = 4, dpi = 600)
}

analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MVA", "MPXV SPox"),
              analyteIn = c("ATI-N"),
              saveModel = c("output/Probit_ATI-N_Antigen_glm.RData"),
              savePlot = c("output/plotProbit_ATI_N_Norm.png"),
              saveROC = c("output/plotROC_ATI_N.png"),
              labelIn = c("Pre/MVA", "MPXV SPox"))


# 2) Delta, D8L, E8L, A33R, A35R, B5R, B6R -> only for age < 40 years
#             Groups  -> Pre 
#                     -> MPXV SPox


analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("Delta"),
              saveModel = c("output/Probit_Delta_glm.RData"),
              savePlot = c("output/plotProbit_Delta_Norm.png"),
              saveROC = c("output/plotROC_Delta.png"),
              labelIn = c("Pre", "MPXV SPox"))


analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("D8L"),
              saveModel = c("output/Probit_D8L_glm.RData"),
              savePlot = c("output/plotProbit_D8L_Norm.png"),
              saveROC = c("output/plotROC_D8L.png"),
              labelIn = c("Pre", "MPXV SPox"))


analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("E8L"),
              saveModel = c("output/Probit_E8L_glm.RData"),
              savePlot = c("output/plotProbit_E8L_Norm.png"),
              saveROC = c("output/plotROC_E8L.png"),
              labelIn = c("Pre", "MPXV SPox"))

analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("A33R"),
              saveModel = c("output/Probit_A33R_glm.RData"),
              savePlot = c("output/plotProbit_A33R_Norm.png"),
              saveROC = c("output/plotROC_A33R.png"),
              labelIn = c("Pre", "MPXV SPox"))

analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("A35R"),
              saveModel = c("output/Probit_A35R_glm.RData"),
              savePlot = c("output/plotProbit_A35R_Norm.png"),
              saveROC = c("output/plotROC_A35R.png"),
              labelIn = c("Pre", "MPXV SPox"))

analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("B5R"),
              saveModel = c("output/Probit_B5R_glm.RData"),
              savePlot = c("output/plotProbit_B5R_Norm.png"),
              saveROC = c("output/plotROC_B5R.png"),
              labelIn = c("Pre", "MPXV SPox"))

analyseProbit(inputData = dataInputProbitAntigens, 
              ageFilter = c("< 30", "< 40"), 
              panelDataIn = c("Pre", "MPXV SPox"),
              analyteIn = c("B6R"),
              saveModel = c("output/Probit_B6R_glm.RData"),
              savePlot = c("output/plotProbit_B6R_Norm.png"),
              saveROC = c("output/plotROC_B6R.png"),
              labelIn = c("Pre", "MPXV SPox"))

## Load Models
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

probit_ATI_N <- loadRData("output/Probit_ATI-N_Antigen_glm.RData")
probit_Delta <- loadRData("output/Probit_Delta_glm.RData")
probit_D8L <- loadRData("output/Probit_D8L_glm.RData")
probit_E8L <- loadRData("output/Probit_E8L_glm.RData")
probit_B5R <- loadRData("output/Probit_B5R_glm.RData")
probit_B6R <- loadRData("output/Probit_B6R_glm.RData")
probit_A33R <- loadRData("output/Probit_A33R_glm.RData")
probit_A35R <- loadRData("output/Probit_A35R_glm.RData")

probs_Delta = predict(probit_Delta, 
                      newdata = data.frame(dataIn = dataInWide$dataIn_Delta), 
                      type = "response", 
                      se.fit = TRUE)

probs_ATI_N = predict(probit_ATI_N, 
                      newdata = data.frame(dataIn = dataInWide$`dataIn_ATI-N`), 
                      type = "response", 
                      se.fit = TRUE)

probs_D8L = predict(probit_D8L, 
                    newdata = data.frame(dataIn = dataInWide$dataIn_E8L), 
                    type = "response", 
                    se.fit = TRUE)

probs_E8L = predict(probit_E8L, 
                    newdata = data.frame(dataIn = dataInWide$dataIn_E8L), 
                    type = "response", 
                    se.fit = TRUE)

probs_B5R = predict(probit_B5R, 
                    newdata = data.frame(dataIn = dataInWide$dataIn_B5R), 
                    type = "response", 
                    se.fit = TRUE)

probs_B6R = predict(probit_B6R, 
                    newdata = data.frame(dataIn = dataInWide$dataIn_B6R), 
                    type = "response", 
                    se.fit = TRUE)

probs_A33R = predict(probit_A33R, 
                     newdata = data.frame(dataIn = dataInWide$dataIn_A33R), 
                     type = "response", 
                     se.fit = TRUE)

probs_A35R = predict(probit_A35R, 
                     newdata = data.frame(dataIn = dataInWide$dataIn_A35R), 
                     type = "response", 
                     se.fit = TRUE)


outDataframe <-
  dataInWide %>%
  add_column(as_tibble(probs_Delta$fit)) %>%
  rename(Fit_Delta = value) %>% 
  add_column(as_tibble(probs_ATI_N$fit)) %>%
  rename(Fit_ATI_N = value) %>% 
  add_column(as_tibble(probs_D8L$fit)) %>%
  rename(Fit_D8L = value) %>% 
  add_column(as_tibble(probs_E8L$fit)) %>%
  rename(Fit_E8L = value) %>% 
  add_column(as_tibble(probs_B5R$fit)) %>%
  rename(Fit_B5R = value) %>% 
  add_column(as_tibble(probs_B6R$fit)) %>%
  rename(Fit_B6R = value) %>% 
  add_column(as_tibble(probs_A33R$fit)) %>%
  rename(Fit_A33R = value) %>% 
  add_column(as_tibble(probs_A35R$fit)) %>% 
  rename(Fit_A35R = value)

##
# 

export(outDataframe,
       file = "output/outDataframe.csv")
##
# Analyse impact of differen 
analyseData <-
  outDataframe %>% 
  rowwise() %>% 
  mutate(sum_pred = sum(!is.na(pred_AllIgGIgM), !is.na(pred_AllIgG), !is.na(pred_EpiIgGIgM), !is.na(pred_EpiIgG)),
         agreement = do.call(`==`, lapply(c(pmin, pmax), do.call, c(across(pred_AllIgGIgM:pred_EpiIgG), na.rm = TRUE))),
         consensus = case_when(!is.na(pred_AllIgGIgM) ~ pred_AllIgGIgM,
                               !is.na(pred_AllIgG) ~ pred_AllIgG,
                               !is.na(pred_EpiIgGIgM) ~ pred_EpiIgGIgM,
                               !is.na(pred_EpiIgG) ~ pred_EpiIgG))
# pred_clear = if_else(pred_AllIgGIgM == pred))

analyseData %>% 
  filter(panel == "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40")) %>% 
  filter(panel_clustering == "Pre positive") %>% 
  ggplot(aes(x = consensus)) +
  geom_bar(position = "stack")

analyseData %>% 
  filter(panel == "SPox") %>% 
  filter(panel_clustering == "Pre positive") %>% 
  group_by(Age_group, panel_clustering, predxgboostAllIgGIgM) %>% 
  summarise(n = n(),
            mean_D8L = mean(Fit_D8L), 
            mean_E8L = mean(Fit_E8L), 
            mean_ATI_N = mean(Fit_ATI_N))
count() #%>% 
ggplot(aes(x = Age_group, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap("panel_clustering") +
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Proportion Self-reported Pox Vaccinated",
                     limits = c(0,1)) +
  theme_pubr()




analyseData %>% 
  filter(panel == "SPox") %>% 
  filter(panel_clustering == "Pre positive") %>% 
  pivot_longer(starts_with("Fit_"), names_to = "antigen", names_prefix = "Fit_", 
               values_to = "probability") %>% 
  ggplot(mapping = aes(x = Age_group, y = probability, fill = as.factor(sum_pred))) +
  geom_boxplot() +
  geom_jitter( aes(color = as.factor(sum_pred)), alpha = 0.5, position = position_dodge(width = 0.75)) +
#  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), 
#                                        c("MVA", "MVA SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
#                     method = "t.test") +
  facet_wrap("antigen") +
  theme_pubr() +
  scale_color_manual(values = colorblind_pal()(8)[2:8]) +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


analyseData %>% 
  filter(panel == "SPox") %>% 
  filter(panel_clustering == "Pre positive") %>% 
  filter(Fit_E8L > 0.95) %>% 
  filter(sum_pred < 4) %>% 
  pivot_longer(starts_with("Fit_"), names_to = "antigen", names_prefix = "Fit_", 
               values_to = "probability") %>% 

  ggplot(mapping = aes(x = Age_group, y = probability, fill = as.factor(predxgboostAllIgGIgM))) +
  geom_boxplot() +
  geom_jitter( aes(color = as.factor(predxgboostAllIgGIgM)), alpha = 0.5, position = position_dodge(width = 0.75)) +
  #  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), 
  #                                        c("MVA", "MVA SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
  #                     method = "t.test") +
  facet_wrap("antigen") +
  theme_pubr() +
  scale_color_manual(values = colorblind_pal()(8)[2:8]) +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

####
# Perform ROC Analysis based on the wide dataframe
# Generate two dataframes: first dataframe for the prediction on the pre-panel
# Second dataframe: for the prediction on the MVA-panel based on ATI-N alonw
dataInWideAntigens <- 
  dataInWide %>% 
  select(-starts_with("pred")) %>% 
  unique() %>% 
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
  select(-starts_with("pred")) %>% 
  unique() %>% 
  filter(panel_detail != "SPox") %>% 
  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  filter(panel_clustering %in% c("Pre", "MVA", "MPXV SPox")) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  mutate(panel_ATI_N = if_else(panel_ATI %in% c("MPXV SPox"), "MPXV SPox", "Pre/MVA"))

##
# Perform ROC-Analysis for the following antigens:
# Delta, D8L, E8L, A33R, A35R, B5R, B6R, ratioAG_seropositive
rocDelta <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_Delta)

rocD8L <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_D8L)

rocE8L <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_E8L)

rocA33R <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_A33R)

rocA35R <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_A35R)

rocDB5R <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_B5R)

rocB6R <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_B6R)

rocDelta <-
  roc(as.character(dataInWideAntigens$panel_ATI), dataInWideAntigens$dataIn_Delta)


##
# Perform ROC-Analysis for ATI-N for the following antigen
# ATI-N


##
# Read Assay parameters and determine cutoff values


##
# Predict on the SPox-Panel Pre and the SPox-Panel MVA












#plotRatiosYoung <-
dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
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
  filter(panel_ATI %in% c("Pre SPox", "MVA SPox", "MPXV SPox")) %>% 
  filter(isotype == "IgG") %>% 
  select(serostatus_D8L, panel_ATI, sampleID_metadata, starts_with("ratio"), predxgboostAllIgGIgM) %>% 
  pivot_longer(starts_with("ratio"), names_to = "analyte", names_prefix = "ratio_",
               values_to = "ratio") %>% 
  mutate(analyte = factor(analyte, levels = c("A29L_A27L", "M1R_L1R", "E8L_D8L",
                                              "A35R_A33R", "B6R_B5R"), 
                          labels = c("A29L/A27L", "M1R/L1R", "E8L/D8L",
                                     "A35R/A33R", "B6R/B5R"))) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = panel_ATI, y = ratio)) +
  geom_boxplot(outlier.shape = NA, aes(fill = as.factor(predxgboostAllIgGIgM))) +
  facet_wrap("analyte") +
  theme_pubr() +
  scale_color_manual(values = colorblind_pal()(8)[2:8], name = "Serostatus") +
  scale_fill_manual(values = colorblind_pal()(8)[2:8], name = "Serostatus") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


dataInWide %>% 
  unique() %>% 
  filter(panel_detail != "SPox") %>% 
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
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MVA", "MVA SPox", "MPXV", "MPXV SPox")) %>% 
  filter(isotype == "IgG") %>% 
  select(serostatus_D8L, panel_ATI, sampleID_metadata, rat_seropositive = ratio_seropositive, starts_with("ratio"), predxgboostAllIgGIgM) %>% 
  pivot_longer(starts_with("ratio"), names_to = "analyte", names_prefix = "ratio_",
               values_to = "ratio") %>% 
  mutate(analyte = factor(analyte, levels = c("A29L_A27L", "M1R_L1R", "E8L_D8L",
                                              "A35R_A33R", "B6R_B5R"), 
                          labels = c("A29L/A27L", "M1R/L1R", "E8L/D8L",
                                     "A35R/A33R", "B6R/B5R"))) %>% 
  unique() %>% 
  ggplot(mapping = aes(x = panel_ATI, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap("analyte") +
  geom_jitter(aes(color = rat_seropositive), width = 0.1) +
  scale_color_viridis_c() +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), 
                                        c("MVA", "MVA SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  #scale_color_manual(values = colorblind_pal()(8)[2:8], name = "Serostatus") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))




dataClustering %>% 
  # select(-starts_with("pred")) %>% 
  unique() %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  # filter(POX_VACC_STATUS_SELF != "Yes" | is.na(POX_VACC_STATUS_SELF)) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("Pre", "Pre SPox", "MPXV", "MPXV SPox")) %>% 
  filter(analyte %in% c("D8L", "E8L", "A33R", "A35R", "B5R", "B6R")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = data)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  geom_jitter(aes(color = as.factor(predxgboostAllIgGIgM)), width = 0.1, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"), c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  theme_pubr() +
  scale_color_manual(values = colorblind_pal()(8)[2:8], name = "Serostatus") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


dataClustering %>%
  filter(panel_detail != "SPox") %>% 
  filter(!sampleID_metadata %in% c("P-22-01610-001-08", "S-20-004-007-02", "S-20-004-007-01")) %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre","Pre SPox",  "MVA",  "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  #  filter(panel_ATI %in% c("MVA", "MVA SPox", "MPXV")) %>% 
  filter(analyte %in% c("ATI-N", "Delta")) %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = panel_ATI, y = dataIn)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap ("analyte") +
  stat_compare_means(comparisons = list(c("Pre", "Pre SPox"),
                                        c("MVA", "MVA SPox"),
                                        c("MPXV", "MPXV SPox")),label = "p.signif",
                     method = "t.test") +
  geom_jitter(aes(color = as.factor(predxgboostAllIgGIgM)), width = 0.1, alpha = 0.5) +
  scale_color_manual(values = colorblind_pal()(8)[2:8], name = "Serostatus") +
  theme_pubr() +
  theme(strip.background = element_blank(),
        # legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))




# Results for younger subjects
dataInWideYoung <-
  dataInWide %>%   mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                                                panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                                                panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                                                panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                                                panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                                                panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
                          panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                                             ordered = TRUE)) %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50"))

dataPCAyoung <- 
  dataInWide %>%   mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                                                panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                                                panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                                                panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                                                panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                                                panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
                          panel_ATI = factor(panel_ATI, levels = c("Pre",  "Pre SPox", "MVA","MVA SPox", "MPXV", "MPXV SPox"),
                                             ordered = TRUE)) %>% 
  filter(Age_group %in% c("< 30", "< 40", "< 50")) %>% 
  select(starts_with("dataIn"))

PCAyoung <- prcomp(dataPCAyoung, scale = TRUE)
PCAyoungVis <- data.frame(PC1 = PCAyoung$x[,1],
                          PC2 = PCAyoung$x[,2], 
                          group = dataInWideYoung$panel_clustering)
plotPCAyoung <- 
  PCAyoungVis %>% 
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  #geom_density_2d(aes(group = group), show.legend = TRUE, linewidth = 1, alpha = 0.5) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = group)) +
  geom_point() +
  theme_bw() +
  scale_color_manual( values = colorblind_pal()(8)[c(2,3,4,8,5)]) +
  scale_fill_manual(values = colorblind_pal()(8)[c(2,3,4,8,5)])


##
# Analyse missclassifications of 
dataInWidePre <- 
  dataInWide %>% 
  filter(panel == "SPox" & panel_detail == "Pre") 




##
# Perform k-means clustering on the dataset
# https://uc-r.github.io/kmeans_clustering
set.seed(123)
dataPCAAll <-
  dataInWide %>% 
  select(starts_with("dataIn"))

res.km <- kmeans(scale(dataPCAyoung), 4, nstart = 25)
# K-means clusters showing the group of each individuals
res.km$cluster

fviz_cluster(res.km, data =dataPCAyoung,
             palette = colorblind_pal()(8)[c(2,3,4,5,8)], 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)



# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(dataPCAyoung, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(dataPCAAll, kmeans, method = "silhouette")

gap_stat <- clusGap(dataPCAyoung, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
