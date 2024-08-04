##
# Compare serological panel seroconversion
# Daniel Stern
# Robert Koch Institute
# ZBS 3
# Version Final
# Last modified 2024-08-04
##

rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggforce)
library(GGally)
library(brms)
library(tidybayes)
library(see)



##
# Load metadata for Age-groups
## Load metadata
## Load datainput generated for clustering
load("../24 Spox Pre Pos Clustering/input/dataClustering.Rdata")
load("../24 Spox Pre Pos Clustering/output/plotCombinedAll.Rdata")
load("../25 Three-Way ANOVA panel/output/boxplotsAnova.Rdata")

## Unify the format for the Age_groups varibale
dataInputMeta <- 
  dataClustering %>% 
  dplyr::filter(isotype != "IgA") %>% 
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
                            ordered = TRUE)) %>% 
  select(-starts_with("pred")) %>%
  select(sampleID_metadata, Age_group) %>% 
  unique()


##
# Load data: All data
dataInput <-
  dataClustering %>% 
  mutate(serostatus_clear = if_else(serostatus_cat == "positive", 1, 0)) %>% 
  mutate(panel_strat = case_when(panel == "SPox" ~ "SPox",
                                 panel %in% c("MPXV", "MVA", "Pre", "Pre_New") ~ "Acute")) %>% 
  mutate(panel_detail = if_else(panel == "Pre_New", "Pre", panel_detail),
         panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "MPXV", 
                                                        "SPox", "SPox_Rep", "CPXV"),
                               labels = c("Pre", "MVA", "Mpox", "SPox", "SPox_Rep", "CPXV"),
                               ordered = TRUE))  %>% 
  mutate(childhoodImmuAge = case_when(Age_group %in% c("18-29", "30-39") ~ "No",
                                      Age_group %in% c("40-49") ~ "Ambiguous",
                                      Age_group %in% c("50-59", "60+") ~ "Yes",
                                      Age_group %in% c("< 30", "< 40") ~ "No",
                                      Age_group %in% c("< 50") ~ "Ambiguous",
                                      Age_group %in% c("< 60", "60+") ~ "Yes"))

save(dataInput, file = "output/dataInputPanelANOVA.Rdata")


##
# Plot panels distribution as used for ML-algorithms
# New plot: parallel plot with three groups
## Stratify for childhoodimmu
plotIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotIgA <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgA")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotStratifiedChildhood <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = c("a", "b"),
            common.legend = TRUE)


####
# Plot ratios of different antigens 
plotRatioIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         A33R, 
         A35R,
         B5R, 
         B6R) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1R/L1R` = M1R/L1R,
         `A29L/A27L` = A29L/A27L,
         `B6R/B5R` = B6R/B5R, 
         `A35R/A33R` = A35R/A33R,
         `E8L/D8L` = E8L/D8L) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1R/L1R`, `A29L/A27L`, `A35R/A33R`, `B6R/B5R`,  `E8L/D8L`) %>% 
  ggparcoord(columns = c(6:10), alphaLines = 0.1, groupColumn = "panel_detail",
             scale = "center", boxplot = F) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)) +
  theme_bw() +
  facet_grid(. ~ childhoodImmuAge, labeller = as_labeller(c("1" = "Yes", "2" = "No", "3" = "Ambiguous"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotRatioIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         A33R, 
         A35R,
         B5R, 
         B6R) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1R/L1R` = M1R/L1R,
         `A29L/A27L` = A29L/A27L,
         `B6R/B5R` = B6R/B5R, 
         `A35R/A33R` = A35R/A33R,
         `E8L/D8L` = E8L/D8L) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(!is.infinite(`M1R/L1R`)) %>% 
  filter(!is.infinite(`A29L/A27L`)) %>% 
  filter(!is.infinite(`B6R/B5R`)) %>% 
  filter(`M1R/L1R` < 3) %>% 
  filter(`A29L/A27L` < 3) %>% 
  filter(`B6R/B5R` < 3) %>% 
  filter(`A35R/A33R` < 3) %>% 
  filter(`E8L/D8L` < 3) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1R/L1R`, `A29L/A27L`, `A35R/A33R`, `B6R/B5R`,  `E8L/D8L`) %>% 
  ggparcoord(columns = c(6:10), alphaLines = 0.1, groupColumn = "panel_detail",
             scale = "center", boxplot = F) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)) +
  theme_bw() +
  facet_grid(. ~ childhoodImmuAge, labeller = as_labeller(c("1" = "Yes", "2" = "No", "3" = "Ambiguous"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotRatios <-
  ggarrange(plotRatioIgG, plotRatioIgM, ncol = 2, align = "hv", labels = c("a", "b"),
            common.legend = TRUE)


# New plot Fig 3 with ANOVA
plotFig3New <- 
  ggarrange(plotStratifiedChildhood, boxplotsAnovaCombined, nrow = 2, align = "hv",
            common.legend = F, heights = c(3, 6))

ggsave("output/plotFig3New.png", plotFig3New, width = 14, height = 18,
       dpi = 600)

ggsave("output/plotFig3New.pdf", plotFig3New, width = 16, height = 22,
       dpi = 600)

# Save Fig with ratios
ggsave("output/Sup_Fig_Ratios.png", plotRatios, width = 10, height = 3, dpi = 600)





