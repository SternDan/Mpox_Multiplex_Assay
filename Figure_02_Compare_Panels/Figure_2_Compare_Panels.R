##
# Figure 2 - Compare Panels
# Daniel Stern
# Robert Koch Institute
# ZBS 3
# Version 1
# Last modified 2025-03-01
##
rm(list = ls(all.names = TRUE))

library(tidyverse)
library(ggthemes)
library(ggpubr)
#library(ggforce)
library(GGally)
#library(see)


##
# Load data input 
load("input/dataInputComparePanels.Rdata")


##
# Plot panels distribution as used for ML-algorithms
# New plot: parallel plot with three groups
## Stratify for childhoodimmu
plotIgG <-
  dataInputComparePanels %>% 
#  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunisation", values = colorblind_pal()(8)) +
  theme_pubclean() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotIgM <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         H3, 
         A33, 
         A35,
         B5, 
         B6, 
         A5, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunisation", values = colorblind_pal()(8)) +
  theme_pubclean() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotStratifiedChildhood <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = c("a", "b"))

ggsave("output/plotStratified.png", plotStratifiedChildhood, width = 16, height = 6,
       dpi = 600)


####
# Plot comparison panels
dataInIgGAll <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         analyte, dataIn) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  rename(status = panel_detail) %>% 
  group_by(analyte) %>% 
  mutate(dataInNorm = (dataIn-min(dataIn))/(max(dataIn)-min(dataIn))) %>% #min max normalization for each antigen
  ungroup() %>% 
  filter(childhoodImmuAge != "Ambiguous")



# Visualize results as box plots
plotIgGSpox <-
  dataInIgGAll %>% 
  filter(analyte %in% c("E8", "A35", "B6", "ATI-N")) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Spox") %>% 
  ggplot( aes(x = status, y = dataInNorm, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotIgGAcute <-
  dataInIgGAll %>% 
  filter(analyte %in% c("E8", "A35", "B6", "ATI-N")) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Acute") %>% 
  ggplot( aes(x = status, y = dataInNorm, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotCombinedSerogroupsIgG <- 
  ggarrange(plotIgGAcute, plotIgGSpox, ncol = 2, align = "hv", common.legend = F,
            labels = c("c", "d"))

#### Supporting figure with other antigens
plotIgGSpoxSup <-
  dataInIgGAll %>% 
  filter(!(analyte %in% c("E8", "A35", "B6", "ATI-N"))) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Spox") %>% 
  ggplot( aes(x = status, y = dataInNorm, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotIgGAcuteSup <-
  dataInIgGAll %>% 
  filter(!(analyte %in% c("E8", "A35", "B6", "ATI-N"))) %>% 
  filter(isotype == "IgG") %>% 
  filter(panel %in% "Acute") %>% 
  ggplot( aes(x = status, y = dataInNorm, fill = status)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serogroup", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  facet_grid(childhoodImmuAge ~ analyte) +
  theme_pubr() +
  stat_compare_means(aes(group = childhoodImmuAge), comparisons = list(c("Pre", "MVA"), c("MVA", "Mpox"),
                                                                       c("Pre", "Mpox")),method = "t.test",
                     label = "p.signif") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotCombinedSerogroupsIgGSup <- 
  ggarrange(plotIgGAcuteSup, plotIgGSpoxSup, nrow = 2, align = "hv", common.legend = FALSE,
            labels = c("a", "b"))

ggsave(file = "output/plotAntigensFigS2.png", plotCombinedSerogroupsIgGSup,
       width = 14, height = 12, dpi = 600)

####
# Plot ratios of different antigens 
plotRatioAcuteEpiYes <-
  dataInputComparePanels  %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(antigens %in% c("A35/A33", "B6/B5", "E8/D8")) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios, fill = panel_detail)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serostatus", values = colorblind_pal()(8)[2:8]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotRatioAcuteEpiNo <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="No") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(antigens %in% c("A35/A33", "B6/B5", "E8/D8")) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios, fill = panel_detail)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serostatus", values = colorblind_pal()(8)[2:8]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotRatios <-
  ggarrange(plotRatioAcuteEpiNo, plotRatioAcuteEpiYes, ncol = 2, align = "hv", labels = c("e", "f"))


# Plot supporting figure ratios
plotRatioAcuteEpiYesSupp <-
  dataInputComparePanels  %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(!(antigens %in% c("A35/A33", "B6/B5", "E8/D8"))) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios, fill = panel_detail)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serostatus", values = colorblind_pal()(8)[2:8]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

plotRatioAcuteEpiNoSupp <-
  dataInputComparePanels %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  filter(childhoodImmuAge =="No") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27, 
         A29,
         L1, 
         M1, 
         D8, 
         E8,
         A33, 
         A35,
         B5, 
         B6) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1/L1` = M1/L1,
         `A29/A27` = A29/A27,
         `B6/B5` = B6/B5, 
         `A35/A33` = A35/A33,
         `E8/D8` = E8/D8) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1/L1`, `A29/A27`, `A35/A33`, `B6/B5`,  `E8/D8`) %>% 
  pivot_longer(cols = c("M1/L1", "A29/A27", "A35/A33", "B6/B5", "E8/D8"), names_to = "antigens", values_to = "ratios") %>% 
  filter(!(antigens %in% c("A35/A33", "B6/B5", "E8/D8"))) %>% 
  ggplot(mapping = aes(x = panel_detail, y = ratios, fill = panel_detail)) +
  geom_boxplot() +
  scale_fill_manual(name = "Serostatus", values = colorblind_pal()(8)[2:8]) +
  theme_pubr() +
  facet_grid(. ~ antigens) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Pre", "MVA"),
                                                                               c("MVA", "Mpox"),
                                                                               c("Pre", "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")


plotRatiosSupp <-
  ggarrange(plotRatioAcuteEpiNoSupp, plotRatioAcuteEpiYesSupp, ncol = 2, align = "hv", labels = c("a", "b"))

ggsave(file = "output/plotRatiosFigS3.png", plotRatiosSupp, width = 8, height = 5,
       dpi = 600)



plotFigLower <-
  ggarrange( plotCombinedSerogroupsIgG,
             plotRatios, nrow = 2, align = "hv",
             heights = c(2,1.4),
             common.legend = T)
plotFig2New <- 
  ggarrange(plotStratifiedChildhood, plotFigLower,
             nrow = 2, align = "hv",
            common.legend = F, heights = c(1,1.7))

# Save figure 2
ggsave("output/Figure_2.png", plotFig2New, width = 14, height = 12,
       dpi = 600)
ggsave("output/Figure_2.pdf", plotFig2New, width = 14, height = 12,
       dpi = 600)

# Save supporting Figures
ggsave("output/plotIgGSpoxSup.png", plotIgGSpoxSup , width = 10, height = 10, dpi = 600)
ggsave("output/plotIgGAcuteSup.png", plotIgGAcuteSup , width = 10, height = 10, dpi = 600)

# Save supporting Figure with ratios
ggsave("output/plotRatiosSupp.png",plotRatiosSupp , width = 8, height = 4, dpi = 600)
