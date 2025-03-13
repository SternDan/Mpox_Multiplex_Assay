####
# Generate Figure 4
# Daniel Stern RKI
# 2025-03-01
# Version 1.0
####

# Clean environment
rm(list = ls(all.names = TRUE))

# Load libraries
library(rio)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(circlize)
library(data.table)
library(plyr)


# Load data
# load("input/dataIn.Rdata")
load("input/dataInMeta.Rdata")
# load("input/dataInputAnalyte.Rdata")
# load("input/dataInRep.Rdata")
load("input/statisticalDataCombined.Rdata")


##
# Plot data
# Plot performance of ML-algorithms
plotPerformance <- 
  statisticalDataCombined %>% 
  filter(parameter %in% c("mean", "stdev")) %>%
  pivot_longer(cols = c(accuracy, precision, recall, fscore), names_to = "score", values_to = "values") %>%
  mutate(score = if_else(score == "fscore", "F1 score", score)) %>% 
  pivot_wider(names_from = parameter, values_from = values) %>% 
  filter(!(grepl("thresh", algorithm))) %>% 
  filter(!(seropanel %in% c("all acute", "all epi"))) %>% 
  mutate(algorithm = factor(algorithm, levels = c("XGBoost", "LDA", "RF", "LDA RF", "FRBC", "LDA FRBC"), 
                            labels = c("GBC","LDA",  "RF", "LDA RF", "FRBC", "LDA FRBC"),
                            ordered = TRUE)) %>% 
  mutate(seropanel = factor(seropanel, levels = c("acute acute", "epi epi",
                                                  "acute epi", "epi acute", "all all"), ordered = TRUE)) %>% 
  mutate(score = factor(score, levels = c("accuracy", "precision", "recall", "F1 score"),
                        ordered = TRUE)) %>% 
  filter(score == "F1 score") %>% 
  ggplot(mapping = aes(x = algorithm, color = antibody)) +
  geom_line(aes(y = mean, group = antibody)) +
  geom_errorbar(aes(group = algorithm, ymin = mean - stdev, ymax = mean + stdev),
                width = 0.5) +
  facet_grid(. ~ seropanel) +
  scale_color_manual(name = "Isotype", values = colorblind_pal()(8)[2:8]) +
  theme_pubclean() +
  scale_x_discrete(name = "Algorithm") +
  scale_y_continuous(name = "Metric", breaks = seq(0.2:1.0, by = 0.1)) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# Generate plot for misclassification in NK panel
plotNK <-
  dataInMeta %>% 
  filter(grepl("NK", sampleID_metadata)) %>% 
  filter(seropanel == "all all") %>% 
  filter(model %in% c("xgboost",
                      "lda")) %>% 
  mutate(model = factor(model, levels = c("xgboost",
                                          "lda"),
                        labels = c("GBC",
                                   "LDA"),
                        ordered = TRUE),
         pred = factor(pred, levels = c("Pre", "MVA", "Mpox"),
                       labels = c("Pre", "MVA", "Mpox"),
                       ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = Age_group)) +
  geom_bar(aes(fill=as.factor(pred)), position="fill") +
  scale_fill_manual(name = "Prediction", values = colorblind_pal()(8)[2:8]) +
  facet_grid(antibody ~ model) +
  theme_pubr()+
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(name = "Frequency Prediction") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Empty plot for later inclusion of cirular plots of misclassifications
plotNKdummy <-
  ggarrange(plotNK, ncol = 2)

plotFig4ab <-
  ggarrange(plotPerformance, plotNKdummy, nrow = 2, align = "hv", labels = c("auto"))

ggsave("output/Fig4ab.pdf", plotFig4ab, width = 12, height = 8, dpi = 600)


####
# Generate supporting table with combined performance data
supTable <- 
  statisticalDataCombined %>% 
  pivot_longer(cols = c(accuracy, precision, recall, fscore), names_to = "score", values_to = "values") %>%
  mutate(score = if_else(score == "fscore", "F1 score", score)) %>% 
  pivot_wider(names_from = parameter, values_from = values) %>% 
  filter(!(grepl("thresh", algorithm))) %>% 
  filter(!(seropanel %in% c("all acute", "all epi"))) %>% 
  mutate(algorithm = factor(algorithm, levels = c("XGBoost", "RF", "LDA RF", "LDA", "FRBC", "LDA FRBC"), 
                            ordered = TRUE)) %>% 
  mutate(seropanel = factor(seropanel, levels = c("acute acute", "epi epi", "all all",
                                                  "acute epi", "epi acute"), ordered = TRUE)) %>% 
  mutate(score = factor(score, levels = c("accuracy", "precision", "recall", "F1 score"),
                        ordered = TRUE)) %>% 
  select(seropanel, algorithm, antibody, score, mean, stdev) %>% 
  pivot_wider(names_from = algorithm, values_from = c("mean", "stdev")) %>% 
  mutate(LDA = paste(round(mean_LDA, 2), round(stdev_LDA, 2), sep = " \u00B1"),
         RF = paste(round(mean_RF, 2), round(stdev_RF, 2), sep = " \u00B1"),
         GBC = paste(round(mean_XGBoost, 2), round(stdev_XGBoost, 2), sep = " \u00B1"),
         FRBC = paste(round(mean_FRBC, 2), round(stdev_FRBC, 2), sep = " \u00B1"),
         `LDA+FRBC` = paste(round(`mean_LDA FRBC`, 2), round(`stdev_LDA FRBC`, 2), sep = " \u00B1"),
         `LDA+RF` = paste(round(`mean_LDA RF`, 2), round(`stdev_LDA RF`, 2), sep = " \u00B1")) %>% 
  select(seropanel, antibody, score, LDA, RF, GBC, FRBC, `LDA+FRBC`, `LDA+RF`)

save(file = "output/supTablePerformanceRev.Rdata", supTable)

####
# Generate circular plot based on frequency tables
# Calculate frequency tables
freqTableAllAllXGBoostIgG <-
  dataInMeta %>% 
  dplyr::select(model, antibody, seropanel, childhoodImmuAge, real, pred, agreement) %>% 
  filter(!is.na(real)) %>% 
  filter(seropanel == "all all" & model == "xgboost" & antibody == "IgG") %>% 
  group_by(model, antibody, seropanel, childhoodImmuAge, real) %>% 
  count() %>%  # Use plyr count here!
  group_by(model, antibody, seropanel, childhoodImmuAge, real) %>% 
  dplyr::mutate(freqSum = sum(freq))


freqTableAll <-
  dataInMeta %>% 
  dplyr::select(model, antibody, seropanel, childhoodImmuAge, real, pred, agreement) %>% 
  filter(!is.na(real)) %>% 
  group_by(model, antibody, seropanel, childhoodImmuAge) %>% 
  count() %>% 
  group_by(model, antibody, seropanel, childhoodImmuAge, real) %>% 
  dplyr::mutate(freqSum = sum(freq)) %>% 
  dplyr::rename(predFreq = freq, realFreq = freqSum) %>% 
  mutate(model = factor(model, levels = c("xgboost",
                                          "rf", 
                                          "lda_rf",
                                          "lda",
                                          "frbc",
                                          "lda_frbc",
                                          "frbc_threshold",
                                          "frbc_threshold_cut",
                                          "lda_threshold", 
                                          "lda_treshold_cut"),
                        labels = c("XGBoost",
                                   "RF", 
                                   "LDA RF",
                                   "LDA",
                                   "FRBC",
                                   "LDA FRBC",
                                   "FRBC threshold",
                                   "FRBC threshold cut",
                                   "LDA threshold", 
                                   "LDA treshold cut"),
                        ordered = TRUE))


save(freqTableAll, file = "output/freqTableAll.Rdata")

freqTableGBCLDA <-
  dataInMeta %>% 
  dplyr::select(model, antibody, seropanel, childhoodImmuAge, real, pred, agreement) %>% 
  filter(!is.na(real)) %>% 
  filter(seropanel == "all all" & model %in% c("xgboost", "lda")) %>% 
  filter(childhoodImmuAge %in% c("Yes", "No")) %>% 
  group_by(model, antibody, seropanel, childhoodImmuAge) %>% 
  count() %>% 
  group_by(model, antibody, seropanel, childhoodImmuAge, real) %>% 
  dplyr::mutate(freqSum = sum(freq)) %>% 
  dplyr::rename(predFreq = freq, realFreq = freqSum) %>% 
  mutate(model = factor(model, levels = c("xgboost",
                                          "rf", 
                                          "lda_rf",
                                          "lda",
                                          "frbc",
                                          "lda_frbc",
                                          "frbc_threshold",
                                          "frbc_threshold_cut",
                                          "lda_threshold", 
                                          "lda_treshold_cut"),
                        labels = c("XGBoost",
                                   "RF", 
                                   "LDA RF",
                                   "LDA",
                                   "FRBC",
                                   "LDA FRBC",
                                   "FRBC threshold",
                                   "FRBC threshold cut",
                                   "LDA threshold", 
                                   "LDA treshold cut"),
                        ordered = TRUE))

export(freqTableGBCLDA, "output/freqTableGBCLDA.xlsx")


### Build a function to generate circular plots
## Generate dataframe for plotting
circPlot <- function(seropanelIn, modelIn, antibodyIn, childhoodImmuAgeIn){
  dataPlot <-
    dataInMeta %>% 
    dplyr::select(model, antibody, seropanel, childhoodImmuAge, real, pred, agreement) %>% 
    filter(!is.na(real)) %>% 
    filter(seropanel %in% seropanelIn) %>% 
    filter(model %in% modelIn) %>% 
    filter(antibody == antibodyIn) %>% 
    group_by(model, antibody, seropanel, childhoodImmuAge) %>% 
    count() %>% 
    group_by(model, antibody, seropanel, childhoodImmuAge, real) %>% 
    dplyr::mutate(freqSum = sum(freq)) %>% 
    filter(childhoodImmuAge %in% childhoodImmuAgeIn) %>% 
    ungroup() %>% 
    select(real, pred, freq) %>% 
    mutate(real = as.factor(real)) %>% 
    mutate(color = case_when(real == "Pre" ~ colorblind_pal()(8)[2],
                             real == "MVA" ~ colorblind_pal()(8)[3],
                             real == "Mpox" ~ colorblind_pal()(8)[4]))
  
  #  pdf(paste("output/",seropanelIn, modelIn, antibodyIn, childhoodImmuAgeIn, 
  #            ".pdf", sep = ""))
  chordDiagram(dataPlot, preAllocateTracks = list(track.height = 0.25),
               annotationTrack = c("grid", "name"),
               grid.col = setNames(dataPlot$color, dataPlot$real))
  # dev.off()
}

line = -2
cex = 1.5
adj  = 0.5#0.025

pdf("output/Fig4_all_xgboost_lda.pdf", width=12, height=6)
par(mfrow = c(2, 4), mai=c(0.1, 0.1, 0.1, 0.1), mar = c(0.1, 0.1, 0.1, 0.1),
    oma = c(0.1, 0.1, 0.1, 0.1))


circPlot(c("all all"), c("xgboost"), c("IgG"), c("Yes"))
title(outer=outer,adj=0.025,main="c",cex.main=2,col="black",font=2,line=-1)
title(outer=outer,adj=adj,main="Childhood vaccinated IgG",cex.main=cex,col="black",font=1,line=line)
circPlot(c("all all"), c("xgboost"), c("IgG"), c("No"))
title(outer=outer,adj=adj,main="Naive IgG",cex.main=cex,col="black",font=1,line=line)

circPlot(c("all all"), c("xgboost"), c("IgM IgG"), c("Yes"))
title(outer=outer,adj=adj,main="Childhood vaccinated IgM IgG",cex.main=cex,col="black",font=1,line=line)
circPlot(c("all all"), c("xgboost"), c("IgM IgG"), c("No"))
title(outer=outer,adj=adj,main="Naive IgM IgG",cex.main=cex,col="black",font=1,line=line)


circPlot(c("all all"), c("lda"), c("IgG"), c("Yes"))
title(outer=outer,adj=0.025,main="d",cex.main=2,col="black",font=2,line=-1)
title(outer=outer,adj=adj,main="Childhood vaccinated IgG",cex.main=cex,col="black",font=2,line=line)
circPlot(c("all all"), c("lda"), c("IgG"), c("No"))
title(outer=outer,adj=adj,main="Naive IgG",cex.main=cex,col="black",font=1,line=line)

circPlot(c("all all"), c("lda"), c("IgM IgG"), c("Yes"))
title(outer=outer,adj=adj,main="Childhood vaccinated IgM IgG",cex.main=cex,col="black",font=2,line=line)
circPlot(c("all all"), c("lda"), c("IgM IgG"), c("No"))
title(outer=outer,adj=adj,main="Naive IgM IgG",cex.main=cex,col="black",font=1,line=line)

dev.off()







