####
# Analyse revision measurements
# Daniel Stern
# 2025-03-18
# Version 2.0
####

rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(gplots)
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(scales)
library(caret)
library(pROC)
library(rstatix)
library(corrplot)
library(cvms)
library(ggimage)
library(rsvg)

# Load datainput
load("input/ensemblePrediction.Rdata")
load("input/heatmap_input.Rdata")
load("input/dataInputComparison.Rdata")



plotPreGrouped <-
  ensemblePrediction %>%
  mutate(serostatus.delta = factor(serostatus.delta, levels = c(0,1),
                                   labels = c("negative", "positive"))) %>% 
  pivot_longer(cols = c(Strat_Pre, Strat_MVA, Strat_MPXV), names_to = "pred",
               values_to = "ensembl") %>% 
  mutate(pred = factor(pred, levels = c("Strat_Pre", "Strat_MVA", "Strat_MPXV"),
                       labels = c("Pre", "MVA", "MPXV"),
                       ordered = TRUE), 
         year_group = ifelse(year_birth < 1972, "<1972", ">=1972"),
         year_group = factor(year_group, levels = c("<1972", ">=1972"),
                             ordered = TRUE)) %>% 
  filter(real == "Pre") %>%
  filter(!is.na(ensembl)) %>% 
  group_by(year_group, pred) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(year_group) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ggplot(aes(x = year_group, y = percent , fill = pred)) +
  geom_bar(stat = "identity", position = "fill") +  # relative Häufigkeit (100%)
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),  # Zahlen mittig im Balken
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) + # Farben für die Kategorien
  labs(
    x = "Year of birth",
    y = "Relative frequency",
    fill = "Prediction"
  ) +
  scale_y_continuous() +
  theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")




####
# Plot heatmap of raw data
set.seed(123)
matIgG <- as.matrix(scale(heatmap_input_IgG[2:16]))
haIgG <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgG$panel),
                           annotation_height = unit(4, "mm"))

pdf(file = "output/heatmap_IgG_white.pdf", width = 14, height = 10)
Heatmap(matIgG, km = 3,
        name = "Multiplex quantfied",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgG$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Serostatus Delta-VACV", col = c(colorblind_pal()(8)[5:6])) +
  Heatmap(heatmap_input_IgG$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c(viridis_pal()(5))) +
  Heatmap(factor(heatmap_input_IgG$panel, levels = c("Pre", "MVA", "MPXV"), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4])) +
  Heatmap(heatmap_input_IgG$GBC_Pre, name = "Pre GBC", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgG$LDA_Pre, name = "Pre LDA", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgG$Strat_Pre, name = "Pre ensemble", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgG$GBC_MVA, name = "MVA GBC", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgG$LDA_MVA, name = "MVA LDA", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgG$Strat_MVA, name = "MVA ensemble", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgG$GBC_MPXV, name = "MPXV GBC", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgG$LDA_MPXV, name = "MPXV LDA", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgG$Strat_MPXV, name = "MPXV ensemble", col = c("white", colorblind_pal()(8)[4])) 
dev.off()


matIgM <- as.matrix(scale(heatmap_input_IgM[2:16]))
ha <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgM$panel),
                        annotation_height = unit(4, "mm"))
pdf(file = "output/heatmap_IgM_white.pdf", width = 10, height = 8)
Heatmap(matIgM, km = 3,#clustering_distance_rows = "pearson",
        #col = c(viridis_pal()(20)),#c(colorblind_pal()(8)[c(6,5,8)]),
        name = "Multiplex quantfied",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgM$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Serostatus Delta-VACV", col = c(colorblind_pal()(8)[5:6])) +
  Heatmap(heatmap_input_IgM$n_MVA_vacc, name = "MVA vaccinations", na_col = "white", col = c(viridis_pal()(5))) +
  Heatmap(factor(heatmap_input_IgM$panel, levels = c("Pre", "MVA", "MPXV"), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4])) +
  Heatmap(heatmap_input_IgM$GBC_Pre, name = "Pre GBC", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgM$LDA_Pre, name = "Pre LDA", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgM$Strat_Pre, name = "Pre ensemble", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgM$GBC_MVA, name = "MVA GBC", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgM$LDA_MVA, name = "MVA LDA", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgM$Strat_MVA, name = "MVA ensemble", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgM$GBC_MPXV, name = "MPXV GBC", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgM$LDA_MPXV, name = "MPXV LDA", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgM$Strat_MPXV, name = "MPXV ensemble", col = c("white", colorblind_pal()(8)[4])) 
dev.off()

####
# ROC for discrimination between different panels based on ATI-N
rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
          "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
rocInput <- 
  heatmap_input_IgG %>% 
  mutate(predictor = if_else(real == "MPXV", "infected", "non-infected"))
roc_ATI <-
  roc(factor(rocInput$predictor, levels = c("infected", "non-infected"), ordered = TRUE), rocInput$`ATI-N`, )
ci_roc_ATI <- ci.coords(roc_ATI, x="best", input = "threshold",  best.method = c("closest.topleft"), ret = rets,
                        best.policy = "random")



####
# Calculate confusion matrix
confusionMatrix <- 
  confusionMatrix(factor(ensemblePrediction$real), factor(ensemblePrediction$max_col))
confusionMatrix[["byClass"]][ , "F1"]


####
# Select only MPXV sera without breakthrough infection
# Select Pre Sera younger than 1983
# Select Pre Sera with HSA signal below 750 MIF
# Calculate fraction of selected sera and assay performance (F1 score)

ensembleSelected <-
  ensemblePrediction %>% 
  filter(!(year_birth < 1972 & real == "Pre"))# %>% 
 # filter(!(real == "Pre" & highBG == TRUE)) %>% 
#  filter((real == "MPXV" & n_MVA_vacc == 0)| real %in% c("Pre", "MVA"))

confusionMatrixSelected <- 
  confusionMatrix(factor(ensembleSelected$real), factor(ensembleSelected$max_col))
confusionMatrixSelected[["byClass"]][ , "F1"]





####
# Calculate performance with threshold of 0.6 on all data
ensembleSelectedThreshold <-
  ensemblePrediction %>% 
 # filter(!(year_birth < 1972 & real == "Pre")) %>% 
#  filter(!(real == "Pre" & highBG == TRUE)) %>% 
  # filter((real == "MPXV" & n_MVA_vacc == 0)| real %in% c("Pre", "MVA")) %>% 
  filter(pred_ensemble > 0.5) %>% 
  filter(mean_mean_conf > 0.5)

confusionMatrixSelectedThreshold <- 
  confusionMatrix(factor(ensembleSelectedThreshold$max_col), factor(ensembleSelectedThreshold $real))
confusionMatrixSelectedThreshold[["byClass"]][ , "F1"]

####
# Export supporting table containing the data

# Select ATI-N results for merge
dataInputJoin <- 
  dataInputComparison %>% 
  select(sampleID_meta, `ATI-N`)

ensemblePredictionExport <-
  ensemblePrediction %>% 
  left_join(dataInputJoin, by = c("sampleID_meta")) %>% 
  mutate(Sample = row_number(),
         mean_mean_conf = round(mean_mean_conf, digits = 2),
         pred_ensemble = round(pred_ensemble, digits = 2),
         background = round(background, digits = 0)) %>% 
  select(Sample, Real = real, Pred = max_col, Pred_LDA = LDA_max_col, 
         Pred_GBC = GBC_max_col,
         `Serostatus (Delta)` = serostatus_cat.delta,
         `Serostatus ATI-N-CPXV` = `ATI-N`,
         `MVA vacc. (n)`= n_MVA_vacc, `Year Birth` = year_birth,
         `Ensemble confidence` = pred_ensemble, 
         `LDA confidence` = mean_mean_conf, `HSA background` = background)

export(ensemblePredictionExport, file = "output/STableEnseblePrediction.xlsx")

ensemblePredictionExport %>% 
  select(Sample, Real, Pred, `Serostatus ATI-N-CPXV`) %>% 
  group_by(Real, Pred, `Serostatus ATI-N-CPXV`) %>% 
  count()


####
# PLot confusion matrices using package cvsm
conf_mat <- confusion_matrix(targets = factor(ensemblePrediction$real),
                             predictions = factor(ensemblePrediction$max_col))

conf_mat_LDA <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                    predictions = factor(ensemblePrediction$LDA_max_col))

conf_mat_GBC <- confusion_matrix(targets = factor(ensemblePrediction$real),
                                           predictions = factor(ensemblePrediction$GBC_max_col))

plotValidation <-
  plot_confusion_matrix(conf_mat,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = TRUE,
                        add_zero_shading = TRUE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges")

plotValidationLDA <-
  plot_confusion_matrix(conf_mat_LDA,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges")
plotValidationGBC <-
  plot_confusion_matrix(conf_mat_GBC,
                        add_sums = FALSE,
                        add_counts = TRUE,
                        add_normalized = FALSE,
                        add_row_percentages = FALSE,
                        add_col_percentages = FALSE,
                        diag_percentages_only = TRUE,
                        rm_zero_percentages = TRUE,
                        rm_zero_text = FALSE,
                        add_zero_shading = FALSE,
                        amount_3d_effect = 1,
                        add_arrows = TRUE,
                        counts_on_top = FALSE,
                        palette = "Oranges")


plotConfMatrixVal <-
  ggarrange(plotValidationLDA, plotValidationGBC, plotValidation, plotPreGrouped, ncol = 4, 
            align = "hv", labels = c("b", "c", "d", "e"))

ggsave("output/plotConfusionVal.pdf", plotConfMatrixVal,
       width = 12, height = 4, dpi = 600)
