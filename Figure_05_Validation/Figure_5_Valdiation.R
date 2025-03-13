####
# Analyse revision measurements
# Daniel Stern
# 2025-03-04
# Version 1.0
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

# Load datainput
load("input/ensemblePrediction.Rdata")
load("input/heatmap_input.Rdata")
load("input/dataInputComparison.Rdata")


####
# Plot ensemble prediction for MPXV after breakthrough infection
plotMVA <-
  ensemblePrediction %>% 
  filter(!is.na(n_MVA_vacc)) %>% 
  filter(real == "MPXV") %>% 
  ggplot(mapping = aes(x = as.factor(n_MVA_vacc), y = (MVA))) +
  geom_boxplot(outliers = F, fill = colorblind_pal()(8)[3]) +
  geom_jitter(width = 0.05, alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_x_discrete(name = "MVA vaccinations (n)") +
  scale_y_continuous(name = "Ensemble prediction MVA") +
  theme_pubr()

plotMPXV <-
  ensemblePrediction %>% 
  filter(!is.na(n_MVA_vacc)) %>% 
  filter(real == "MPXV") %>% 
  ggplot(mapping = aes(x = as.factor(n_MVA_vacc), y = (MPXV))) +
  geom_boxplot(outliers = F, fill = colorblind_pal()(8)[4]) +
  geom_jitter(width = 0.05, alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_x_discrete(name = "MVA vaccinations (n)") +
  scale_y_continuous(name = "Ensemble prediction MPXV") +
  theme_pubr()

# Plot ensemble prediction for MVA after MVA vaccination
plotMVA_MVA <-
  ensemblePrediction %>% 
  filter(!is.na(n_MVA_vacc)) %>% 
  filter(real == "MVA") %>% 
  ggplot(mapping = aes(x = as.factor(n_MVA_vacc), y =(MVA))) +
  geom_boxplot(outliers = F, fill = colorblind_pal()(8)[3]) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  scale_x_discrete(name = "MVA vaccinations (n)") +
  scale_y_continuous(name = "Ensemble prediction MVA",
                     limits = c(0.5, 1)) +
  theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right")

# Plot Figure 5 c
plotMPXV_MVA_Ensemble <-
  ggarrange(plotMVA, plotMPXV, plotMVA_MVA, ncol = 3,align = "hv", 
            common.legend = F)


# Plot Figure 5 d 
plotPreEnsemble <-
  ensemblePrediction %>% 
  mutate(serostatus.delta = factor(serostatus.delta, levels = c(0,1),
                                   labels = c("negative", "positive"))) %>% 
  pivot_longer(cols = c(Pre, MVA, MPXV), names_to = "pred",
               values_to = "ensembl") %>% 
  mutate(pred = factor(pred, levels = c("Pre", "MVA", "MPXV"), 
                       ordered = TRUE)) %>% 
  filter(real == "Pre") %>% 
  ggplot(mapping = aes(x = (year_birth), y = (ensembl), color = pred)) +
  geom_point(aes(shape = serostatus.delta)) +
  labs(shape = "Serostatus",
       colour = "Prediction") +
  scale_shape_manual(values = c(16, 15)) +
  geom_smooth(method = "lm", se = F, linetype = "longdash", alpha = 0.2) +
  scale_color_manual(values = colorblind_pal()(8)[2:8])+
  scale_x_continuous(name = "Year of birth", breaks = c(1940, 1950, 1960, 1970, 1980, 1990, 2000)) +
  scale_y_continuous(name = "Ensemble prediction",
                     breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_pubr()  +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right")


####
# Calcuatie correlation and linear regression to quantify in manuscript
pred_pre <-
  ensemblePrediction %>% 
  filter(real == "Pre") %>% 
  filter(multiple_max == FALSE)

cor.test(pred_pre$mean_mean_conf[pred_pre$max_col == "MVA"], pred_pre$year_birth[pred_pre$max_col == "MVA"])
cor.test(pred_pre$mean_mean_conf[pred_pre$max_col == "Pre"], pred_pre$year_birth[pred_pre$max_col == "Pre"])


# Combine plots
plotMPXVMVAPre <-
  ggarrange(plotMPXV_MVA_Ensemble , plotPreEnsemble, nrow = 2, align = "hv",
            labels = c("c", "d"))

# Save to combine in Illustrator
ggsave("output/plotMPXVMVAPre.pdf", plotMPXVMVAPre, width = 6, height = 6)



####
# Plot heatmap of raw data
set.seed(123)
matIgG <- as.matrix(scale(heatmap_input_IgG[3:17]))
haIgG <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgG$panel),
                           annotation_height = unit(4, "mm"))

pdf(file = "output/heatmap_IgG_white.pdf", width = 10, height = 8)
Heatmap(matIgG, km = 3,
        name = "Multiplex quantfied",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2") +
  Heatmap(factor(heatmap_input_IgG$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Serostatus", col = c(colorblind_pal()(8)[5:6])) +
  Heatmap(factor(heatmap_input_IgG$panel, levels = c("Pre", "MVA", "MPXV"), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4])) +
  Heatmap(heatmap_input_IgG$Pre, name = "Pre ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgG$MVA, name = "MVA ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgG$MPXV, name = "MPXV ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgG$n_MVA_vacc, name = "MVA vaccinations", col = c(viridis_pal()(5)))
dev.off()


matIgM <- as.matrix(scale(heatmap_input_IgM[3:17]))
ha <- HeatmapAnnotation(df = data.frame(panel = heatmap_input_IgM$panel),
                        annotation_height = unit(4, "mm"))
pdf(file = "output/heatmap_IgM_white.pdf", width = 10, height = 8)
Heatmap(matIgM, km = 3,#clustering_distance_rows = "pearson",
        #col = c(viridis_pal()(20)),#c(colorblind_pal()(8)[c(6,5,8)]),
        name = "Multiplex quantfied",
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2") +
  #  Heatmap(heatmap_input_IgM$mean_mean_conf, name = "Confidence", col = c("white", "black")) +
  Heatmap(factor(heatmap_input_IgM$serostatus.delta, levels = c(0,1),
                 labels = c("negative", "positive"),
                 ordered = TRUE), name = "Serostatus", col = c(colorblind_pal()(8)[5:6])) +
  Heatmap(factor(heatmap_input_IgM$panel, levels = c("Pre", "MVA", "MPXV"), 
                 ordered = TRUE), name = "Panel", col = c(colorblind_pal()(8)[2:4])) +
  Heatmap(heatmap_input_IgM$Pre, name = "Pre ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[2])) +
  Heatmap(heatmap_input_IgM$MVA, name = "MVA ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[3])) +
  Heatmap(heatmap_input_IgM$MPXV, name = "MPXV ensemble", na_col = "white", col = c("white", colorblind_pal()(8)[4])) +
  Heatmap(heatmap_input_IgM$n_MVA_vacc, name = "MVA vaccinations", col = c(viridis_pal()(5))) 
dev.off()

####
# ROC for discrimination between different panels based on ATI-N
rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
          "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
rocInput <- 
  heatmap_inputIgG %>% 
  mutate(predictor = if_else(real == "MPXV", "infected", "non-infected"))
roc_ATI <-
  roc(factor(rocInput$predictor, levels = c("infected", "non-infected"), ordered = TRUE), rocInput$`ATI-N`, )
ci_roc_ATI <- ci.coords(roc_ATI, x="best", input = "threshold",  best.method = c("closest.topleft"), ret = rets,
                        best.policy = "random")




# Helper function 
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  dimnames(m)[1]<-x[,1]
  m
}



# Plots Multiplex IgG
dataInputCorIgG <-
  heatmap_input_IgG %>%
  dplyr::select(c(D8:mean_mean_conf)) %>% 
  select(-real, -mean_mean_conf)

dataInputCorIgG[is.na(dataInputCorIgG)] <- 0
# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixneg <- cor_pmat(dataInputCorIgG)
matrixneg <- cor_mat(dataInputCorIgG)

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationAllIgG.pdf"
pdf(height=4, width=4, file=file_path)
corrplot(matrix.please(matrixneg), p.mat = matrix.please(pmatrixneg), method = "circle", order = "hclust",diag = F,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05,
         type = "lower")
dev.off()


####
# Select only MPXV sera without breakthrough infection
# Select Pre Sera younger than 1983
# Select Pre Sera with HSA signal below 750 MIF
# Calculate fraction of selected sera and assay performance (F1 score)

ensembleSelected <-
  ensemblePrediction %>% 
  filter(!(year_birth < 1972 & real == "Pre")) %>% 
  filter(!(real == "Pre" & highBG == TRUE)) %>% 
  filter((real == "MPXV" & n_MVA_vacc == 0)| real %in% c("Pre", "MVA"))

confusionMatrixSelected <- 
  confusionMatrix(factor(ensembleSelected$real), factor(ensembleSelected$max_col))
confusionMatrixSelected[["byClass"]][ , "F1"]


####
# Calculate performance with threshold of 0.6 on all data
ensembleSelectedThreshold <-
  ensemblePrediction %>% 
  #filter(!(year_birth < 1972 & real == "Pre")) %>% 
  #filter(!(real == "Pre" & highBG == TRUE)) %>% 
  # filter((real == "MPXV" & n_MVA_vacc == 0)| real %in% c("Pre", "MVA")) %>% 
  filter(pred_ensemble > 0.6)

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
  select(Sample, Real = real, Pred = max_col,
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

conf_mat_select <- confusion_matrix(targets = factor(ensembleSelected$real),
                                    predictions = factor(ensembleSelected$max_col))

conf_mat_select_thresh <- confusion_matrix(targets = factor(ensembleSelectedThreshold $real),
                                           predictions = factor(ensembleSelectedThreshold$max_col))
plotValidation <-
  plot_confusion_matrix(conf_mat)
plotValidationSelected <-
  plot_confusion_matrix(conf_mat_select)
plotValidationSelectedThresh <-
  plot_confusion_matrix(conf_mat_select_thresh)

plotConfMatrixVal <-
  ggarrange(plotValidation, plotValidationSelected, plotValidationSelectedThresh, ncol = 3, 
            align = "hv", labels = "auto")

ggsave("output/plotConfusionVal.pdf", plotConfMatrixVal,
       width = 8, height = 3, dpi = 600)
