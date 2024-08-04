##
# Correlate NT - Multiplex-Assay
# Daniel Stern
# Robert Koch Institute
# 11.04.2023
# Version 2.0
##

##
# Load libraries
library(rio)
library(tidyverse)
library(pROC)
library(ggpubr)
library(lubridate)
library(boot)
library(purrr)
library(corrplot)
library(rstatix)
library(ggthemes)
library(glmnet)
library(caret)

##
# Clean envrionment
rm(list = ls(all.names = TRUE))

# Source functions
source("../4 Method Comparison/determineCutoffFunction.R", encoding = "UTF-8")
source("../4 Method Comparison/rocFunction.R", encoding = "UTF-8")

##
# Load data
load("../7 Patient panel merge/output/dataInputMPXVmeta.Rdata")
dataInputNT <- import("input/Ergebnisse_NT MPXV_für ZBS3.xlsx",
                      skip = 6) %>% 
  mutate(NT_titer = as.numeric(str_remove(`NT Titer (1:x)`, "<"))) %>% 
  mutate(sampleID_core = str_remove(Probennummer, "-001-[0-9][0-9]"),
         sampleID_core = str_remove(sampleID_core, "P-22-"),
         sampleID_core = str_remove(sampleID_core, "^0+"),
         sampleID_core = str_remove(sampleID_core, "-0")) %>% 
  mutate(last_sampleID = sub('.*(?=.{2}$)', '', Probennummer, perl=T)) 

dataInputNTImvanex <- import("input/2022_10_21_Imvanex_Studie_anonym_Titer für Daniel.xlsx") %>% 
  select(sampleID_metadata = MaterialID, meas_NT = `MVA Titer`) %>% 
  mutate(meas_NT = case_when(meas_NT == "neg" ~ 14,
                             TRUE ~ as.numeric(meas_NT)))

##
# Load color palette for plotting 
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##
# Join Data and select data with NT results only
dataInput <-
  dataInputMPXVmeta %>% 
  mutate(last_sampleID = sub('.*(?=.{2}$)', '', sampleID_metadata, perl=T)) %>% 
  left_join(dataInputNT, by = c("sampleID_core", "last_sampleID")) %>% 
  dplyr::filter(!is.na(NT_titer)) %>% 
  mutate(childhoodImmu = if_else(date_birthday < 1975, 1, 0))

save(dataInput , file = "../4 Method Comparison/input/dataInputNT.R")

plotRegression <- 
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = log2(NT_titer), y = log10(dataIn)) ) +
  geom_point(aes(color = childhoodImmu)) +
  geom_smooth(aes(color = childhoodImmu),method = 'lm', formula = 'y ~ log(x)', se = FALSE) +
  theme_bw() +
  scale_x_continuous(name = "NT Titer (log2)") +
  scale_y_continuous(name = "Analyte (log10 µg/mL VIG)") +
  coord_cartesian(ylim = c(0.5,5.5)) +
  scale_color_manual(name = "Year of birth", values = colorBlindBlack8[2:3]) +
  facet_wrap("analyte", scales = "free_y") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "top")

ggsave("output/plotRegression.png", plotRegression, width = 8, height = 8, dpi = 300)
##
# Plot correlation
## Plot Korrelation zwischen den verschiedenen Assays
####
# Plot Correlation for normalized results for different antigens
# Function to transform dataframe (output from cor_pmat and cor_mat) to matrix
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  dimnames(m)[1]<-x[,1]
  m
}


# Plots Correlation between multiplex and NT IgG
dataInputCorIgG <-
  dataInput %>%
  dplyr::filter(dilution_assay == 100
                & isotype == "IgG" &
                  assaytype == "Multiplex" &
                  #                  childhoodImmu == 1 &
                  analyte != "VACV") %>%
  dplyr::select(sampleID_metadata, analyte, dataIn, NT_titer) %>%
  unique() %>%
  mutate(dataIn = log10(dataIn),
         NT_titer = log2(NT_titer)) %>% 
  pivot_wider(names_from = "analyte", values_from = "dataIn") %>%
  dplyr::select(-sampleID_metadata)

# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixneg <- cor_pmat(dataInputCorIgG, method = "spearman")
matrixneg <- cor_mat(dataInputCorIgG, method = "spearman")

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationNTIgG.pdf"
pdf(height=10, width=10, file=file_path)
corrplot(matrix.please(matrixneg), p.mat = matrix.please(pmatrixneg), method = "square", order = "hclust",diag = T,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05)
dev.off()

dataInputCorIgM <-
  dataInput %>%
  dplyr::filter(dilution_assay == 100
                & isotype == "IgM" &
                  assaytype == "Multiplex" &
                  analyte != "VACV") %>%
  dplyr::select(sampleID_metadata, analyte, dataIn, NT_titer) %>%
  unique() %>%
  mutate(dataIn = log10(dataIn),
         NT_titer = log2(NT_titer)) %>% 
  pivot_wider(names_from = "analyte", values_from = "dataIn") %>%
  dplyr::select(-sampleID_metadata)

# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixnegNI <- cor_pmat(dataInputCorIgM, method = "spearman")
matrixnegNI <- cor_mat(dataInputCorIgM, method = "spearman")

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationIgM.pdf"
pdf(height=10, width=10, file=file_path)
corrplot(matrix.please(matrixnegNI), p.mat = matrix.please(pmatrixnegNI), method = "square", order = "hclust",diag = T,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05)
dev.off()


plotCorrelationNT <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(!(analyte %in% c("VACV", "Delta"))) %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29","L1R", "M1",
                                              "D8L", "E8",
                                              "H3L", 
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ paste0(intToUtf8(8804), "14"),
                                  NT_titer <= 25 ~ paste0(intToUtf8(8804), "25"),
                                  NT_titer <= 50 ~ paste0(intToUtf8(8804), "50"),
                                  NT_titer <= 75 ~ paste0(intToUtf8(8804), "75"),
                                  NT_titer <= 100 ~ paste0(intToUtf8(8804), "100")),
         NT_titer_bin = factor(NT_titer_bin,
                               levels = c(paste0(intToUtf8(8804), "14"),
                                          paste0(intToUtf8(8804), "25"), 
                                          paste0(intToUtf8(8804), "50"),
                                          paste0(intToUtf8(8804), "75"),
                                          paste0(intToUtf8(8804), "100")), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = (NT_titer_bin), y = log10(dataIn), fill = isotype) ) +
  geom_boxplot()+
  geom_point() +
  theme_bw() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  #  coord_cartesian(ylim = c(0.5,5.5)) +
  scale_fill_manual(name = "Year of birth", values = colorBlindBlack8[4]) +
  facet_wrap("analyte", scales = "free_y") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")

plotCorrelationNTIgM <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgM") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(!(analyte %in% c("VACV", "Delta"))) %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29","L1R", "M1",
                                              "D8L", "E8",
                                              "H3L", 
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ paste0(intToUtf8(8804), "14"),
                                  NT_titer <= 25 ~ paste0(intToUtf8(8804), "25"),
                                  NT_titer <= 50 ~ paste0(intToUtf8(8804), "50"),
                                  NT_titer <= 75 ~ paste0(intToUtf8(8804), "75"),
                                  NT_titer <= 100 ~ paste0(intToUtf8(8804), "100")),
         NT_titer_bin = factor(NT_titer_bin,
                               levels = c(paste0(intToUtf8(8804), "14"),
                                          paste0(intToUtf8(8804), "25"), 
                                          paste0(intToUtf8(8804), "50"),
                                          paste0(intToUtf8(8804), "75"),
                                          paste0(intToUtf8(8804), "100")), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = (NT_titer_bin), y = log10(dataIn), fill = isotype) ) +
  geom_boxplot()+
  geom_point() +
  theme_bw() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  #  coord_cartesian(ylim = c(0.5,5.5)) +
  scale_fill_manual(name = "Year of birth", values = colorBlindBlack8[4]) +
  facet_wrap("analyte", scales = "free_y") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")



plotCorrelationNTDelta <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte == "Delta") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29","L1R", "M1",
                                              "D8L", "E8",
                                              "H3L", 
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ paste0(intToUtf8(8804), "14"),
                                  NT_titer <= 25 ~ paste0(intToUtf8(8804), "25"),
                                  NT_titer <= 50 ~ paste0(intToUtf8(8804), "50"),
                                  NT_titer <= 75 ~ paste0(intToUtf8(8804), "75"),
                                  NT_titer <= 100 ~ paste0(intToUtf8(8804), "100")),
         NT_titer_bin = factor(NT_titer_bin,
                               levels = c(paste0(intToUtf8(8804), "14"),
                                          paste0(intToUtf8(8804), "25"), 
                                          paste0(intToUtf8(8804), "50"),
                                          paste0(intToUtf8(8804), "75"),
                                          paste0(intToUtf8(8804), "100")), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = (NT_titer_bin), y = log10(dataIn), fill = isotype) ) +
  geom_boxplot()+
  geom_point() +
  theme_bw() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  #  coord_cartesian(ylim = c(0.5,5.5)) +
  scale_fill_manual(name = "Year of birth", values = colorBlindBlack8[4]) +
  # facet_wrap("analyte", scales = "free_y") +
  labs(subtitle = "IgG") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")

plotCorrelationNTDeltaIgM <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgM") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte == "Delta") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29","L1R", "M1",
                                              "D8L", "E8",
                                              "H3L", 
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ paste0(intToUtf8(8804), "14"),
                                  NT_titer <= 25 ~ paste0(intToUtf8(8804), "25"),
                                  NT_titer <= 50 ~ paste0(intToUtf8(8804), "50"),
                                  NT_titer <= 75 ~ paste0(intToUtf8(8804), "75"),
                                  NT_titer <= 100 ~ paste0(intToUtf8(8804), "100")),
         NT_titer_bin = factor(NT_titer_bin,
                               levels = c(paste0(intToUtf8(8804), "14"),
                                          paste0(intToUtf8(8804), "25"), 
                                          paste0(intToUtf8(8804), "50"),
                                          paste0(intToUtf8(8804), "75"),
                                          paste0(intToUtf8(8804), "100")), ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = (NT_titer_bin), y = log10(dataIn), fill = isotype) ) +
  geom_boxplot()+
  geom_point() +
  theme_bw() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  #  coord_cartesian(ylim = c(0.5,5.5)) +
  scale_fill_manual(name = "Year of birth", values = colorBlindBlack8[4]) +
  # facet_wrap("analyte", scales = "free_y") +
  labs(subtitle = "IgM") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")


##
# Calculate spearman rank correlation between NT-Titer and log10(Analyte) multiplex
correlation <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  #  dplyr::filter(analyte == "Delta") %>%
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ 1,
                                  NT_titer <= 25 ~ 2,
                                  NT_titer <= 50 ~ 3,
                                  NT_titer <= 75 ~ 4,
                                  NT_titer <= 100 ~ 5)) %>% 
  group_by(analyte) %>% 
  summarise(correlation = list(cor.test(NT_titer_bin, log10(dataIn), method = 'spearman')),
            pearson = list(cor.test(NT_titer, log10(dataIn), method = 'pearson')))
#corr <- cor.test(as.numeric(levels(x), y, method = 'spearman')


p <- vector()
estimate <- vector()
for(i in 1:nrow(correlation)){
  p[i] <-correlation$correlation[[i]]$p.value
  estimate[i] <-correlation$correlation[[i]]$estimate
}

p_p <- vector()
estimate_p <- vector()
for(i in 1:nrow(correlation)){
  p_p[i] <-correlation$pearson[[i]]$p.value
  estimate_p[i] <-correlation$pearson[[i]]$estimate
}


correlationNT <-
  data.frame(analyte = correlation$analyte,
             rho_estimate = estimate,
             p_estimate = p,
             rho_estimate_pearson = estimate_p,
             p_estimate_pearson = p_p)



##
# Calculate spearman rank correlation between NT-Titer and log10(Analyte) multiplex
# for IgM results
correlationIgM <-
  dataInput %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgM") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  #  dplyr::filter(analyte == "Delta") %>%
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  mutate(NT_titer_bin = case_when(NT_titer == 14 ~ 1,
                                  NT_titer <= 25 ~ 2,
                                  NT_titer <= 50 ~ 3,
                                  NT_titer <= 75 ~ 4,
                                  NT_titer <= 100 ~ 5)) %>% 
  group_by(analyte) %>% 
  summarise(correlation = list(cor.test(NT_titer_bin, log10(dataIn), method = 'spearman')),
            pearson = list(cor.test(NT_titer, log10(dataIn), method = 'pearson')))
#corr <- cor.test(as.numeric(levels(x), y, method = 'spearman')


p_IgM <- vector()
estimate_IgM <- vector()
for(i in 1:nrow(correlationIgM)){
  p_IgM[i] <-correlationIgM$correlation[[i]]$p.value
  estimate_IgM[i] <-correlationIgM$correlation[[i]]$estimate
}

p_p_IgM <- vector()
estimate_p_IgM <- vector()
for(i in 1:nrow(correlationIgM)){
  p_p_IgM[i] <-correlationIgM$pearson[[i]]$p.value
  estimate_p_IgM[i] <-correlationIgM$pearson[[i]]$estimate
}


correlationNTIgM <-
  data.frame(analyte = correlationIgM$analyte,
             rho_estimate = estimate_IgM,
             p_estimate = p_IgM,
             rho_estimate_pearson = estimate_p_IgM,
             p_estimate_pearson = p_p_IgM)



##
# Export results for import in methodcomparison
save(plotCorrelationNT, correlationNT, plotCorrelationNTDelta, 
     plotCorrelationNTIgM, plotCorrelationNTDeltaIgM, correlationNTIgM,
     file = "../4 Method Comparison/input/corrNT.Rdata")


##
# Perform ROC in comparison to NT assay
## Perform ROC Analysis in comparison to IFA
# Select data and transform to wide format for IgG
rocDataframeNT <- 
  dataInput %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, panel, NT_titer, assaytype) %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  filter(!is.na(NT_titer)) %>% 
  dplyr::select(-dilution_assay, -isotype) %>%
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  filter(assaytype == "Multiplex") %>% 
  mutate(analyte = str_replace(analyte, "-", ".")) %>%
  mutate(serostatus = if_else(NT_titer == 14, 0, 1)) %>%
  #group_by(sampleID_metadata, panel, titer, serostatus, analyte) %>% 
  #summarise(meandataIn = mean(dataIn, na.rm = TRUE)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn) %>% 
  # filter(!is.na(A27L)) %>% 
  ungroup()

# Perform ROC for IgG
outplot_listNT <- list()
j = 1
for(i in c("A27L", "A29",
           "D8L", "E8",
           "H3L", "L1R", "M1",
           "A33R", "A35R",
           "B5R", "B6",
           "A5L", "ATI.C", "ATI.N", "Delta")){
  outplot_listNT[[j]] <- rocFunction(rocDataframeNT, "serostatus",i)
  assign(paste("roc",  i, sep = "_"), outplot_listNT[[j]])
  j = j+1
}

plot(roc_E8)
# Combined ROC for all multiplex-antigens by fitting a generalized linear model

aucCombinedNT <-glm(serostatus ~ 
                      A27L + 
                      D8L +  
                      H3L + 
                      L1R +
                      A33R +
                      B5R +
                      A5L +
                      ATI.C +
                      ATI.N +
                      Delta,
                    data=rocDataframeNT, family="binomial")  #Logistic regression model
roc_CombinedNT <- roc(serostatus ~ predict(aucCombinedNT), data=rocDataframeNT)

summary(aucCombinedNT)

plot(roc_CombinedNT)

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
paramDataFrameNT <- rocParamFunc(param_list)
paramDataFrame_transpose_sep <- transformRocParamFunc(param_list, "NT")
sensNT <- prepDataFrameFunc(paramDataFrame_transpose_sep, "sensitivity", "ci_roc_") %>%
  filter(antigene != "CombinedNT") %>% 
  mutate(isotype = "IgG") %>%
  mutate(parameter = "sensitivity") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype)
specNT <- prepDataFrameFunc(paramDataFrame_transpose_sep, "specificity", "ci_roc_") %>%
  filter(antigene != "CombinedNT") %>% 
  mutate(isotype = "IgG") %>%
  mutate(parameter = "specificity") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype)
accNT <- prepDataFrameFunc(paramDataFrame_transpose_sep, "accuracy", "ci_roc_") %>%
  filter(antigene != "CombinedNT") %>% 
  mutate(isotype = "IgG") %>%
  mutate(parameter = "accuracy") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype)

##
# Export results for method comparison NT test
export(sensNT, file = "output/sensNT.xlsx")
export(specNT, file = "output/specNT.xlsx")
export(accNT, file = "output/accNT.xlsx")


##
# Plot ROC 
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

plotrocA27L <- 
  plotROCfunction(list("A27L" = `roc_A27L`,
                       "A29" = `roc_A29`), "IgG")
plotrocA33R <- 
  plotROCfunction(list("A33R" = `roc_A33R`,
                       "A35R" = `roc_A35R`), "IgG")
plotrocB5 <- 
  plotROCfunction(list("B5R" = `roc_B5R`,
                       "B6" = `roc_B6`), "IgG")
plotrocA5L <- 
  plotROCfunction(list("A5L" = `roc_A5L`,
                       "ATI-C" = `roc_ATI.C`,
                       "ATI-N" = `roc_ATI.N`), "IgG")
plotrocL1R <- 
  plotROCfunction(list("L1R" = `roc_L1R`,
                       "M1" = `roc_M1`), "IgG")
plotrocD8L <- 
  plotROCfunction(list("D8L" = `roc_D8L`,
                       "E8" = `roc_E8`,
                       "H3L" = `roc_H3L`), "IgG")
plotrocLysate <- 
  plotROCfunction(list("Delta" = `roc_Delta`), "IgG")

plotROCNT <-
  ggarrange(plotrocLysate,
            plotrocA27L,
            plotrocL1R,
            plotrocD8L,
            plotrocA33R,
            plotrocB5,
            plotrocA5L, 
            ncol = 4, nrow = 2,
            align = "hv")

ggsave("output/plotROCNT.png", plotROCNT, width = 10, height = 7, dpi = 300)


##
# Perform lasso regression for NT dataset
# https://www.statology.org/lasso-regression-in-r/
# https://www.datacareer.ch/blog/ridge-and-lasso-in-r/

set.seed(123)

# Define response variable
y <- log10(rocDataframeNT$NT_titer)

# Define matrix of predictor variables
X <- 
  dataInput %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, panel, NT_titer, assaytype) %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  filter(!is.na(NT_titer)) %>% 
  dplyr::select(-dilution_assay, -isotype) %>%
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  filter(assaytype == "Multiplex") %>% 
  mutate(analyte = str_replace(analyte, "-", ".")) %>%
  mutate(serostatus = if_else(NT_titer == 14, 0, 1)) %>%
  mutate(dataIn = log10(dataIn)) %>% 
  #group_by(sampleID_metadata, panel, titer, serostatus, analyte) %>% 
  #summarise(meandataIn = mean(dataIn, na.rm = TRUE)) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn) %>% 
  # filter(!is.na(A27L)) %>% 
  ungroup() %>% 
  select(A29:Delta) %>% 
  as.matrix()

##
# Split dataset
index <- createDataPartition(y, p=0.7, list=FALSE)
X_train <- X[ index, ]
X_test <- X[-index, ]
y_train <- y[index]
y_test<-y[-index]

#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(X_train, y_train, alpha = 1)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda

#produce plot of test MSE by lambda value
plot(cv_model) 

#find coefficients of best model
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
coef(best_model)

#use fitted best model to make predictions
y_predicted <- predict(best_model, s = best_lambda, newx = X)


#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)


#find R-Squared
rsq <- 1 - sse/sst
rsq


## 
# Save model
save(best_model, best_lambda, file = "output/NTLassoReg.Rdata")


##
# Load patient data
load("input/lassoData.Rdata")
dataInputLassoMPXV <-
  dataInputMPXVexport %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  mutate(analyte = str_replace(analyte, "-", "."), 
         dataIn = log10(dataIn)) %>%
  select(sampleID_metadata, caseID_metadata, time_class, childhoodImmu, 
         vaccination_pox, analyte, dataIn) %>% 
  unique() %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")

dataInputLassoMVA <-
  dataInputMVAexport %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  mutate(analyte = str_replace(analyte, "-", "."), 
         dataIn = log10(dataIn)) %>%
  select(sampleID_metadata, caseID_metadata, childhoodImmu, 
         serostatus_factor, analyte, dataIn) %>% 
  unique() %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = "mean")


##
# Predict NT titer in Mpox panel based on lasso regression
mpox_pred_in <- 
  dataInputLassoMPXV %>% 
  select(A29:Delta) %>% 
  as.matrix()

predNT <-
  predict(best_model, s = best_lambda, newx = mpox_pred_in) %>% 
  as.data.frame()

dataPredMPXV <-
  dataInputLassoMPXV %>% 
  select(sampleID_metadata) %>% 
  mutate(pred_NT = predNT$s1)

##
# Predict NT titer in MVA panel based on lasso regression
mva_pred_in <- dataInputLassoMVA %>% 
  select(A29:Delta) %>% 
  as.matrix()

predNTMVA<-
  predict(best_model, s = best_lambda, newx = mva_pred_in) %>% 
  as.data.frame()

dataPredMVA <-
  dataInputLassoMVA %>% 
  select(sampleID_metadata) %>% 
  mutate(pred_NT = predNTMVA$s1)

##
# Compare input with prediction
plotLassoComparison <-
  dataInput %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  select(sampleID_metadata, analyte, NT_titer, childhoodImmu) %>% 
  unique() %>% 
  left_join(dataPredMPXV, by = c("sampleID_metadata")) %>% 
  ggplot(mapping = aes(x = log10(NT_titer), y = (pred_NT))) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  scale_x_continuous(name = "NT (Titer)") +
  scale_y_continuous(name = "NT (Titer) Predicted") +
  theme(strip.background = element_blank(),
        legend.position = "none")


data_pred_NT_MPXV <-
  dataInputMPXVexport %>% 
  left_join(dataPredMPXV, by = c("sampleID_metadata"))

plotDeltaNTPred <-
  data_pred_NT_MPXV %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  #  filter(!(analyte %in% c("Delta"))) %>%
  dplyr::filter(analyte != "VACV") %>% 
  #  dplyr::filter(analyte == "Delta") %>%
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  ggplot(mapping = aes(x = log10(dataIn), y = (pred_NT), color = "Panel")) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  facet_wrap("analyte", scales = "free_x") +
  scale_color_manual(values = colorBlindBlack8[4]) +
  scale_x_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  scale_y_continuous(name = "NT (Titer) Predicted", limits = c(0.9, 2)) +
  theme(strip.background = element_blank(),
        legend.position = "none")


##
# Plot MVA panel
data_pred_NT_MVA <-
  dataInputMVAexport %>% 
  left_join(dataPredMVA, by = c("sampleID_metadata"))

plotDeltaNTPredMVA <-
  data_pred_NT_MVA %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  #  filter(!(analyte %in% c("Delta"))) %>%
  dplyr::filter(analyte != "VACV") %>% 
  #  dplyr::filter(analyte == "Delta") %>%
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  ggplot(mapping = aes(x = log10(dataIn), y = (pred_NT), color = "Panel")) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  scale_color_manual(values = colorBlindBlack8[2]) +
  facet_wrap("analyte", scales = "free_x") +
  scale_x_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  scale_y_continuous(name = "NT (Titer) Predicted", limits = c(0.9, 2)) +
  theme(strip.background = element_blank(),
        legend.position = "none")

plot_pred_comparison <- 
  ggarrange(plotDeltaNTPredMVA, plotDeltaNTPred, ncol = 2, align = "hv")


ggsave("output/plot_pred_comparison.png", plot_pred_comparison, width = 14, height = 9, dpi = 300)

##
# Compare seroconversion rate
data_pred_1 <-
  data_pred_NT_MVA %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  select(sampleID_metadata, caseID_metadata, childhoodImmu, 
         serostatus = serostatus_factor, pred_NT) %>% 
  unique() %>% 
  mutate(panel = "Imvanex")

data_pred_2 <-
  data_pred_NT_MPXV %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  dplyr::filter(analyte != "VACV") %>% 
  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
                                              "D8L", "E8",
                                              "H3L", "L1R", "M1",
                                              "A33R", "A35R",
                                              "B5R", "B6",
                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  select(sampleID_metadata, caseID_metadata, childhoodImmu, 
         serostatus = time_class, pred_NT) %>% 
  unique() %>% 
  mutate(panel = "Mpox")

data_pred <-
  rbind(data_pred_1, data_pred_2) %>% 
  filter(!is.na(serostatus)) %>% 
  filter(!is.na(childhoodImmu)) %>% 
  mutate(NT_status = if_else(pred_NT <= log10(14), 0, 1))

data_pred %>% 
  ggplot(mapping = aes(x = serostatus, y = pred_NT, fill = panel)) +
  geom_boxplot() +
  geom_point(aes(color = as.factor(NT_status))) +
  facet_wrap("childhoodImmu") +
  theme_bw() +
  scale_color_manual(values = c(colorBlindBlack8[5], colorBlindBlack8[3]) ) +
  scale_fill_manual(values = c(colorBlindBlack8[2], colorBlindBlack8[4])) +
  #  scale_x_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  #  scale_y_continuous(name = "NT (Titer) Predicted", limits = c(0.9, 2)) +
  theme(strip.background = element_blank(),
        legend.position = "none")


##
# Compare predicted NT-titer with measured titers for Imvanex sera
plotNTMeasPrep <-
  data_pred %>% 
  left_join(dataInputNTImvanex, by = "sampleID_metadata") %>% 
  filter(!is.na(meas_NT)) %>% 
  filter(!(sampleID_metadata %in% c("S-20-004-005-02", "S-20-004-005-03"))) %>%  # Exclude samples that have been mixed up
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("\u2265 1975", "< 1975"))) %>% 
  mutate(pred_NT = if_else(pred_NT <= log10(14), log10(14), pred_NT)) %>% 
  ggplot(mapping = aes(y = pred_NT, x = log10(meas_NT))) +
  geom_point(aes(color = as.factor(childhoodImmu))) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  # facet_wrap("childhoodImmu") +
  theme_bw() +
  scale_x_continuous(name = "NT Titer measured (MVA)") +
  scale_y_continuous(name = "NT Titer predicted (MPox)") +
  scale_color_manual(values = c(colorBlindBlack8[2], colorBlindBlack8[3]) ) +
  theme(strip.background = element_blank(),
        legend.position = "none")

ggsave("output/plotNTMeasPrep.png", plotNTMeasPrep, width = 4, height = 4, dpi = 300)

dataNTcomparison <-
  data_pred %>% 
  left_join(dataInputNTImvanex, by = "sampleID_metadata") %>% 
  filter(!is.na(meas_NT)) 

correlationNT <-
  data_pred %>% 
  left_join(dataInputNTImvanex, by = "sampleID_metadata") %>% 
  filter(!is.na(meas_NT)) %>% 
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("\u2265 1975", "< 1975"))) %>% 
  mutate(pred_NT = if_else(pred_NT <= log10(14), log10(14), pred_NT)) %>% 
  group_by(childhoodImmu) %>% 
  summarise(correlation = list(cor.test(pred_NT, log10(meas_NT), method = 'pearson')))

p <- vector()
estimate <- vector()
for(i in 1:nrow(correlationNT)){
  p[i] <-correlationNT$correlation[[i]]$p.value
  estimate[i] <-correlationNT$correlation[[i]]$estimate
}


correlationNTdataframe <-
  data.frame(childhoodImmu = correlationNT$childhoodImmu,
             rho =  round(estimate, digits = 2),
             p =  round(p, digits = 5)) %>% 
  arrange(desc(rho))

export(correlationIFAdataframe, file = ("output/correlationIFA.xlsx"))


##
# Export correlation data
NT_prediction_measured <-
  data_pred %>% 
  left_join(dataInputNTImvanex, by = "sampleID_metadata") %>% 
  filter(!is.na(meas_NT)) %>% 
  mutate(pred_NT = 10^pred_NT) %>% 
  arrange(childhoodImmu, meas_NT)

export(NT_prediction_measured, file = "output/NT_prediction_measured.xlsx")

discrepantSamplesLow <-
  data_pred_NT_MVA %>% 
  filter(sampleID_metadata == "S-20-004-005-03") %>% 
  filter(dilution_assay == 100)

