##
# Method comparison for establishent of Mpox-Multiplex-Assays
# Daniel Stern
# 2023/05/22
# Version 1.1
##

##
# This script analyses the assay performance of the Mpox-Multiplex-Assay
# To this aim, several method comparisons are performed:
# 1) Method comparison with a VACV-lysate based ELISA
# - To this aim, linear regression between the quantified results is performce
# - Both IgG and IgM detection is compared
# 2) Method comparison with IFA
# - Based on a IFA-Cutoff of 1:80, ROC analysis is performed for all analytes
# - Both IgG and IgM detection is compared
##

## 
# Clean environment
rm(list = ls(all.names = TRUE))

## 
# Load packages
library(rio)
library(tidyverse)
library(pROC)
library(mcr)
library(ggpubr)
library(ggthemes)

# Source functions
source("determineCutoffFunction.R", encoding = "UTF-8")
source("rocFunction.R", encoding = "UTF-8")

# Load data and metadata
load("../2 Import metadata/output/metadata_IFA.Rdata")
load("../3 Quantify Normalize data/output/dataInputQuant.Rdata")
load("../3 Quantify Normalize data/output/dataInputQuantNK.Rdata")
load("input/corrNT.Rdata")
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
metadata_IFA <-
  metadata_IFA %>% 
  mutate(sampleID_metadata_mut = if_else(panel == "MPXV",
                                         str_replace(sampleID_metadata, 
                                                     "P-22-", "P-22-00"),
                                         sampleID_metadata))
# Generate new variable dataIn
dataInputQuant <-
  dataInputQuant %>% 
  mutate(dataIn = case_when(assaytype == "Multiplex" ~ quant_conc_fix,
                            assaytype == "ELISA" ~ quant_conc_fix)) 

##
# 1) Method comparison with ELISA
# Select only data for which both ELISA and Multiplex-data is available
samples_multiplex <-
  dataInputQuant %>% 
  filter(assaytype == "Multiplex") %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()

samples_ELISA <-
  dataInputQuant %>% 
  filter(assaytype == "ELISA") %>% 
  select(sampleID_metadata) %>% 
  unique() %>% 
  pull()


comp_panel <-
  intersect(samples_ELISA, samples_multiplex)

comp_panel <-
  comp_panel[2:length(comp_panel)] # Exclude NC-Sample

rm(samples_ELISA, samples_multiplex)

# Prepare dataframe for plotting
dataInputMultiplex <-
  dataInputQuant %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, remarks, panel) %>% 
  filter(sampleID_metadata %in% comp_panel) %>% 
  filter(dilution_assay == 100) %>% 
  filter(assaytype == "Multiplex") %>% 
  group_by(analyte) %>% 
  select(-assaytype, -dilution_assay)

dataInputELISA <-
  dataInputQuant %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, remarks) %>% 
  filter(sampleID_metadata %in% comp_panel) %>% 
  filter(analyte == "Delta") %>% 
  filter(dilution_assay == 100) %>% 
  filter(assaytype == "ELISA") %>% 
  group_by(analyte) %>% 
  select(-assaytype, -dilution_assay)

dataInputComparison <-
  dataInputMultiplex %>% 
  left_join(dataInputELISA, by = c("isotype", "sampleID_metadata"), suffix = c("", "_ELISA"))

rm(dataInputELISA, dataInputMultiplex)

# Select IgG and Delta for comparison and exclude deviant samples which 
# must be reanalysed
dataDeltaIgG <-
  dataInputComparison %>% 
  filter(isotype == "IgG" & analyte == "Delta") %>% 
  filter(sampleID_metadata != "P-22-01381-001-18") 

# Perform Passing Bablok Regression for IgG data
pbBaRegDeltaIgG <- mcreg((log10(dataDeltaIgG$dataIn_ELISA)), 
                         log10(dataDeltaIgG$dataIn), method.reg = "PaBa")

plot(pbBaRegDeltaIgG)

# Select IgM and Delta for comparison and exclude deviant samples which 
# must be reanalysed
dataDeltaIgM <-
  dataInputComparison %>% 
  filter(isotype == "IgM" & analyte == "Delta") %>% 
  filter(sampleID_metadata != "P-22-01381-001-18") %>% 
  filter(sampleID_metadata !=	"P-22-01133-001-14")

# Perform Passing Bablok Regression for IgM data
pbBaRegDeltaIgM <- mcreg(log10(dataDeltaIgM$dataIn_ELISA), 
                         log10(dataDeltaIgM$dataIn), method.reg = "PaBa")

plot(pbBaRegDeltaIgM)

# Write parameter for export in dataframe
parampaBa <- data.frame(isotype = c("IgG", "IgM"),
                        pearsonsR = c(0.933, 0.846),
                        slope = c(0.87, 0.89),
                        intercept = c(0.13, -0.63),
                        n = c(75, 75))

export(parampaBa, "output/parampaBa.xlsx")

## 
# Generate ggplots for Passing Bablok Regression Fis
plotpaBaIgG <-
  dataDeltaIgG %>% 
  ggplot(aes(x = log10(dataIn_ELISA), y = log10(dataIn)))+ 
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = colorBlindBlack8[4])+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,1], intercept = pbBaRegDeltaIgG@para[1,1], color = colorBlindBlack8[3], size = 1, alpha = 0.5)+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,3], intercept = pbBaRegDeltaIgG@para[1,3], lty = "dashed", color = colorBlindBlack8[2])+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,4], intercept = pbBaRegDeltaIgG@para[1,4], lty = "dashed", color = colorBlindBlack8[2])+ 
  scale_x_continuous(name = "ELISA (log10 VIG µg/mL)", limits = c(0,5)) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,5)) +
  labs(subtitle = "IgG") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5)) 

plotpaBaIgM <-
  dataDeltaIgM %>% 
  ggplot(aes(x = log10(dataIn_ELISA), y = log10(dataIn)))+ 
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = colorBlindBlack8[4])+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,1], intercept = pbBaRegDeltaIgM@para[1,1], color = colorBlindBlack8[3], size = 1, alpha = 0.5)+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,3], intercept = pbBaRegDeltaIgM@para[1,3], lty = "dashed", color = colorBlindBlack8[2])+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,4], intercept = pbBaRegDeltaIgM@para[1,4], lty = "dashed", color = colorBlindBlack8[2])+ 
  scale_x_continuous(name = "ELISA (log10 VIG µg/mL)", limits = c(0,5)) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,5)) +
  labs(subtitle = "IgM") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5))

plotpaBacombined <-
  ggarrange(plotpaBaIgG, plotpaBaIgM, ncol = 2, align = "hv")

##
# 2) Method comparison with IFA
plotIFAIgG <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata" = "sampleID_metadata_mut", "isotype", "panel")) %>% 
  select(assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(analyte == "Delta") %>% 
  filter(isotype == "IgG") %>% 
  filter(assaytype == "Multiplex") %>% 
  mutate(titer = if_else(titer < 320, (paste0(intToUtf8(8804), "80")), as.character(titer)),
         titer = factor(titer, levels = c(paste0(intToUtf8(8804), "80"), 
                                          "320", "1280", "5120", "20480"),
                        ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = as.factor((titer)), y = log10(dataIn), fill = assaytype)) +
  geom_boxplot(notch = F, alpha = 1, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), alpha = 1, aes(group = assaytype), color = "black") +
  theme_bw()+
  scale_fill_manual(name = "Assay", values = colorBlindBlack8[2:3]) + 
  theme(strip.background = element_blank()) +
  xlab("IFA (Titer)")+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,5)) +
  labs(subtitle = "IgG") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

plotIFAIgM <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata", "isotype", "panel")) %>% 
  select(assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(analyte == "Delta") %>% 
  filter(isotype == "IgM") %>% 
  filter(assaytype == "Multiplex") %>% 
  mutate(titer = if_else(titer < 320, (paste0(intToUtf8(8804), "80")), as.character(titer)),
         titer = factor(titer, levels = c(paste0(intToUtf8(8804), "80"), 
                                          "320", "1280", "5120", "20480"),
                        ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = as.factor((titer)), y = log10(dataIn), fill = assaytype)) +
  geom_boxplot(notch = F, alpha = 1, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), alpha = 1, aes(group = assaytype), color = "black") +
  theme_bw()+
  scale_fill_manual(name = "Assay", values = colorBlindBlack8[3], drop = F) + 
  theme(strip.background = element_blank()) +
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,5)) +
  labs(subtitle = "IgM") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(name = "IFA (Titer)")

plotIFAcombined <-
  ggarrange(plotIFAIgG , plotIFAIgM, plotCorrelationNTDelta, ncol = 3, 
            heights = c(4,4,4), align = "hv", common.legend = F)


plotMultiplexIgG <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata", "isotype", "panel")) %>% 
  select(assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(!(analyte %in% c("Delta", "VACV"))) %>% 
  filter(isotype == "IgG") %>% 
  mutate(titer = if_else(titer < 320, (paste0(intToUtf8(8804), "80")), as.character(titer)),
         titer = factor(titer, levels = c(paste0(intToUtf8(8804), "80"), 
                                          "320", "1280", "5120", "20480"),
                        ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = as.factor((titer)), y = log10(dataIn), fill = isotype)) +
  geom_boxplot(notch = F, alpha = 1, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), alpha = 1, aes(group = isotype), color = "black") +
  theme_bw()+
  scale_fill_manual(name = "Isotype", values = colorBlindBlack8[2]) + 
  facet_wrap("analyte", scales = "free_y") +
  xlab("IFA (Titer)")+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")

plotMultiplexIgM <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata", "isotype", "panel")) %>% 
  select(assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(!(analyte %in% c("Delta", "VACV"))) %>% 
  filter(isotype == "IgM") %>% 
  mutate(titer = if_else(titer < 320, (paste0(intToUtf8(8804), "80")), as.character(titer)),
         titer = factor(titer, levels = c(paste0(intToUtf8(8804), "80"), 
                                          "320", "1280", "5120", "20480"),
                        ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = as.factor((titer)), y = log10(dataIn), fill = isotype)) +
  geom_boxplot(notch = F, alpha = 1, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), alpha = 1, aes(group = isotype), color = "black") +
  theme_bw()+
  scale_fill_manual(name = "Isotype", values = colorBlindBlack8[3]) + 
  facet_wrap("analyte", scales = "free_y") +
  xlab("IFA (Titer)")+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  theme_bw() +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")


plotDelta <-
  ggarrange(plotpaBacombined, plotIFAcombined, nrow = 2, align = "hv", labels = "auto",
            widths = c(2,3))
##
# Export plots
ggsave("output/plotDelta_Fig_2.jpg", plotDelta, width = 8, height = 6, dpi = 600)


##
# Export plots supporting figures
plotMethCompSupp <- 
  ggarrange(plotMultiplexIgG, plotMultiplexIgM, plotCorrelationNT, nrow = 3,
            labels = "auto")

ggsave("output/plotMethCompSupp.jpg", plotMethCompSupp, width = 8, height = 15, dpi = 600)
ggsave("output/SFig_Multi_IgG.jpg", plotMultiplexIgG, width = 8, height = 6,
       dpi = 600)
ggsave("output/SFig_Multi_IgM.jpg", plotMultiplexIgM, width = 8, height = 6,
       dpi = 600)
ggsave("output/SFig_Multi_NT.jpg", plotCorrelationNT, width = 8, height = 6,
       dpi = 600)

ggsave("output/SFig_Multi_NT_IgM.jpg", plotCorrelationNTIgM, width = 8, height = 6,
       dpi = 600)
ggsave("output/SFig_NT_IgM.jpg", plotCorrelationNTDeltaIgM, width = 4, height = 3,
       dpi = 600)


##
# Calculate correlation coefficients between IFA and Multiplex assay
correlationIFA <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata", "isotype", "panel")) %>% 
  select(assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(!(analyte %in% c("VACV"))) %>% 
  # filter(isotype == "IgG") %>% 
  mutate(titer = if_else(titer < 320, (paste0(intToUtf8(8804), "80")), as.character(titer)),
         titer = factor(titer, levels = c(paste0(intToUtf8(8804), "80"), 
                                          "320", "1280", "5120", "20480"),
                        ordered = TRUE)) %>% 
  mutate(titer_bin = case_when(titer == paste0(intToUtf8(8804), "80") ~ 1,
                               titer == "320" ~ 2,
                               titer == "1280" ~ 3,
                               titer == "5120" ~ 4,
                               titer == "20480" ~ 5)) %>% 
  group_by(analyte, isotype) %>% 
  summarise(correlation = list(cor.test(titer_bin, log10(dataIn), method = 'spearman')))

p <- vector()
estimate <- vector()
for(i in 1:nrow(correlationIFA)){
  p[i] <-correlationIFA$correlation[[i]]$p.value
  estimate[i] <-correlationIFA$correlation[[i]]$estimate
}


correlationIFAdataframe <-
  data.frame(Analyte = correlationIFA$analyte,
             Isotype = correlationIFA$isotype,
             rho =  round(estimate, digits = 2),
             p =  round(p, digits = 5)) %>% 
  arrange(Isotype, desc(rho))

export(correlationIFAdataframe, file = ("output/correlationIFA.xlsx"))


correlationNTExportIgG <- 
correlationNT %>% 
  select(Analyte = analyte,
         rho = rho_estimate,
         p = p_estimate,
         pearson_r = rho_estimate_pearson,
         pearson_p = p_estimate_pearson) %>% 
  mutate(rho =  round(rho, digits = 2),
         p =  round(p, digits = 5),
         pearson_r =  round(pearson_r, digits = 2),
         pearson_p =  round(pearson_p, digits = 5))  %>% 
  mutate(Isotype = "IgG") %>% 
  arrange(Isotype, desc(rho)) 

correlationNTExportIgM <- 
  correlationNTIgM %>% 
  select(Analyte = analyte,
         rho = rho_estimate,
         p = p_estimate,
         pearson_r = rho_estimate_pearson,
         pearson_p = p_estimate_pearson) %>% 
  mutate(rho =  round(rho, digits = 2),
         p =  round(p, digits = 5),
         pearson_r =  round(pearson_r, digits = 2),
         pearson_p =  round(pearson_p, digits = 5))  %>% 
  mutate(Isotype = "IgM") %>% 
  arrange(Isotype, desc(rho)) 

correlationNTExport <-
  rbind(correlationNTExportIgG, correlationNTExportIgM)
export(correlationNTExport, file = ("output/correlationNTExport.xlsx"))

## Perform ROC Analysis in comparison to IFA
# Select data and transform to wide format for IgG
rocDataframeIgG <- 
  dataInputQuant %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, panel) %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata" = "sampleID_metadata_mut", "isotype", "panel")) %>% 
  filter(dilution_assay == 100) %>% 
  filter(!is.na(titer)) %>% 
  dplyr::select(-dilution_assay, -sampleID_metadata.y) %>%
  filter(isotype == "IgG") %>%
  #  filter(assaytype == "Multiplex") %>% 
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  mutate(analyte = str_replace(analyte, "-", ".")) %>%
  mutate(analyte = if_else(assaytype == "ELISA", paste(analyte, "ELISA", sep = "_"), analyte)) %>%
  mutate(serostatus = if_else(titer < 320, 0, 1)) %>%
  group_by(sampleID_metadata, panel, titer, serostatus, analyte) %>% 
  summarise(meandataIn = mean(dataIn, na.rm = TRUE)) %>% 
  pivot_wider(names_from = analyte, values_from = meandataIn) %>% 
  filter(!is.na(A27L)) %>% 
  ungroup()

# Perform ROC for IgG
outplot_listIgG <- list()
j = 1
for(i in c("A27L", "A29",
           "D8L", "E8",
           "H3L", "L1R", "M1",
           "A33R", "A35R",
           "B5R", "B6",
           "A5L", "ATI.C", "ATI.N", "Delta", "Delta_ELISA")){
  outplot_listIgG[[j]] <- rocFunction(rocDataframeIgG, "serostatus",i)
  assign(paste("roc",  i, sep = "_"), outplot_listIgG[[j]])
  j = j+1
}

# Combined ROC for all multiplex-antigens by fitting a generalized linear model
rocCombinedGLMIgG <- 
  rocDataframeIgG %>%
  dplyr::select(-Delta_ELISA)
aucCombinedIgG <-glm(serostatus ~ 
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
                     data=rocCombinedGLMIgG, family="binomial")  #Logistic regression model
roc_CombinedIgG <- roc(serostatus ~ predict(aucCombinedIgG), data=rocCombinedGLMIgG)

summary(aucCombinedIgG)

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
paramDataFrameIFAIgG <- rocParamFunc(param_list)
paramDataFrame_transpose_sep <- transformRocParamFunc(param_list, "IFA IgG")
sensIFAIgG <- prepDataFrameFunc(paramDataFrame_transpose_sep, "sensitivity", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "sensitivity")
specIFAIgG <- prepDataFrameFunc(paramDataFrame_transpose_sep, "specificity", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "specificity")
accIFAIgG <- prepDataFrameFunc(paramDataFrame_transpose_sep, "accuracy", "ci_roc_") %>%
  mutate(isotype = "IgG") %>%
  mutate(parameter = "accuracy")

## Berechnung der ROC-Kurven für den IgM-Nachweis
## Auswahl der Daten für den Vergleich und Umformen in wide Format
rocDataframeIgM <- 
  dataInputQuant %>% 
  select(assaytype, dilution_assay, isotype, analyte,
         sampleID_metadata, dataIn, panel) %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata" = "sampleID_metadata_mut", "isotype", "panel")) %>% 
  filter(dilution_assay == 100) %>% 
  filter(!is.na(titer)) %>% 
  dplyr::select(-dilution_assay, -sampleID_metadata.y) %>%
  filter(isotype == "IgM") %>%
  filter(!(analyte %in% c("VACV", "VACV_ELISA"))) %>%
  mutate(analyte = str_replace(analyte, "-", ".")) %>%
  mutate(analyte = if_else(assaytype == "ELISA", paste(analyte, "ELISA", sep = "_"), analyte)) %>%
  mutate(serostatus = if_else(titer < 320, 0, 1)) %>%
  group_by(sampleID_metadata, panel, titer, serostatus, analyte) %>% 
  summarise(meandataIn = mean(dataIn, na.rm = TRUE)) %>% 
  pivot_wider(names_from = analyte, values_from = meandataIn) %>% 
  filter(!is.na(A27L)) %>% 
  ungroup()


# 1) Calculate ROCs compared to NT-Test ####
#names(dataInput)
outplot_listIgM <- list()
j = 1
for(i in c("A27L", "A29",
           "D8L", "E8",
           "H3L", "L1R", "M1",
           "A33R", "A35R",
           "B5R", "B6",
           "A5L", "ATI.C", "ATI.N", "Delta", "Delta_ELISA")){
  outplot_listIgM[[j]] <- rocFunction(rocDataframeIgM, "serostatus",i)
  assign(paste("roc_IgM",  i, sep = "_"), outplot_listIgM[[j]])
  j = j+1
}

# Combined ROC for all multiplex-antigens by fitting a generalized linear model
rocCombinedGLMIgM <- 
  rocDataframeIgM %>%
  dplyr::select(-Delta_ELISA) %>%
  filter(!is.na(A27L))
aucCombinedIgM <-glm(serostatus ~ 
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
                     data=rocCombinedGLMIgM, family="binomial")  #Logistic regression model
roc_IgM_Combined <- roc(serostatus ~ predict(aucCombinedIgM), data=rocCombinedGLMIgM)

# Calculate 95% confidence interval
k = 1
ci_list <- list()
for(k in c(1:length(ls(all.names = TRUE, pattern = ("^roc_IgM"))))){
  ci_list[[k]] <- ci.coords(get(ls(all.names = TRUE, pattern = ("^roc_IgM"))[k]), x="best", input = "threshold", ret=rets,  best.method = c("closest.topleft"),
                            best.policy = "random")
  assign(paste("ci",  ls(all.names = TRUE, pattern = ("^roc_IgM"))[k], sep = "_"), ci_list[[k]])
  k = k+1
}

# Export results of roc-analysis
k = 1
param_listIgM <- list()
for(k in c(1:length(ls(all.names = TRUE, pattern = ("^ci_roc_IgM"))))){
  param_listIgM[[k]] <- readRocParam(ls(all.names = TRUE, pattern = ("^ci_roc_IgM"))[k])
  k = k+1
}
# Save parameters in a dataframe
paramDataFrameIFAIgM <- rocParamFunc(param_listIgM)
paramDataFrame_transpose_sep_IgM <- transformRocParamFunc(param_listIgM, "IFA IgM")
sensIFAIgM <- prepDataFrameFunc(paramDataFrame_transpose_sep_IgM, "sensitivity", "ci_roc_IgM") %>%
  mutate(isotype = "IgM") %>%
  mutate(parameter = "sensitivity")
specIFAIgM <- prepDataFrameFunc(paramDataFrame_transpose_sep_IgM, "specificity", "ci_roc_IgM") %>%
  mutate(isotype = "IgM") %>%
  mutate(parameter = "specificity")
accIFAIgM <- prepDataFrameFunc(paramDataFrame_transpose_sep_IgM, "accuracy", "ci_roc_IgM") %>%
  mutate(isotype = "IgM") %>%
  mutate(parameter = "accuracy")
#########


threshold_IFA_IgG <-
  paramDataFrame_transpose_sep %>%
  mutate(assay = "IFA IgG") %>%
  dplyr::select(assay, antigene, CI, threshold) %>%
  # filter(grepl("SARS.CoV.2", antigene)) %>%
  mutate(antigene = str_remove(antigene, "ci_roc_")) 


threshold_IFA_IgM <-
  paramDataFrame_transpose_sep_IgM %>%
  mutate(assay = "IFA IgM") %>%
  dplyr::select(assay, antigene, CI, threshold) %>%
  # filter(grepl("SARS.CoV.2", antigene)) %>%
  mutate(antigene = str_remove(antigene, "ci_roc_IgM_")) 

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


# Plot ROC IgM
plotrocA27LIgM <- 
  plotROCfunction(list("A27L" = `roc_IgM_A27L`,
                       "A29" = `roc_IgM_A29`), "IgM")
plotrocA33RIgM <- 
  plotROCfunction(list("A33R" = `roc_IgM_A33R`,
                       "A35R" = `roc_IgM_A35R`), "IgM")
plotrocB5IgM <- 
  plotROCfunction(list("B5R" = `roc_IgM_B5R`,
                       "B6" = `roc_IgM_B6`), "IgM")
plotrocA5LIgM <- 
  plotROCfunction(list("A5L" = `roc_IgM_A5L`,
                       "ATI-C" = `roc_IgM_ATI.C`,
                       "ATI-N" = `roc_IgM_ATI.N`), "IgM")
plotrocL1RIgM <- 
  plotROCfunction(list("L1R" = `roc_IgM_L1R`,
                       "M1" = `roc_IgM_M1`), "IgM")
plotrocD8LIgM <- 
  plotROCfunction(list("D8L" = `roc_IgM_D8L`,
                       "E8" = `roc_IgM_E8`,
                       "H3L" = `roc_IgM_H3L`), "IgM")
plotrocLysateIgM <- 
  plotROCfunction(list("Delta" = `roc_IgM_Delta`), "IgM")



plotROCLysate <-
  ggarrange(plotrocLysate, plotrocLysateIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCA27LCombined <-
  ggarrange(plotrocA27L, plotrocA27LIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCL1RCombined <-
  ggarrange(plotrocL1R, plotrocL1RIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCD8LRCombined <-
  ggarrange(plotrocD8L, plotrocD8LIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCA33RCombined <-
  ggarrange(plotrocA33R, plotrocA33RIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCB5Combined <-
  ggarrange(plotrocB5, plotrocB5IgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROCA5LCombined <-
  ggarrange(plotrocA5L, plotrocA5LIgM, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

plotROC <-
  ggarrange(plotROCLysate, plotROCA27LCombined, plotROCL1RCombined,
            plotROCD8LRCombined, plotROCA33RCombined, plotROCB5Combined, plotROCA5LCombined,
            ncol = 7,
            align = "hv")

plotROCIgG <-
  ggarrange(plotrocLysate,
            plotrocA27L,
            plotrocL1R,
            plotrocD8L,
            plotrocA33R,
            plotrocB5,
            plotrocA5L, 
            ncol = 4, nrow = 2,
            align = "hv")

plotROCIgM <-
  ggarrange(plotrocLysateIgM,
            plotrocA27LIgM,
            plotrocL1RIgM,
            plotrocD8LIgM,
            plotrocA33RIgM,
            plotrocB5IgM,
            plotrocA5LIgM, 
            ncol = 4, nrow = 2,
            align = "hv")


ggsave("output/plotROC.png", plotROC, width = 22, height = 7, dpi = 300)
ggsave("output/plotROCIgG.png", plotROCIgG, width = 10, height = 7, dpi = 300)
ggsave("output/plotROCIgM.png", plotROCIgM, width = 10, height = 7, dpi = 300)

##
# Combine Assay performance parameters in a dataframe and plot values
sensSpec <-
  sensIFAIgG %>% 
  add_row(sensIFAIgM) %>% 
  add_row(specIFAIgG) %>% 
  add_row(specIFAIgM) %>% 
  add_row(accIFAIgG) %>% 
  add_row(accIFAIgM) %>% 
  mutate(antigene = trimws(antigene),
         antigene = case_when(antigene == "CombinedIgG" ~ "Combined",
                              TRUE ~ antigene)) %>% 
  mutate(antigene_class = case_when(antigene %in% c("A27L", "A29", "D8L", "E8",
                                                    "H3L", "L1R", "M1") ~ "IMV",
                                    antigene %in% c("A33R", "A35R", "B5R", "B6") ~ "EEV",
                                    antigene %in% c("A5L") ~ "Core",
                                    antigene %in% c("ATI-C", "ATI-N") ~ "Cytosol",
                                    antigene %in% c("Delta", "Delta ELISA",
                                                    "Combined") ~ "Combined")) %>% 
  mutate(antigene = factor(antigene, levels = c("A27L",
                                                "A29", 
                                                "L1R", 
                                                "M1",
                                                "D8L",
                                                "E8",
                                                "H3L",
                                                "A33R", 
                                                "A35R",
                                                "B5R", 
                                                "B6", 
                                                "A5L",
                                                "ATI-C", 
                                                "ATI-N",
                                                "Delta", 
                                                "Delta ELISA",
                                                "Combined"), ordered = TRUE))

#### Plot parameters specific ####
plotParameters <-
  sensSpec %>%
  group_by(isotype, parameter) %>%
  arrange(median) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM"))) %>%
  mutate(antigene_class = if_else(antigene_class == "SARS-CoV-2", "Multiplex", antigene_class)) %>%
  mutate(antigene_class = factor(antigene_class, levels = c("IMV", 
                                                            "Combined",
                                                            "Cytosol",
                                                            "Core",
                                                            "EEV"), ordered = TRUE)) %>%
  mutate(median = median * 100,
         low = low * 100,
         high = high * 100) %>%
  ggplot(mapping = aes(x = (antigene), y = median, ymin = low, ymax = high, color = antigene_class)) +
  geom_pointrange() +
  scale_y_continuous(name = "Value (%)")+
  scale_x_discrete(name = "")+
  labs(color = "")+
  # coord_flip() +
  scale_color_colorblind()+
  theme_bw() +
  theme(legend.position="top") +
  facet_grid(parameter ~ isotype)+
  theme(strip.background = element_blank())


plotMethCompParam <- 
  ggarrange(plotMethComp, plotParameters, nrow = 2,
            labels = c("", "c", heigths = c(3, 1)))

ggsave("output/plotMethCompParam .png", plotMethCompParam , width = 14, height = 12, dpi = 300)

## 
# Export dataframe with assay parameters and cutoff values
sensSpecSens <- 
  sensSpec %>% 
  filter(!(grepl("ELISA", antigene))) %>% 
  filter(antigene != "Combined") %>% 
  filter(parameter == "sensitivity") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype) %>% 
  mutate(low = round(low, digits = 2),
         median = round(median, digits = 2),
         high = round(high, digits = 2)) %>% 
  arrange(Isotype, Analyte)

sensSpecSpec <- 
  sensSpec %>% 
  filter(!(grepl("ELISA", antigene))) %>% 
  filter(antigene != "Combined") %>% 
  filter(parameter == "specificity") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype) %>% 
  mutate(low = round(low, digits = 2),
         median = round(median, digits = 2),
         high = round(high, digits = 2)) %>% 
  arrange(Isotype, Analyte)


sensSpecAccu <- 
  sensSpec %>% 
  filter(!(grepl("ELISA", antigene))) %>% 
  filter(antigene != "Combined") %>% 
  filter(parameter == "accuracy") %>% 
  select(Analyte = antigene, low, median, high, Isotype = isotype) %>% 
  mutate(low = round(low, digits = 2),
         median = round(median, digits = 2),
         high = round(high, digits = 2)) %>% 
  arrange(Isotype, Analyte)

export(sensSpec, "output/sensSpec.xlsx")
export(sensSpecSens, "output/sensSpecSens.xlsx")
export(sensSpecSpec, "output/sensSpecSpec.xlsx")
export(sensSpecAccu, "output/sensSpecAccu.xlsx")
export(threshold_IFA_IgG, "output/threshold_IFA_IgG.xlsx")
export(threshold_IFA_IgM, "output/threshold_IFA_IgM.xlsx")

##
# Determne Serostatus based on cutoff values
# 

threshold <-
  rbind(threshold_IFA_IgG, threshold_IFA_IgM) %>% 
  mutate(assay = str_remove(assay, "IFA ")) %>% 
  mutate(antigene = str_replace(antigene, "\\.", "-")) %>% 
  mutate(assaytype = if_else(grepl("ELISA", antigene), "ELISA", "Multiplex"),
         antigene = str_remove(antigene, "_ELISA")) %>% 
  pivot_wider(names_from = CI, values_from = threshold)

dataInputQuantCat <-
  dataInputQuant %>% 
  left_join(threshold, by = c("isotype" = "assay", "analyte" = "antigene", "assaytype")) %>% 
  mutate(serostatus = if_else(dataIn > median, 1, 0),
         serostatus_cat = case_when(dataIn > high ~ "positive",
                                    dataIn > median & dataIn <= high ~ "borderline positive",
                                    dataIn > low & dataIn <= median ~ "borderline negative",
                                    dataIn <= low ~ "negative"),
         serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"),
                                 ordered = TRUE))

save(dataInputQuantCat, file = "output/dataInputQuantCat.Rdata")

dataInputQuantCatNK <-
  dataInputQuantNK %>% 
  mutate(dataIn = case_when(assaytype == "Multiplex" ~ quant_conc_fix,
                            assaytype == "ELISA" ~ quant_conc_fix)) %>% 
  left_join(threshold, by = c("isotype" = "assay", "analyte" = "antigene", "assaytype")) %>% 
  mutate(serostatus = if_else(dataIn > median, 1, 0),
         serostatus_cat = case_when(dataIn > high ~ "positive",
                                    dataIn > median & dataIn <= high ~ "borderline positive",
                                    dataIn > low & dataIn <= median ~ "borderline negative",
                                    dataIn <= low ~ "negative"),
         serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                            "borderline positive", "positive"),
                                 ordered = TRUE))

save(dataInputQuantCat, file = "output/dataInputQuantCat.Rdata")
save(dataInputQuantCatNK, file = "output/dataInputQuantCatNK.Rdata")

save(threshold, file = "output/threshold.Rdata")

export(threshold, "output/threshold.xlsx")

##
# Analyse and export ELISA results only
dataInputELISA <-
  dataInputQuantCat %>% 
  filter(assaytype == "ELISA") %>% 
  filter(dilution_assay == 100) %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata" = "sampleID_metadata_mut", "isotype", "panel")) %>% 
  filter(!is.na(titer))

export(dataInputELISA, "output/dataInputELISA.xlsx")

samplesELISA <- unique(dataInputELISA$sampleID_metadata)
samplesIFA <- unique(dataInputELISA$sampleID_metadata[!is.na(dataInputELISA$titer)])
samplesRica <- c("P-22-01266-001-16", 
                 "P-22-1037-001-02",
                 "P-22-049-001-14",
                 "P-22-01432-001-10",
                 "P-22-01095-001-15",
                 "P-20-002-001-02",
                 "P-20-162-001-01",
                 "P-18-058-001-01",
                 "P-18-052-002-02", 
                 "P-17-077-001-01",
                 "P-19-079-001-02",
                 "S-20-004-001-03",
                 "S-20-004-010-03",
                 "S-20-004-011-02",
                 "S-20-004-012-03",
                 "S-20-004-015-03")
samplesRica %in% samplesELISA
samplesRica %in% samplesIFA

mfiCutoffHigh <-
dataInputQuantCat %>% 
  filter(analyte == "Delta" & dilution_assay == 100 & assaytype == "Multiplex") %>% 
  filter(log10(dataIn) >= 2.08 & log10(dataIn) < 2.1) %>% 
  summarise(mean_data = mean(data),
            sd_data = sd(data),
            n_data = length(data)) %>% 
  mutate(cutoff = "2.08")


mfiCutoffPositive <-
  dataInputQuantCat %>% 
  filter(analyte == "Delta" & dilution_assay == 100 & assaytype == "Multiplex") %>% 
  filter(log10(dataIn) >= 1.48199104 & log10(dataIn) < 1.5) %>% 
  summarise(mean_data = mean(data),
            sd_data = sd(data),
            n_data = length(data))%>% 
  mutate(cutoff = "1.48199104")

mfiCutoff <-
  rbind(mfiCutoffHigh, mfiCutoffPositive)

export(mfiCutoff, "output/mfiCutoff.xlsx")
