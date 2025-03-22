##
# Method comparison for establishent of Mpox-Multiplex-Assays
# Daniel Stern
# 2025/03/01
# Version 1
##

## 
# Clean environment
rm(list = ls(all.names = TRUE))

## 
# Load packages
library(rio)
library(tidyverse)
library(mcr)
library(ggpubr)
library(ggthemes)
library(corrplot)
library(rstatix)


# Load data and metadata
load("input/metadata_IFA.Rdata")
load("input/dataInputNT.Rdata")
load("input/dataInputELISAMultiplex.Rdata")
load("input/dataInputQuant.Rdata")

# Define functions
# Helper functions to transform dataframe into matrix
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  dimnames(m)[1]<-x[,1]
  m
}


####
# 1) Plot correlations between reference methods and multiplex antigens
# Fig 1b: Correlation between ELISA and mutliplex data
dataInputCorrMultiplexIgG <- 
  dataInputELISAMultiplex %>%
  filter(isotype == "IgG") %>% 
  filter(analyte != "VACV") %>% 
  select(sampleID_metadata, analyte, dataIn, dataIn_ELISA) %>% 
  mutate(dataIn = log10(dataIn),
         dataIn_ELISA = log10(dataIn_ELISA)) %>% 
  rename(ELISA = dataIn_ELISA) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(-sampleID_metadata)

# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixMultiIgG <- cor_pmat(dataInputCorrMultiplexIgG, method = "pearson")
matrixMultiIgG <- cor_mat(dataInputCorrMultiplexIgG, method = "pearson")

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationMultiIgG.pdf"
pdf(height=4, width=4, file=file_path)
corrplot(matrix.please(matrixMultiIgG), p.mat = matrix.please(pmatrixMultiIgG), method = "circle", order = "hclust",diag = F,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05,
         type = "lower")
dev.off()


# Fig 1c: Plot correlation coefficients between IFA titers (log2) and Multiplex results 
dataInputCorIFAIgG <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata", "isotype", "panel")) %>% 
  dplyr::select(sampleID_metadata, assaytype, dilution_assay, isotype, analyte, dataIn, titer, remarks) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(dilution_assay == 100, !is.na(titer)) %>% 
  filter(!(analyte %in% c("VACV"))) %>% 
  filter(isotype == "IgG") %>%
  mutate(dataIn = log10(dataIn),
         titer = if_else(titer == 0, 0, log2(titer))) %>% 
  dplyr::select(sampleID_metadata, analyte, dataIn, titer) %>%
  unique() %>%
  pivot_wider(names_from = "analyte", values_from = "dataIn", values_fn = mean) %>%
  dplyr::select(-sampleID_metadata) %>% 
  rename(IFA = titer)

# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixnegIFAIgG <- cor_pmat(dataInputCorIFAIgG, method = "spearman")
matrixnegIFAIgG <- cor_mat(dataInputCorIFAIgG, method = "spearman")

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationIFAIgG.pdf"
pdf(height=4, width=4, file=file_path)
corrplot(matrix.please(matrixnegIFAIgG), p.mat = matrix.please(pmatrixnegIFAIgG), method = "circle", order = "hclust",diag = F,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05,
         type = "lower")
dev.off()


# Fig 1d: Correlation between NT and multiplex data
dataInputCorNTIgG <-
  dataInputNT %>%
  dplyr::filter(dilution_assay == 100 
                & assaytype == "Multiplex"
                & isotype == "IgG" &
                  analyte != "VACV") %>%
  dplyr::select(sampleID_metadata, analyte, dataIn, NT_titer) %>%
  unique() %>%
  mutate(dataIn = log10(dataIn),
         NT_titer = log2(NT_titer)) %>% 
  pivot_wider(names_from = "analyte", values_from = "dataIn") %>%
  dplyr::select(-sampleID_metadata) %>% 
  rename(NT = NT_titer)

# Calculate correlation matrix and matrix of p-values (package rstatix)
pmatrixNTIgG <- cor_pmat(dataInputCorNTIgG)
matrixNTIgG <- cor_mat(dataInputCorNTIgG)

# Save corrplot to file -> combination in Adobe Illustrator if suggestion is accepted
file_path= "output/plotCorrelationNTIgG.pdf"
pdf(height=4, width=4, file=file_path)
corrplot(matrix.please(matrixNTIgG), p.mat = matrix.please(pmatrixNTIgG), method = "circle", order = "hclust",diag = F,
         hclust.method = c("ward.D2"), tl.col = "black", insig = "pch", pch.col = "grey", sig.level = 0.05,
         type = "lower")
dev.off()


####
# 2) Plot regression between Delta results from multiplex assay and 
# reference assays
# Select IgG and Delta 
dataDeltaIgG <-
  dataInputELISAMultiplex %>% 
  filter(isotype == "IgG" & analyte == "Delta")

# Perform Passing Bablok Regression for IgG data
pbBaRegDeltaIgG <- mcreg((log10(dataDeltaIgG$dataIn_ELISA)), 
                         log10(dataDeltaIgG$dataIn), method.reg = "PaBa")

# Select IgM and Delta 
dataDeltaIgM <-
  dataInputELISAMultiplex %>% 
  filter(isotype == "IgM" & analyte == "Delta") 

# Perform Passing Bablok Regression for IgM data
pbBaRegDeltaIgM <- mcreg(log10(dataDeltaIgM$dataIn_ELISA), 
                         log10(dataDeltaIgM$dataIn), method.reg = "PaBa")

# Write parameter for export in dataframe
parampaBa <- data.frame(isotype = c("IgG", "IgM"),
                        pearsonsR = c(0.933, 0.846),
                        slope = c(0.87, 0.89),
                        intercept = c(0.13, -0.63),
                        n = c(75, 75))

export(parampaBa, "output/parampaBa.xlsx")

## 
# Generate ggplots for Passing Bablok Regression 
plotpaBaIgG <-
  dataDeltaIgG %>% 
  ggplot(aes(x = log10(dataIn_ELISA), y = log10(dataIn)))+ 
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = colorblind_pal()(8)[4])+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,1], intercept = pbBaRegDeltaIgG@para[1,1], color = "grey25", size = 1, alpha = 0.5)+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,3], intercept = pbBaRegDeltaIgG@para[1,3], lty = "dashed", color = "grey25")+
  geom_abline(slope = pbBaRegDeltaIgG@para[2,4], intercept = pbBaRegDeltaIgG@para[1,4], lty = "dashed", color = "grey25")+ 
  scale_x_continuous(name = "ELISA quant.", limits = c(0,5)) +
  scale_y_continuous(name = "Multiplex Delta quant.", limits = c(0,5)) +
  labs(subtitle = "IgG") +
  theme_pubr() +
  theme(plot.subtitle=element_text(hjust=0.5)) 

plotpaBaIgM <-
  dataDeltaIgM %>% 
  ggplot(aes(x = log10(dataIn_ELISA), y = log10(dataIn)))+ 
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = colorblind_pal()(8)[4])+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,1], intercept = pbBaRegDeltaIgM@para[1,1], color = "grey25", size = 1, alpha = 0.5)+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,3], intercept = pbBaRegDeltaIgM@para[1,3], lty = "dashed", color = "grey25")+
  geom_abline(slope = pbBaRegDeltaIgM@para[2,4], intercept = pbBaRegDeltaIgM@para[1,4], lty = "dashed", color = "grey25")+ 
  scale_x_continuous(name = "ELISA quant.", limits = c(0,5)) +
  scale_y_continuous(name = "Multiplex Delta quant.", limits = c(0,5)) +
  labs(subtitle = "IgM") +
  theme_pubr() +
  theme(plot.subtitle=element_text(hjust=0.5))

plotpaBacombined <-
  ggarrange(plotpaBaIgG, plotpaBaIgM, ncol = 2, align = "hv")


##
# Plots for methodcomparison with IFA
plotIFAIgG <-
  dataInputQuant %>% 
  left_join(metadata_IFA, by = c("sampleID_metadata" = "sampleID_metadata", "isotype", "panel")) %>% 
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
  scale_fill_manual(name = "Assay", values = colorblind_pal()(8)[2:3]) + 
  theme(strip.background = element_blank()) +
  xlab("IFA (Titer)")+
  scale_y_continuous(name = "Multiplex Delta quant.", limits = c(0,5)) +
  labs(subtitle = "IgG") +
  theme_pubr() +
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
  scale_fill_manual(name = "Assay", values = colorblind_pal()(8)[3], drop = F) + 
  theme(strip.background = element_blank()) +
  scale_y_continuous(name = "Multiplex Delta quant.", limits = c(0,5)) +
  labs(subtitle = "IgM") +
  theme_pubr() +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(name = "IFA (Titer)")

# Plots for method comparison with NT
plotCorrelationNTDeltaIgG <-
  dataInputNT %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgG") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte == "Delta") %>% 
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
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
  theme_pubr() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex Delta quant.") +
  scale_fill_manual(name = "Year of birth", values = colorblind_pal()(8)[4]) +
  labs(subtitle = "IgG") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")

plotCorrelationNTDeltaIgM <-
  dataInputNT %>% 
  dplyr::filter(dilution_assay == 100) %>% 
  dplyr::filter(isotype == "IgM") %>% 
  dplyr::filter(assaytype == "Multiplex") %>% 
  dplyr::filter(analyte == "Delta") %>% 
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM", "IgA"))) %>%
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
  theme_pubr() +
  scale_x_discrete(name = "NT (Titer)") +
  scale_y_continuous(name = "Multiplex Delta quant.") +
  scale_fill_manual(name = "Year of birth", values = colorblind_pal()(8)[4]) +
  labs(subtitle = "IgM") +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")


##
# Plot comparison between delta antigen and reference methods
Fig1_plotDelta <-
  ggarrange(plotpaBaIgG, plotpaBaIgM ,
            plotIFAIgG, plotIFAIgM,
            plotCorrelationNTDeltaIgG, 
            plotCorrelationNTDeltaIgM, ncol = 6, align = "hv")

ggsave("output/Fig1_plotDelta.pdf", Fig1_plotDelta, width = 14, height = 3)


####
# Generate supporting figures
# Figure S1: Plot IgG and IgM multiplex antigens in comparison to ELISA
plotSFigELISAMultiplex <-
  dataInputELISAMultiplex %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  ggplot(mapping = aes(x = log10(dataIn_ELISA), y = log10(dataIn), color = isotype)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colorblind_pal()(8)[2:3]) +
  scale_x_continuous(name = "ELISA quant.") +
  scale_y_continuous(name = "Multiplex quant.", breaks = c(0,1,2,3,4,5,6)) +
  facet_grid(isotype  ~ analyte) +
  theme_pubr() +
  theme(#plot.subtitle=element_text(hjust=0.5),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none")

ggsave("output/FigS1_plotELISAMultiplex.png", plotSFigELISAMultiplex, 
       width = 14, height = 6)


