####
# Perform classical ROC analysis 
# Daniel Stern RKI
# 2025/02/25
####

####
library(rio)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(pROC)


rm(list = ls(all.names = TRUE))

# Source functions
source("functions/determineCutoffFunction.R", encoding = "UTF-8")
source("functions/rocFunction.R", encoding = "UTF-8")

## Load datainput generated for clustering in script "23 Spox Mpox Pos"
load("input/dataInWide.Rdata")

# Load plot for improvement of ML Predictions
# load("input/plotImprovePerfomance.Rdata")


####
# Perform ROC Analysis based on the wide dataframe
# Generate two dataframes: first dataframe for the prediction on the pre-panel
# Second dataframe: for the prediction on the MVA-panel based on ATI-N alone
# Comparison between "Pre" from acute panel only, younger population
# and "MPXV/MVA" for both panels
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
  filter(panel_ATI %in% c("Pre", "MPXV", "MPXV SPox", "MVA", "MVA SPox", "Pre SPox")) %>% 
  mutate(panel_Seropos = if_else(panel_ATI %in% c("MPXV SPox", "MPXV", "MVA", "MVA SPox"), "MPXV/MVA", "Pre"))


dataInWideAntigens %>% 
  select(sampleID_metadata, panel_ATI, panel_Seropos) %>% 
  group_by(panel_Seropos, panel_ATI) %>% 
  count()

# Comparison for all subjecte of all age_groups
# between MVA from the acute panel and MPXV from the acute and epi panel
dataInWideATI_N <- 
  dataInWide %>% 
  filter(panel_detail != "SPox") %>% 
  #  filter((Age_group %in% c("< 30", "< 40"))) %>% 
  mutate(panel_ATI = case_when(panel == "MPXV" ~ "MPXV",
                               panel == "MVA" & panel_detail == "MVA" ~ "MVA",
                               panel == "SPox" & panel_detail == "MPXV" ~ "MPXV SPox",
                               panel == "SPox" & panel_detail == "MVA" ~ "MVA SPox",
                               panel == "SPox" & panel_detail == "Pre" ~ "Pre SPox",
                               panel == "Pre_New" | (panel == "MVA" & panel_detail == "Pre") ~ "Pre"),
         panel_ATI = factor(panel_ATI, levels = c("Pre", "MVA", "Pre SPox", "MVA SPox", "MPXV", "MPXV SPox"),
                            ordered = TRUE)) %>% 
  filter(panel_ATI %in% c("MVA", "MPXV SPox", "MPXV", "MVA SPox")) %>% 
  mutate(panel_ATI_N = if_else(panel_ATI %in% c("MPXV SPox", "MPXV"), "MPXV", "MVA"))

dataInWideATI_N %>% 
  select(sampleID_metadata, panel_ATI, panel_ATI_N) %>% 
  group_by(panel_ATI_N) %>% 
  count()


##
# Perform ROC-Analysis for all antigens:
# Determine serostatus in younger subjects
rm(list = ls(pattern = c("roc_")))

roc_Delta <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_Delta)

roc_A27L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_A27L)

roc_A29L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_A29L)

roc_ATI_N_Serostatus <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$`dataIn_ATI-N`)

roc_ATI_C_Serostatus <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$`dataIn_ATI-C`)

roc_D8L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_D8L)

roc_E8L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_E8L)

roc_L1R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_L1R)

roc_M1R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_M1R)

roc_A5L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_A5L)

roc_H3L <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_H3L)

roc_A33R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_A33R)

roc_A35R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_A35R)

roc_B5R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_B5R)

roc_B6R <-
  roc(factor(dataInWideAntigens$panel_Seropos, levels = c("Pre", "MPXV/MVA"),
             ordered = TRUE), dataInWideAntigens$dataIn_B6R)


##
# Perform ROC-Analysis for ATI-N for the following antigen
# ATI-N
roc_ATI_N <-
  roc(factor(as.character(dataInWideATI_N$panel_ATI_N), levels = c("MVA", "MPXV"),
             ordered = TRUE),dataInWideATI_N$`dataIn_ATI-N`)

roc_E8L_Infected <-
  roc(factor(as.character(dataInWideATI_N$panel_ATI_N), levels = c("MVA", "MPXV"),
             ordered = TRUE),dataInWideATI_N$`dataIn_E8L`)

roc_A35R_Infected <-
  roc(factor(as.character(dataInWideATI_N$panel_ATI_N), levels = c("MVA", "MPXV"),
             ordered = TRUE),dataInWideATI_N$`dataIn_A35R`)

roc_B6R_Infected <-
  roc(factor(as.character(dataInWideATI_N$panel_ATI_N), levels = c("MVA", "MPXV"),
             ordered = TRUE),dataInWideATI_N$`dataIn_B6R`)


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

# Export AUC of roc-analysis
k = 1
auc_list <- list()
for(k in c(1:length(ls(all.names = TRUE, pattern = ("^roc_"))))){
  auc_list[[k]] <- as.numeric(get(ls(all.names = TRUE, pattern = ("^roc_"))[k])$auc)
  k = k+1
}


AUC = as_tibble(do.call(rbind, auc_list))
Antigen = as_tibble(ls(all.names = TRUE, pattern = ("^roc_")))

add_column(AUC, Antigen) %>% 
  rename(AUC = V1, Antigen = value) %>% 
  mutate(Antigen = str_remove(Antigen, "roc_"),
         Antigen = str_replace(Antigen, "ATI_", "ATI-")) %>% 
  separate(Antigen, into = c("Antigen", "Class"), sep = "_") %>% 
  mutate(Class = if_else(is.na(Class) & Antigen == "ATI-N", "Infected", Class),
         Class = if_else(is.na(Class), "Serostatus", Class))

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

export(file = "output/2025-02-25_Populations_based_Cutoff.xlsx", paramDataFrame_transpose_sep)

# Determine threshold 
threshold <-
  paramDataFrame_transpose_sep %>%
  dplyr::select(antigene, CI, threshold) %>%
  mutate(antigene = str_remove(antigene, "ci_roc_")) %>% 
  mutate(antigene = if_else(antigene == "ATI_N", "ATI-N", antigene))

threshold_median <-
  threshold %>% 
  filter(CI == "median")

ROC_parameters <-
  paramDataFrame_transpose_sep %>% 
  filter(CI == "median") %>% 
  select(Antigen = antigene, threshold, sensitivity, specificity, accuracy, threshold, 
         tn:fp) %>% 
  rowwise() %>% 
  mutate(n = sum(tn, tp, fn, fp),
         Antigen = str_remove(Antigen, "ci_roc_"),
         Antigen = str_replace(Antigen, "ATI_", "ATI-")) %>% 
  separate(Antigen, into = c("Antigen", "Class"), sep = "_") %>% 
  mutate(Class = if_else(is.na(Class) & Antigen == "ATI-N", "Infected", Class),
         Class = if_else(is.na(Class), "Serostatus", Class),
         sensitivity = round(sensitivity, digits = 2),
         specificity = round(specificity, digits = 2),
         accuracy = round(accuracy, digits = 2),
         tn = as.integer(tn),
         tp = as.integer(tp),
         fn = as.integer(fn),
         fp = as.integer(fp),
         n = as.integer(n)) %>% 
  arrange(Class)

export(ROC_parameters, file = "output/ROC_parameters_classic.xlsx")


threshold_population <- 
  threshold 

# Export results of threshold calculation for import in determination of specificities
save(threshold_population, file = "output/treshold_population_Rev.RData")


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
    theme_pubr() +
    theme(plot.subtitle=element_text(hjust=0.5),
          legend.position = "top") 
} 

plotrocAntigensGood <- 
  plotROCfunction(list("Delta" = roc_Delta,
                       "D8" = roc_D8L,
                       "E8" = roc_E8L,
                       "A33" = roc_A33R,
                       "A35" = roc_A35R,
                       "B5" = roc_B5R,
                       "B6" = roc_B6R), "IgG")


plotrocAntigensBad <- 
  plotROCfunction(list("A27" = roc_A27L,
                       "A29" = roc_A29L,
                       "L1" = roc_L1R,
                       "M1" = roc_M1R,
                       "H3" = roc_H3L,
                       "A5" = roc_A5L,
                       "ATI-N" = roc_ATI_N_Serostatus,
                       "ATI-C" = roc_ATI_C_Serostatus), "IgG")

plotrocMPXVMVA <- 
  plotROCfunction(list("ATI-N" = roc_ATI_N,
                       "A35" = roc_A35R_Infected,
                       "B6" = roc_B6R_Infected,
                       "E8" = roc_E8L_Infected), "IgG")

plotRocCombined <-
  ggarrange(plotrocAntigensGood, plotrocAntigensBad, plotrocMPXVMVA,
            ncol = 3, align = "hv", labels = "auto")

plotFig6 <-
  ggarrange(plotRocCombined, 
            align = "hv")

ggsave(filename = "output/plotFig6.png", plotFig6 ,
       width = 10, height = 4, dpi = 600)




