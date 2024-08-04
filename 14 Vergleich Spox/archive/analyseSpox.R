##
# Analysis SPox study
# Assay validation
# Compare results from predictions with metadata 
# Daniel Stern RKI
# 26.07.2023
##

# Prepare environment and load packages
rm(list = ls(all.names = TRUE))
library("rio")
library("tidyverse")
library("ggthemes")
library("ggpubr")
library("caret")
library("pROC")

# Source functions
source("functions/determineCutoffFunction.R")
source("functions/rocFunction.R")

# Import data
dataInSpox1 <- import("Zusammenfassung Ergebnisse SPox 1.xlsx")
dataInSpox2 <- import("Zusammenfassung Ergebnisse SPox 2.xlsx")

dataInMeasured <-
  rbind(dataInSpox1, dataInSpox2)

dataInMeta <- import("MPOX_clinic2.xlsx") %>% 
  rename(sampleID_meta = `LAB ID`)

dataIn <-
  dataInMeasured %>% 
  left_join(dataInMeta, by= c("sampleID_meta")) %>% 
  mutate(across(starts_with("dataIn_"), ~ log10(.x)))


unique(dataIn$Final_combined)
unique(dataIn$MPX_diagnose_final)

dataInConfusionMPXV <- 
  dataIn %>% 
  select(Final_combined, MPX_diagnose_final) %>% 
  mutate(MPXV_confusion = case_when(MPX_diagnose_final == "Yes" ~ "MPXV",
                                    MPX_diagnose_final == "No" ~ "Other"),
         Final_confusion = case_when(Final_combined == "MPXV" ~ "MPXV",
                                     Final_combined == "seronegative" ~ "seronegative",
                                     TRUE ~ "Other")) %>% 
  filter(MPXV_confusion %in% c("MPXV", "Other")) %>% 
  filter(Final_confusion %in% c("MPXV", "Other"))

confusionMatrix(data = as.factor(dataInConfusionMPXV$Final_confusion),
                reference = as.factor(dataInConfusionMPXV$MPXV_confusion))



dataInConfusionMVA <- 
  dataIn %>% 
  select(Final_combined, MPX_VACC_STATUS_1) %>% 
  mutate(MVA_confusion = case_when(MPX_VACC_STATUS_1 == "Yes" ~ "MVA",
                                   MPX_VACC_STATUS_1 == "No" ~ "Other"),
         Final_confusion = case_when(Final_combined == "MVA" ~ "MVA",
                                     Final_combined == "seronegative" ~ "seronegative",
                                     TRUE ~ "Other")) %>% 
  filter(MVA_confusion %in% c("MVA", "Other")) %>% 
  filter(Final_confusion %in% c("MVA", "Other"))


confusionMatrix(data = as.factor(dataInConfusionMVA$Final_confusion),
                reference = as.factor(dataInConfusionMVA$MVA_confusion))


# Generate plot for different groups

dataIn %>% 
  pivot_longer(starts_with("MPXV_"), names_to = "LDA", values_to = "Probability", 
               names_prefix = "MPXV_") %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  ggplot(mapping = aes(x = LDA, y = Probability, fill = MPX_diagnose_final)) +
  geom_boxplot(notch = F, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75
  ), alpha = 0.3) +
  theme_bw()+
  # facet_grid(. ~ Antigene)+
  #  scale_x_discrete(labels = c("Admission", "d7", "d10", "d14", "d21"))+
  scale_fill_manual(name = "MPXV Diagnose", values = colorblind_pal()(9)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("MPXV Probability") #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

dataIn %>% 
  pivot_longer(starts_with("MVA_"), names_to = "LDA", values_to = "Probability", 
               names_prefix = "MVA_") %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  ggplot(mapping = aes(x = LDA, y = Probability, fill = MPX_diagnose_final)) +
  geom_boxplot(notch = F, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75
  ), alpha = 0.3) +
  theme_bw()+
  # facet_grid(. ~ Antigene)+
  #  scale_x_discrete(labels = c("Admission", "d7", "d10", "d14", "d21"))+
  scale_fill_manual(name = "MPXV Diagnose", values = colorblind_pal()(9)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("MVA Probability")


dataIn %>% 
  pivot_longer(starts_with("Pre_"), names_to = "LDA", values_to = "Probability", 
               names_prefix = "Pre_") %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  ggplot(mapping = aes(x = LDA, y = Probability, fill = MPX_diagnose_final)) +
  geom_boxplot(notch = F, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75
  ), alpha = 0.3) +
  theme_bw()+
  # facet_grid(. ~ Antigene)+
  #  scale_x_discrete(labels = c("Admission", "d7", "d10", "d14", "d21"))+
  scale_fill_manual(name = "MPXV Diagnose", values = colorblind_pal()(9)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("Pre Probability")


dataIn %>% 
  pivot_longer(starts_with("ratio_"), names_to = "LDA", values_to = "Probability", 
               names_prefix = "ratio_") %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  mutate(Probability = log10(Probability)) %>% 
  ggboxplot(x = "MPX_diagnose_final", y = "Probability", color = "LDA",
            add = "point") +
  # geom_boxplot(notch = F, outlier.shape = NA) +
  #  geom_point(position = position_dodge(width = 0.75
  # ), alpha = 0.3) +
  theme_bw()+
  stat_compare_means(comparisons = list(c("No", "Yes")
  )) +
  facet_grid(.~LDA) +
  scale_color_manual(name = "MPXV Diagnose", values = colorblind_pal()(8)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("Ratio") 


dataIn %>% 
  pivot_longer(starts_with("MPXV_"), names_to = "LDA", values_to = "Probability", 
               names_prefix = "MPXV_") %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  mutate(Probability = (Probability)) %>% 
  ggboxplot(x = "MPX_diagnose_final", y = "Probability", color = "LDA",
            add = "point") +
  # geom_boxplot(notch = F, outlier.shape = NA) +
  #  geom_point(position = position_dodge(width = 0.75
  # ), alpha = 0.3) +
  theme_bw()+
  stat_compare_means(comparisons = list(c("No", "Yes")
  )) +
  facet_grid(.~LDA) +
  scale_color_manual(name = "MPXV Diagnose", values = colorblind_pal()(8)[2:7]) + 
  theme(strip.background = element_blank()) +
  xlab("")+
  ylab("Probability") 


##  
# Perform ROC Analysis with different data inputs

rocDataframe <-
  dataIn %>% 
  filter(serostatus.delta == 1) %>% 
  dplyr::select(dataIn_D8L:dataIn_Delta,
                MPXV_all, MPXV_delta, MPXV_high,
                ratio_A29_A27L, 
                ratio_E8_D8L,
                ratio_A35R_A33R, 
                ratio_B6_B5R,
                GLM_value, MPX_diagnose_final) %>% 
  filter(!is.na(MPX_diagnose_final)) %>% 
  mutate(MPX_diagnose_final = if_else(MPX_diagnose_final == "Yes", 1, 0))

names(rocDataframe) <- str_replace(names(rocDataframe), "-", "_")

names(rocDataframe)[1:23]

# Perform ROC for All sera
outplot_listIgG <- list()
j = 1
for(i in names(rocDataframe)[1:23]){
  outplot_listIgG[[j]] <- rocFunction(rocDataframe, "MPX_diagnose_final",i)
  assign(paste("roc",  i, sep = "_"), outplot_listIgG[[j]])
  j = j+1
}

aucCombined <-glm(MPX_diagnose_final ~ 
                    dataIn_D8L+ dataIn_L1R   +  dataIn_B5R    +  dataIn_A33R  +  
                    dataIn_M1  +     `dataIn_ATI_N`   + dataIn_B6  +     dataIn_E8 +      
                    dataIn_A5L +      dataIn_A29    +  dataIn_H3L   +   dataIn_A35R +   
                    `dataIn_ATI_C` +   dataIn_A27L +    dataIn_Delta,
                  data=rocDataframe, family="binomial")  #Logistic regression model
roc_Combined <- roc(MPX_diagnose_final ~ predict(aucCombined), data=rocDataframe)


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
paramDataFrame <- rocParamFunc(param_list)
paramDataFrame_transpose_sep <- transformRocParamFunc(param_list, "All")
sens <- prepDataFrameFunc(paramDataFrame_transpose_sep, "sensitivity", "ci_roc_") %>%
  mutate(panel = "All") %>%
  mutate(parameter = "sensitivity")
spec <- prepDataFrameFunc(paramDataFrame_transpose_sep, "specificity", "ci_roc_") %>%
  mutate(panel = "All") %>%
  mutate(parameter = "specificity")
acc <- prepDataFrameFunc(paramDataFrame_transpose_sep, "accuracy", "ci_roc_") %>%
  mutate(panel = "All") %>%
  mutate(parameter = "accuracy")

threshold_IgGAll <-
  paramDataFrame_transpose_sep %>%
  dplyr::select(assay, antigene, CI, threshold) %>%
  # filter(grepl("SARS.CoV.2", antigene)) %>%
  mutate(antigene = str_remove(antigene, "ci_rocAll_")) 


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

plotrocAllIMV <- 
  plotROCfunction(list("A27L" = roc_dataIn_A27L,
                       "A29" = roc_dataIn_A29, 
                       "L1R" = roc_dataIn_L1R,
                       "M1" = roc_dataIn_M1,
                       "D8L" = roc_dataIn_D8L,
                       "E8L" = roc_dataIn_E8,
                       "H3L" = roc_dataIn_H3L),  "All")

plotrocAllEEV <- 
  plotROCfunction(list("A33R" = roc_dataIn_A33R,
                       "A35R" = roc_dataIn_A35R, 
                       "B5R" = roc_dataIn_B5R,
                       "B6" = roc_dataIn_B6),  "All")
plotrocAllCore <- 
  plotROCfunction(list("A5L" = roc_dataIn_A5L,
                       "ATI-C" = roc_dataIn_ATI_C, 
                       "ATI-N" = roc_dataIn_ATI_N),  "All")

plotrocAllDelta <- 
  plotROCfunction(list("Delta" = roc_dataIn_Delta,
                       "Combined" = roc_Combined),"All")

plotrocRatiosAll <- 
  plotROCfunction(list("A29/A27L" = roc_ratio_A29_A27L,
                       "A35R/A33R" = roc_ratio_A35R_A33R,
                       "B6/B5R" = roc_ratio_B6_B5R,
                       "E8/D8L" = roc_ratio_E8_D8L),"All")


plotrocCombined <-
  ggarrange(plotrocAllIMV, plotrocAllEEV, plotrocAllCore, plotrocAllDelta, plotrocRatiosAll,
            ncol = 5, nrow = 1, align = "hv",
            labels = c("a", "", "", "", ""))


ggsave("output/plotrocCombined.png", plotrocCombined, width = 12, height = 3, dpi = 600)
