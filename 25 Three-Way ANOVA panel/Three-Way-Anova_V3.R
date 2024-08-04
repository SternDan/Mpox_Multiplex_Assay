####
# Compare different panels by Three-Way-ANOVA
# Daniel Stern, RKI
# 2024-04-12
# Version 3.0
####

# Reference to method used: https://www.datanovia.com/en/lessons/anova-in-r/


# Clean environment
rm(list = ls(all.names = TRUE))

# Load packages   
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(rstatix)
library(lme4)
library(ggeffects)

# Load data together with metadata 
load("../16 Compare Panels/output/dataInputPanelANOVA.Rdata")

set.seed(123)

# Generate function for analysis of different panels and isotypes

anovaFunction <- function(panelIn, isotypeIn){
  outlist <- list()
  ####
  # Prepare and select IgG dataframe and relevant data
  dataInIgG <-
    dataInput %>% 
    filter(!is.na(panel_strat)) %>% 
    filter(!is.na(childhoodImmuAge)) %>% 
    filter(panel_detail != "SPox") %>% 
    filter(analyte != "VACV") %>% 
    mutate(panel = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
    filter(panel %in% panelIn) %>% 
    select(panel, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
    select(panel, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
           analyte, dataIn) %>% 
    mutate(panel_detail = factor(panel_detail), 
           isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
           childhoodImmuAge = factor(childhoodImmuAge,
                                     levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
    filter(isotype %in% isotypeIn) %>% 
    rename(status = panel_detail) %>% 
    group_by(analyte) %>% 
    mutate(dataInNorm = (dataIn-min(dataIn))/(max(dataIn)-min(dataIn))) %>% #min max normalization for each antigen
    ungroup()
  
  # Check for outliers
  outliers <- 
    dataInIgG %>% 
    group_by(status, childhoodImmuAge, analyte) %>% 
    identify_outliers(dataInNorm)
  
  # Remove all outliers
  outliers <-
    outliers %>% 
    filter(is.outlier == TRUE) %>% 
    select(sampleID_metadata) %>% 
    pull()
  
  dataInIgGnonOut <- 
    dataInIgG %>% 
    filter(!(sampleID_metadata %in% outliers)) %>% 
    filter(childhoodImmuAge != "Ambiguous")
  
  # Groups for Three-Way-ANOVA:
  # panel: Combine date from acute and epi panel in a first run
  # Then repeat analysis for acute and epi panel separate
  # # Var 1: panel
  # Var 2: status
  # Var 3: chilhoodImmuAge
  # Var 4: analyte
  # # Var 5: isotype: analyse IgG and IgM seperately
  
  # Compute summary statistics (mean and sd) of dataInNorm by groups
  summaryStatistics <- 
    dataInIgGnonOut %>% 
    group_by(status, childhoodImmuAge, analyte) %>% 
    get_summary_stats(dataInNorm, type = "mean_sd")
  
  # Visualize results as box plots
  boxplot <-
    ggboxplot(dataInIgGnonOut, x = "status", y = "dataInNorm", color = "childhoodImmuAge",
              facet.by = "analyte") +
    scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)[c(2:8)]) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Min Max norm.") +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Check normality assumption
  model  <- lm(dataInNorm ~ status*childhoodImmuAge*analyte, data = dataInIgGnonOut)
  # Create a QQ plot of residuals
  qqplotAlldata <- 
    ggqqplot(residuals(model)) +
    theme_bw()
  # Compute Shapiro-Wilk test of normality
  shapiroOut <- 
    shapiro_test(sample(residuals(model), 1000))
  
  # Normality assumption is violated for all data
  
  # Check normality assumptions by groups
  normalityGrouped <-
    dataInIgGnonOut %>% 
    group_by(status, childhoodImmuAge, analyte) %>% 
    shapiro_test(dataInNorm) %>% 
    mutate(normDist = if_else(p >= 0.05, "normal distributed", "deviation"))
  
  # Calculate ratio of normal distribute antigens and deviation
  distributionIgG <-
    normalityGrouped %>% 
    count(normDist)
  
  # Violated for ~ 80 sets, normality for  55 sets
  qqplotAll <- 
    dataInIgGnonOut %>% 
    left_join(normalityGrouped, by = c("status", "childhoodImmuAge", "analyte")) %>% 
    ggqqplot("dataInNorm", ggtheme = theme_bw(), color = "normDist") +
    facet_grid(childhoodImmuAge + status ~ analyte, labeller = "label_both")+
    scale_y_continuous(name = "Min Max norm.") +
    guides(fill = "none") +
    scale_color_manual(name = "Normal distribution", values = colorblind_pal()(8)[c(2:8)]) +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # ggsave("output/SupFig_qqplotIgG.png", qqplotAll, width = 18, height = 14, dpi = 600)
  
  # homogneity of variance assumption
  leveneOut <-
    dataInIgGnonOut %>% levene_test(dataInNorm ~ status*childhoodImmuAge*analyte)
  
  # Calculate anova test
  res.aov <- dataInIgGnonOut %>% anova_test(dataInNorm ~ status*childhoodImmuAge*analyte)
  
  # There was a statistically significant three-way interaction between status, childhoodImmuAge and 
  # analyte F(56, 19515) = 2.95
  
  
  # Generate interaction plot based on summary statistics all
  interactionPlot <-
    summaryStatistics %>% 
    ggplot(mapping = aes(x = status, y = mean, color = childhoodImmuAge)) +
    geom_line(aes(group = childhoodImmuAge)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5, alpha = 0.5) +
    facet_wrap("analyte") +
    theme_bw() +
    scale_color_manual(name = "Smallpox vaccination", values = colorblind_pal()(8)[c(2:8)]) +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ####
  # Implement post-hoc tests
  
  # Compute simple two-way interactions
  # In this case: evaluate effect of analyte*status on each level of childhoodImmuAge
  # Use group-specific error terms, as assumptions have been violated
  # Group the data by childhoodImmuAge and 
  # fit simple two-way interaction 
  # modelTwoWayIgG  <- lm(pain_score ~ gender*risk*treatment, data = headache)
  
  twoWay <-
    dataInIgGnonOut %>%
    group_by(childhoodImmuAge) %>%
    anova_test(dataInNorm ~ status*analyte)
  # There was a statistically significant simple two-way interaction between
  # status and analyse irrespective of the childhoodImmuAge status. However, the
  # F value was largest in the unvaccinated group
  
  # Compute simple simple main effects
  # Investigate the effect of status at every level of analyte and 
  # investigate the effect of analyte at every level of status on dataInNorm
  # Group the data by childHoodImmuAge and status, and fit anova
  simpleSimpleMainStatus <-
    dataInIgGnonOut %>%
    group_by(childhoodImmuAge, status) %>%
    anova_test(dataInNorm ~ analyte)
  
  simpleSimpleMainAnalyte <-
    dataInIgGnonOut %>%
    group_by(childhoodImmuAge, analyte) %>%
    anova_test(dataInNorm ~ status)
  
  # Compute simple simple comparisons
  pwc <- dataInIgGnonOut %>%
    group_by(childhoodImmuAge, analyte) %>%
    emmeans_test(dataInNorm ~ status)# %>%
  #select(-df, -statistic, -p) # Remove details
  # Show comparison results for male at high risk
  meanDiffIgG <-
    get_emmeans(pwc)
  
  # Report draft
  # A three-way ANOVA was conducted to determine the effects of a possible childhood immunization against smallpox, the serostatus and the analyte tested on 
  # the normlized data measured by the multiplex assay.
  # Residual analysis was performed to test for the assumptions of the three-way ANOVA.
  # Normality was assessed using Shapiro-Wilk’s normality test and homogeneity of variances was assessed by Levene’s test.
  # Residuals were normally distributed (p > 0.05) and there was homogeneity of variances (p > 0.05).
  # There was a statistically significant three-way interaction between gender, risk and treatment, F(2, 60) = 7.41, p = 0.001.
  # Statistical significance was accepted at the p < 0.025 level for simple two-way interactions and simple simple main effects. There was a statistically significant simple two-way interaction between risk and treatment for males, F(2, 60) = 5.2, p = 0.008, but not for females, F(2, 60) = 2.8, p = 0.065.
  # There was a statistically significant simple simple main effect of treatment for males at high risk of migraine, F(2, 60) = 14.8, p < 0.0001), but not for males at low risk of migraine, F(2, 60) = 0.66, p = 0.521.
  # All simple simple pairwise comparisons, between the different treatment groups, were run for males at high risk of migraine with a Bonferroni adjustment applied.
  # There was a statistically significant mean difference between treatment X and treatment Y. However, the difference between treatment Y and treatment Z, was not statistically significant.
  
  # Visualization: box plots with p-values
  pwc <- pwc %>% add_xy_position(x = "status")
  # pwcIgG.filtered <- pwcIgG %>% filter(gender == "male", risk == "high")
  boxPlotPairwise <-
    boxplot +
    stat_pvalue_manual(
      pwc, color = "childhoodImmuAge", hide.ns = TRUE,
      tip.length = 0.025, step.increase = 0.05, step.group.by = "analyte"
    ) #+
  # labs(
  #    subtitle = get_test_label(res.aov, detailed = TRUE),
  #    caption = get_pwc_label(pwc)
  #)
  
  # Generate Output as named List
  # 1) datainput
  # 2) dataInIgGnonOut
  # 3) summaryStatistics
  # 4) qqplotAlldata
  # 5) shapiroOut
  # 6) normalityGrouped
  # 7) distributionIgG
  # 8) qqplotAll
  # 9) leveneOut
  # 10) res.aov
  # 11) interactionPlot
  # 12) twoWay
  # 13) simpleSimpleMainStatus
  # 14) simpleSimpleMainAnalyte
  # 15) pwc
  # 16) meanDiffIgG
  # 17) boxPlotPairwise
  outlist$dataInput <- dataInput
  outlist$dataInIgGnonOut <- dataInIgGnonOut
  outlist$summaryStatistics <- summaryStatistics
  outlist$qqplotAlldata <- qqplotAlldata
  outlist$shapiroOut <- shapiroOut
  outlist$normalityGrouped <- normalityGrouped
  outlist$distributionIgG <- distributionIgG
  outlist$qqplotAll <- qqplotAll
  outlist$leveneOut <- leveneOut
  outlist$res.aov <- res.aov
  outlist$interactionPlot <- interactionPlot
  outlist$twoWay <- twoWay
  outlist$simpleSimpleMainStatus <- simpleSimpleMainStatus
  outlist$simpleSimpleMainAnalyte <- simpleSimpleMainAnalyte
  outlist$pwc <- pwc
  outlist$meanDiffIgG <- meanDiffIgG
  outlist$boxPlotPairwise <- boxPlotPairwise
  return(outlist)
}

anovaResultsAllIgG <-
  anovaFunction(c("Acute", "Spox"), c("IgG"))

anovaResultsAllIgM <-
  anovaFunction(c("Acute", "Spox"), c("IgM"))

anovaResultsAcuteIgG <-
  anovaFunction(c("Acute"), c("IgG"))

anovaResultsAcuteIgM <-
  anovaFunction(c("Acute"), c("IgM"))

anovaResultsSpoxIgG <-
  anovaFunction(c("Spox"), c("IgG"))

anovaResultsSpoxIgM <-
  anovaFunction(c("Spox"), c("IgM"))


# qqplot All data 
plotqqAll <-
  ggarrange(
    anovaResultsAllIgG$qqplotAlldata,
    anovaResultsAllIgM$qqplotAlldata,
    anovaResultsAcuteIgG$qqplotAlldata,
    anovaResultsAcuteIgM$qqplotAlldata,
    anovaResultsSpoxIgG$qqplotAlldata,
    anovaResultsSpoxIgM$qqplotAlldata, ncol = 2, nrow = 3, align = "hv",
    labels = "auto")

ggsave("output/SFig_plotqqAll.png", plotqqAll, width = 8, height = 6, 
       dpi = 600)

plotqqAcute <-
  ggarrange(
    anovaResultsAcuteIgG$qqplotAll,
    anovaResultsAcuteIgM$qqplotAll, nrow = 2, align = "hv",
    labels = "auto", common.legend = TRUE)

ggsave("output/SFig_plotqqAcute.png", plotqqAcute, width = 16, height = 20, 
       dpi = 600)


plotqqSpox <-
  ggarrange(
    anovaResultsSpoxIgG$qqplotAll,
    anovaResultsSpoxIgM$qqplotAll, nrow = 2, align = "hv",
    labels = "auto", common.legend = TRUE)

ggsave("output/SFig_plotqqSpox.png", plotqqSpox, width = 16, height = 20, 
       dpi = 600)


boxplotsAnovaAll <- 
  ggarrange(anovaResultsAllIgG$boxPlotPairwise, anovaResultsAllIgM$boxPlotPairwise,
            ncol = 2, align = "hv", labels = "auto", common.legend = TRUE)

boxplotsAnovaAcute <- 
  ggarrange(anovaResultsAcuteIgG$boxPlotPairwise, anovaResultsAcuteIgM$boxPlotPairwise,
            ncol = 2, align = "hv", labels = "auto", common.legend = TRUE)

boxplotsAnovaSpox <- 
  ggarrange(anovaResultsSpoxIgG$boxPlotPairwise, anovaResultsSpoxIgM$boxPlotPairwise,
            ncol = 2, align = "hv", labels = "auto", common.legend = TRUE)



ggsave("output/SFig_boxplotsAnovaAll.png", boxplotsAnovaAll, width = 12, height = 8,
       dpi = 600)

ggsave("output/SFig_boxplotsAnovaAcute.png", boxplotsAnovaAcute, width = 12, height = 8,
       dpi = 600)

ggsave("output/SFig_boxplotsAnovaSpox.png", boxplotsAnovaSpox, width = 12, height = 8,
       dpi = 600)

## Generate qqplots
boxplotsAnovaAcute <- 
  ggarrange(anovaResultsAcuteIgG$boxPlotPairwise, anovaResultsAcuteIgM$boxPlotPairwise,
            ncol = 2, align = "hv", labels = "auto")

boxplotsAnovaCombined <- 
  ggarrange(anovaResultsAcuteIgG$boxPlotPairwise, anovaResultsAcuteIgM$boxPlotPairwise,
            anovaResultsSpoxIgG$boxPlotPairwise, anovaResultsSpoxIgM$boxPlotPairwise,
            ncol = 2, nrow = 2, align = "hv", common.legend = TRUE,
            labels = c("c", "d", "e", "f"))

save(boxplotsAnovaCombined, file = "output/boxplotsAnova.Rdata")

## 
# Generate output for results of ANOVA analyis

##
# Generate results for Three Way ANOVA
anovaAcuteIgG <- anovaResultsAcuteIgG$res.aov
anovaEpiIgG <- anovaResultsSpoxIgG$res.aov
anovaAcuteIgM <- anovaResultsAcuteIgM$res.aov
anovaEpiIgM <- anovaResultsSpoxIgM$res.aov


# Generate out with results for pairwise comparison
simplesimpleAcuteIgG <- anovaResultsAcuteIgG$simpleSimpleMainAnalyte
simplesimpleEpiIgG <- anovaResultsSpoxIgG$simpleSimpleMainAnalyte
simplesimpleAcuteIgM <- anovaResultsAcuteIgM$simpleSimpleMainAnalyte
simplesimpleEpiIgM <- anovaResultsSpoxIgM$simpleSimpleMainAnalyte

# Generate results for pairwise comparison
pairwiseAcuteIgG <- anovaResultsAcuteIgG$pwc  %>% 
  select(c(1:11)) %>% 
  select(-term, -.y.)
pairwiseEpiIgG <- anovaResultsSpoxIgG$pwc %>% 
  select(c(1:11)) %>% 
  select(-term, -.y.)
pairwiseAcuteIgM <- anovaResultsAcuteIgM$pwc %>% 
  select(c(1:11)) %>% 
  select(-term, -.y.)
pairwiseEpiIgM <- anovaResultsSpoxIgM$pwc %>% 
  select(c(1:11)) %>% 
  select(-term, -.y.)

save(anovaAcuteIgG, anovaEpiIgG, anovaAcuteIgM, anovaEpiIgM,
     simplesimpleAcuteIgG, simplesimpleEpiIgG, simplesimpleAcuteIgM, simplesimpleEpiIgM,
     pairwiseAcuteIgG, pairwiseEpiIgG, pairwiseAcuteIgM, pairwiseEpiIgM,
     file = "output/tables.RData" )

##
# Fit a glm model to analyse the impact of childhoodImmu on the normalized data
# in dependence of the serostatus

####
# Prepare and select IgG dataframe and relevant data
dataInNorm <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  #  filter(panel %in% panelIn) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         analyte, dataIn) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  #  filter(isotype %in% isotypeIn) %>% 
  rename(status = panel_detail) %>% 
  group_by(analyte) %>% 
  mutate(dataInNorm = (dataIn-min(dataIn))/(max(dataIn)-min(dataIn))) %>% #min max normalization for each antigen
  ungroup()

# Check for outliers
outliers <- 
  dataInNorm %>% 
  group_by(status, childhoodImmuAge, analyte) %>% 
  identify_outliers(dataInNorm)

# Remove all outliers
outliers <-
  outliers %>% 
  filter(is.outlier == TRUE) %>% 
  select(sampleID_metadata) %>% 
  pull()

dataInNormnonOutIgG <- 
  dataInNorm %>% 
  filter(!(sampleID_metadata %in% outliers)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(isotype == "IgG")


modelIgG <- lmer(dataInNorm ~ childhoodImmuAge + (1|analyte) + (1|sampleID_metadata) + (1|status) + (1|panel), data = dataInNormnonOutIgG)
summary(modelIgG)
plot(modelIgG)
qqnorm(resid(modelIgG))
qqline(resid(modelIgG)) 

predIgG <- ggpredict(modelIgG, terms = c("childhoodImmuAge"))
plotIgG <-
  ggplot(predIgG) + 
  geom_violin(data = dataInNormnonOutIgG,                      # adding the raw data (scaled values)
              aes(x = childhoodImmuAge, y = dataInNorm), width = 0.75) +
  geom_line(aes(x = x, y = predicted, group = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error,
                  group = group), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  stat_compare_means(data = dataInNormnonOutIgG,                      # adding the raw data (scaled values)
                     aes(x = childhoodImmuAge, y = dataInNorm),
                     method = "wilcox.test") +
  labs(x = "Smallpox vaccination", y = "Min Max norm.") + 
  theme_pubr()


# Repeat analysis for IgM
dataInNormnonOutIgM <- 
  dataInNorm %>% 
  filter(!(sampleID_metadata %in% outliers)) %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
  filter(isotype == "IgM")


modelIgM <- lmer(dataInNorm ~ childhoodImmuAge + (1|analyte) + (1|sampleID_metadata) + (1|status) + (1|panel), data = dataInNormnonOutIgM)
summary(modelIgM)
plot(modelIgM)
qqnorm(resid(modelIgM))
qqline(resid(modelIgM)) 

predIgM <- ggpredict(modelIgM, terms = c("childhoodImmuAge"))

plotIgM <- 
  ggplot(predIgM) + 
  geom_violin(data = dataInNormnonOutIgM,                      # adding the raw data (scaled values)
              aes(x = childhoodImmuAge, y = dataInNorm), width = 0.75) + 
  geom_line(aes(x = x, y = predicted, group = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error,
                  group = group), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  stat_compare_means(data = dataInNormnonOutIgM,                      # adding the raw data (scaled values)
                     aes(x = childhoodImmuAge, y = dataInNorm),
                     method = "wilcox.test") +
  labs(x = "Smallpox vaccination", y = "Min Max norm.") + 
  theme_pubr()
plotMeanIgGIgM <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = "auto")

ggsave("output/SFig_compMeans.png", plotMeanIgGIgM, width = 6, height = 3, dpi = 600)




# Model over panel
modelIgMPanel <- lmer(dataInNorm ~ childhoodImmuAge + panel + (1|analyte) + (1|sampleID_metadata) + (1|status), data = dataInNormnonOutIgM)
summary(modelIgMPanel)
plot(modelIgM)
qqnorm(resid(modelIgM))
qqline(resid(modelIgM)) 

predIgMPanel <- ggpredict(modelIgMPanel, terms = c("childhoodImmuAge", "panel"))

#plotIgMPanel <- 
  ggplot(predIgMPanel) + 
  geom_violin(data = dataInNormnonOutIgM,                      # adding the raw data (scaled values)
              aes(x = childhoodImmuAge, y = dataInNorm, fill = panel), width = 0.75) + 
  geom_line(aes(x = x, y = predicted, group = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, color = group), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  stat_compare_means(data = dataInNormnonOutIgM,                      # adding the raw data (scaled values)
                     aes(x = childhoodImmuAge, y = dataInNorm),
                     method = "wilcox.test") +
  labs(x = "Smallpox vaccination", y = "Min Max norm.") + 
  theme_pubr()
plotMeanIgGIgM <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = "auto")

ggsave("output/SFig_compMeans.png", plotMeanIgGIgM, width = 6, height = 3, dpi = 600)


