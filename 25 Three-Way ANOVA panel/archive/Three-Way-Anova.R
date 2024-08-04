####
# Compare different panels by Three-Way-ANOVA
# Daniel Stern, RKI
# 2024-04-12
# Version 1.0
####

# Reference to method used: https://www.datanovia.com/en/lessons/anova-in-r/


# Clean environment
rm(list = ls(all.names = TRUE))

# Load packages   
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(rstatix)
# library(emmeans)

# Load data together with metadata 
load("../16 Compare Panels/output/dataInputPanelANOVA.Rdata")

# Prepare and select IgG dataframe and relevant data
dataInIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
 # filter(panel == "Spox") %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  select(panel, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         analyte, dataIn) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  filter(isotype %in% c("IgG")) %>% 
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
outliersExtreme <-
  outliers %>% 
  filter(is.outlier == TRUE) %>% 
  select(sampleID_metadata) %>% 
  pull()

dataInIgGnonOut <- 
  dataInIgG %>% 
  filter(!(sampleID_metadata %in% outliersExtreme)) %>% 
  filter(childhoodImmuAge != "Ambiguous")


# Groups for Three-Way-ANOVA:
# panel: Combine date from acute and epi panel in a first run
# Then repeat analysis for acute and epi panel separate
# # Var 1: panel
# Var 2: status
# Var 3: chilhoodImmuAge
# Var 4: analyte
# # Var 5: isotype: analyse IgG and IgM seperately

set.seed(123)

# Compute summary statistics (mean and sd) of dataInNorm by groups
summaryStatisticsAll <- 
  dataInIgGnonOut %>% 
  group_by(status, childhoodImmuAge, analyte) %>% 
  get_summary_stats(dataInNorm, type = "mean_sd")

# Visualize results as box plots
boxplotAll <-
  ggboxplot(dataInIgGnonOut, x = "status", y = "dataInNorm", color = "childhoodImmuAge",
            facet.by = "analyte") +
  scale_color_manual(name = "Smallpox vaccination", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# Check normality assumption
model  <- lm(dataInNorm ~ status*childhoodImmuAge*analyte, data = dataInIgGnonOut)
# Create a QQ plot of residuals
ggqqplot(residuals(model)) +
  theme_bw()
# Compute Shapiro-Wilk test of normality
shapiro_test(sample(residuals(model), 5000))

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

ggsave("output/SupFig_qqplotIgG.png", qqplotAll, width = 18, height = 14, dpi = 600)


# homogneity of variance assumption
dataInIgGnonOut %>% levene_test(dataInNorm ~ status*childhoodImmuAge*analyte)


# Calculate anova test
res.aov <- dataInIgGnonOut %>% anova_test(dataInNorm ~ status*childhoodImmuAge*analyte)

# There was a statistically significant three-way interaction between status, childhoodImmuAge and 
# analyte F(56, 19515) = 2.95


interaction.plot(dataInIgGnonOut$status, dataInIgGnonOut$childhoodImmuAge,
                 dataInIgGnonOut$dataInNorm)

# Generate interaction plot based on summary statistics all
summaryStatisticsAll %>% 
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

twoWayIgG <-
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
simpleSimpleMainStatusIgG <-
  dataInIgGnonOut %>%
  group_by(childhoodImmuAge, status) %>%
  anova_test(dataInNorm ~ analyte)

simpleSimpleMainAnalyeIgG <-
  dataInIgGnonOut %>%
  group_by(childhoodImmuAge, analyte) %>%
  anova_test(dataInNorm ~ status)

# Compute simple simple comparisons
pwcIgG <- dataInIgGnonOut %>%
  group_by(childhoodImmuAge, analyte) %>%
  emmeans_test(dataInNorm ~ status)# %>%
#select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk
meanDiffIgG <-
  get_emmeans(pwcIgG)

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
pwcIgG <- pwcIgG %>% add_xy_position(x = "status")
# pwcIgG.filtered <- pwcIgG %>% filter(gender == "male", risk == "high")
boxPlotPairwiseIgG <-
  boxplotAll +
  stat_pvalue_manual(
    pwcIgG, color = "childhoodImmuAge", hide.ns = TRUE,
    tip.length = 0.025, step.increase = 0.05, step.group.by = "analyte"
  ) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwcIgG)
  )




####
# Repeat analysis for IgM data
# Prepare and select IgM dataframe and relevant
dataInIgM <-
  dataInput %>% 
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
  filter(isotype %in% c("IgM")) %>% 
  rename(status = panel_detail) %>% 
  group_by(analyte) %>% 
  mutate(dataInNorm = (dataIn-min(dataIn))/(max(dataIn)-min(dataIn))) %>% #min max normalization for each antigen
  ungroup()


# Check for outliers
outliersIgM <- 
  dataInIgM %>% 
  group_by(status, childhoodImmuAge, analyte) %>% 
  identify_outliers(dataInNorm)

# Remove all outliers
outliersExtremeIgM <-
  outliers %>% 
  filter(is.outlier == TRUE) %>% 
  select(sampleID_metadata) %>% 
  pull()

dataInIgMnonOut <- 
  dataInIgM %>% 
  filter(!(sampleID_metadata %in% outliersExtremeIgM)) %>% 
  filter(childhoodImmuAge != "Ambiguous")




# Groups for Three-Way-ANOVA:
# panel: Combine date from acute and epi panel in a first run
# Then repeat analysis for acute and epi panel separate
# # Var 1: panel
# Var 2: status
# Var 3: chilhoodImmuAge
# Var 4: analyte
# # Var 5: isotype: analyse IgG and IgM seperately

set.seed(123)

# Compute summary statistics (mean and sd) of dataInNorm by groups
summaryStatisticsAllIgM <- 
  dataInIgMnonOut %>% 
  group_by(status, childhoodImmuAge, analyte) %>% 
  get_summary_stats(dataInNorm, type = "mean_sd")

# Visualize results as box plots
boxplotAllIgM <-
  ggboxplot(dataInIgMnonOut, x = "status", y = "dataInNorm", color = "childhoodImmuAge",
            facet.by = "analyte") +
  scale_color_manual(name = "Smallpox vaccination", values = colorblind_pal()(8)[c(2:8)]) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Min Max norm.") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# Check normality assumption
modelIgM  <- lm(dataInNorm ~ status*childhoodImmuAge*analyte, data = dataInIgMnonOut)
# Create a QQ plot of residuals
ggqqplot(residuals(modelIgM))
# Compute Shapiro-Wilk test of normality
shapiro_test(sample(residuals(modelIgM), 5000))

# Normality assumption is violated for all data

# Check normality assumptions by groups
normalityGroupedIgM <-
  dataInIgMnonOut %>% 
  group_by(status, childhoodImmuAge, analyte) %>% 
  shapiro_test(dataInNorm)  %>% 
  mutate(normDist = if_else(p >= 0.05, "normal distributed", "deviation"))


distributionIgM <-
  normalityGroupedIgM %>% 
  count(normDist)

# Violated for ~ 80 sets, normality for  55 sets
qqplotAllIgM <- 
  dataInIgMnonOut %>% 
  left_join(normalityGroupedIgM, by = c("status", "childhoodImmuAge", "analyte")) %>% 
  ggqqplot("dataInNorm", ggtheme = theme_bw(), color = "normDist") +
  facet_grid(childhoodImmuAge + status ~ analyte, labeller = "label_both")+
  scale_y_continuous(name = "Min Max norm.") +
  guides(fill = "none") +
  scale_color_manual(name = "Normal distribution", values = colorblind_pal()(8)[c(2:8)]) +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


ggsave("output/SupFig_qqplotIgM.png", qqplotAllIgM, width = 18, height = 14, dpi = 600)


# homogneity of variance assumption
dataInIgMnonOut %>% levene_test(dataInNorm ~ status*childhoodImmuAge*analyte)


# Calculate anova test
res.aovIgM <- dataInIgMnonOut %>% anova_test(dataInNorm ~ status*childhoodImmuAge*analyte)
summary(res.aovIgM)
interaction.plot(dataInIgMnonOut$status, dataInIgMnonOut$childhoodImmuAge,
                 dataInIgMnonOut$dataInNorm)

# Generate interaction plot based on summary statistics all
summaryStatisticsAllIgM %>% 
  filter(childhoodImmuAge != "Ambiguous") %>% 
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

twoWayIgM <-
  dataInIgMnonOut %>%
  group_by(childhoodImmuAge) %>%
  anova_test(dataInNorm ~ status*analyte)
# There was a statistically significant simple two-way interaction between
# status and analyse irrespective of the childhoodImmuAge status. However, the
# F value was largest in the unvaccinated group

# Compute simple simple main effects
# Investigate the effect of status at every level of analyte and 
# investigate the effect of analyte at every level of status on dataInNorm
# Group the data by childHoodImmuAge and status, and fit anova
simpleSimpleMainStatusIgM <-
  dataInIgMnonOut %>%
  group_by(childhoodImmuAge, status) %>%
  anova_test(dataInNorm ~ analyte)

simpleSimpleMainAnalyeIgM <-
  dataInIgMnonOut %>%
  group_by(childhoodImmuAge, analyte) %>%
  anova_test(dataInNorm ~ status)

# Compute simple simple comparisons
pwcIgM <- dataInIgMnonOut %>%
  group_by(childhoodImmuAge, analyte) %>%
  emmeans_test(dataInNorm ~ status)# %>%
#select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk
meanDiffIgM <-
  get_emmeans(pwcIgM)

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
pwcIgM <- pwcIgM %>% add_xy_position(x = "status")
# pwcIgG.filtered <- pwcIgG %>% filter(gender == "male", risk == "high")
boxPlotPairwiseIgM <-
  boxplotAllIgM +
  stat_pvalue_manual(
    pwcIgM, color = "childhoodImmuAge", hide.ns = TRUE,
    tip.length = 0.025, step.increase = 0.05, step.group.by = "analyte"
  ) +
  labs(
    subtitle = get_test_label(res.aovIgM, detailed = TRUE),
    caption = get_pwc_label(pwcIgM)
  )









