##
# Compare serological panel seroconversion
# Daniel Stern
# Robert Koch Institute
# ZBS 3
# Version 1.2
# Last modified 2024-03-05
##
rm(list = ls(all.names = TRUE))

library(rio)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggforce)
library(GGally)
library(brms)
library(tidybayes)
library(see)



##
# Load metadata for Age-groups
## Load metadata
## Load datainput generated for clustering
load("../24 Spox Pre Pos Clustering/input/dataClustering.Rdata")

## Unify the format for the Age_groups varibale
dataInputMeta <- 
  dataClustering %>% 
  dplyr::filter(isotype != "IgA") %>% 
  mutate(Age_group = case_when(Age_group == "18-29" ~ "< 30",
                               Age_group == "30-39" ~ "< 40",
                               Age_group == "40-49" ~ "< 50",
                               Age_group == "50-59" ~ "< 60",
                               TRUE ~ Age_group),
         Age_group = factor(Age_group, levels = c("< 30", 
                                                  "< 40", 
                                                  "< 50",
                                                  "< 60", 
                                                  "60+"),
                            ordered = TRUE)) %>% 
  select(-starts_with("pred")) %>%
  select(sampleID_metadata, Age_group) %>% 
  unique()


##
# Load data: All data
dataInput <- # import("../15 Generate Dataframe ML/output/dataInputAll.csv") %>% 
  dataClustering %>% 
  mutate(serostatus_clear = if_else(serostatus_cat == "positive", 1, 0)) %>% 
  mutate(panel_strat = case_when(panel == "SPox" ~ "SPox",
                                 panel %in% c("MPXV", "MVA", "Pre", "Pre_New") ~ "Acute")) %>% 
  mutate(#analyte = factor(analyte, levels = c("A27L", "A29",
    #                                    "L1R", "M1", 
    #                                   "D8L", "E8",
    #                                  "H3L", "A33R", "A35R",
    #                                 "B5R", "B6", "A5L", "ATI-C", 
    #                                "ATI-N", "Delta", "VACV"), 
    #           labels = c("A27L", "A29L",
    #                     "L1R", "M1R", 
    #                    "D8L", "E8L",
    #                   "H3L", "A33R", "A35R",
    #                  "B5R", "B6R", "A5L", "ATI-C", 
    #                 "ATI-N", "Delta", "VACV"), ordered = TRUE),
    #  panel_strat = if_else(panel_detail == "Pre_New", "Pre_New", panel_strat), 
    panel_detail = if_else(panel == "Pre_New", "Pre", panel_detail),
    panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "MPXV", 
                                                   "SPox", "SPox_Rep", "CPXV"),
                          labels = c("Pre", "MVA", "Mpox", "SPox", "SPox_Rep", "CPXV"),
                          ordered = TRUE))  %>% 
  # left_join(dataInputMeta, by = c("sampleID_metadata")) %>% 
  mutate(childhoodImmuAge = case_when(Age_group %in% c("18-29", "30-39") ~ "No",
                                      Age_group %in% c("40-49") ~ "Ambiguous",
                                      Age_group %in% c("50-59", "60+") ~ "Yes",
                                      Age_group %in% c("< 30", "< 40") ~ "No",
                                      Age_group %in% c("< 50") ~ "Ambiguous",
                                      Age_group %in% c("< 60", "60+") ~ "Yes"))




unique(dataClustering$Age_group)
##
# Plot panels distribution as used for ML-algorithms
# New plot: parallel plot with three groups
## Stratify for childhoodimmu
plotIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotIgA <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgA")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         H3L, 
         A33R, 
         A35R,
         B5R, 
         B6R, 
         A5L, 
         `ATI-C`, 
         `ATI-N`, 
         Delta) %>% 
  mutate(panel_detail = factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE)) %>% 
  ggparcoord(columns = c(6:20), alphaLines = 0.3, groupColumn = "childhoodImmuAge",
             scale = "uniminmax", boxplot = F) +
  scale_color_manual(name = "Childhood Immunization", values = colorblind_pal()(8)) +
  theme_bw() +
  coord_polar() +
  facet_grid(panel_plot ~ panel_detail, labeller = as_labeller(c("1" = "Pre", "2" = "MVA",
                                                                 "3" = "Mpox",
                                                                 "Pre" = "Pre",
                                                                 "MVA" = "MVA", "Mpox" = "Mpox"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotStratifiedChildhood <-
  ggarrange(plotIgG, plotIgM, ncol = 2, align = "hv", labels = c("c", "d"),
            common.legend = TRUE)

ggsave("output/plotStratified.png", plotStratifiedChildhood, width = 14, height = 6,
       dpi = 600)

####
# Plot ratios of different antigens 
plotRatioIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         A33R, 
         A35R,
         B5R, 
         B6R) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1R/L1R` = M1R/L1R,
         `A29L/A27L` = A29L/A27L,
         `B6R/B5R` = B6R/B5R, 
         `A35R/A33R` = A35R/A33R,
         `E8L/D8L` = E8L/D8L) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1R/L1R`, `A29L/A27L`, `A35R/A33R`, `B6R/B5R`,  `E8L/D8L`) %>% 
  ggparcoord(columns = c(6:10), alphaLines = 0.1, groupColumn = "panel_detail",
             scale = "center", boxplot = F) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)) +
  theme_bw() +
  facet_grid(. ~ childhoodImmuAge, labeller = as_labeller(c("1" = "Yes", "2" = "No", "3" = "Ambiguous"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotRatioIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  pivot_wider(names_from = analyte, values_from = dataIn, values_fn = mean) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         A27L, 
         A29L,
         L1R, 
         M1R, 
         D8L, 
         E8L,
         A33R, 
         A35R,
         B5R, 
         B6R) %>% 
  mutate(panel_detail = as.factor(panel_detail), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         `M1R/L1R` = M1R/L1R,
         `A29L/A27L` = A29L/A27L,
         `B6R/B5R` = B6R/B5R, 
         `A35R/A33R` = A35R/A33R,
         `E8L/D8L` = E8L/D8L) %>% 
  filter(!is.na(childhoodImmuAge)) %>% 
  filter(!is.infinite(`M1R/L1R`)) %>% 
  filter(!is.infinite(`A29L/A27L`)) %>% 
  filter(!is.infinite(`B6R/B5R`)) %>% 
  filter(`M1R/L1R` < 3) %>% 
  filter(`A29L/A27L` < 3) %>% 
  filter(`B6R/B5R` < 3) %>% 
  filter(`A35R/A33R` < 3) %>% 
  filter(`E8L/D8L` < 3) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, childhoodImmuAge,
         `M1R/L1R`, `A29L/A27L`, `A35R/A33R`, `B6R/B5R`,  `E8L/D8L`) %>% 
  ggparcoord(columns = c(6:10), alphaLines = 0.1, groupColumn = "panel_detail",
             scale = "center", boxplot = F) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)) +
  theme_bw() +
  facet_grid(. ~ childhoodImmuAge, labeller = as_labeller(c("1" = "Yes", "2" = "No", "3" = "Ambiguous"))) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


plotRatios <-
  ggarrange(plotRatioIgG, plotRatioIgM, ncol = 2, align = "hv", labels = c("e", "f"),
            common.legend = TRUE)

ggsave("output/plotRatios.png", plotRatios, width = 14, height = 3,
       dpi = 600)


plotDensScaledIgG <- 
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  # filter(analyte == "ATI-N") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  #  filter(panel_plot %in% c("Acute")) %>% 
  #  filter(childhoodImmuAge == "Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "Mpox"),
                               ordered = TRUE), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         analyte = factor(analyte, levels = c("A27L", "A29L",
                                              "L1R", "M1R", 
                                              "D8L", "E8L",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6R", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), ordered = TRUE),) %>% 
  group_by(analyte) %>% 
  mutate(dataInMinMax = minMax(dataIn)) %>% 
  ggplot(mapping = aes(x = dataInMinMax, color = panel_detail, ..scaled..)) +
  geom_density(aes(fill = panel_detail), alpha = 0.1) +
  theme_bw() +
  scale_x_continuous(name = "Data", limits = c(0,1), labels = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(name = "Density", limits = c(0,1),  labels = c(0,1),  breaks = c(0,1)) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)[2:8]) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:8]) +
  facet_grid(panel_plot ~ analyte) +
  theme(strip.background = element_blank(),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotDensScaledIgM <- 
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  # filter(analyte == "ATI-N") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
  #  filter(panel_plot %in% c("Acute")) %>% 
  #  filter(childhoodImmuAge == "Yes") %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "Mpox"),
                               ordered = TRUE), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         analyte = factor(analyte, levels = c("A27L", "A29L",
                                              "L1R", "M1R", 
                                              "D8L", "E8L",
                                              "H3L", "A33R", "A35R",
                                              "B5R", "B6R", "A5L", "ATI-C", 
                                              "ATI-N", "Delta", "VACV"), ordered = TRUE),) %>% 
  group_by(analyte) %>% 
  mutate(dataInMinMax = minMax(dataIn)) %>% 
  ggplot(mapping = aes(x = dataInMinMax, color = panel_detail, ..scaled..)) +
  geom_density(aes(fill = panel_detail), alpha = 0.1) +
  theme_bw() +
  scale_x_continuous(name = "Data", limits = c(0,1), labels = c(0,1), breaks = c(0,1)) +
  scale_y_continuous(name = "Density", limits = c(0,1),  labels = c(0,1),  breaks = c(0,1)) +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)[2:8]) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:8]) +
  facet_grid(panel_plot ~ analyte) +
  theme(strip.background = element_blank(),
        #   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotDens <-
  ggarrange(plotDensScaledIgG, plotDensScaledIgM, nrow = 2, align = "hv", common.legend = TRUE, labels = c("a", "b"))


ggsave("output/plotDens.png", plotDens, width = 14, height = 5, dpi = 600)


plotPanelsCombined <- 
  ggarrange(plotDens, plotStratifiedChildhood, plotRatios, nrow = 3, align = "hv",
            common.legend = F, heights = c(2, 2,1.3))

ggsave("output/plotPanelsCombined.png", plotPanelsCombined, width = 12, height = 12,
       dpi = 600)

####
# Bayes statistical model to compare impact of different groups on 
# normalized binding response (min max for each antigen)


# custom function to implement min max scaling
minMax <- function(x) {
  (x - min(x)) / (max(x) - (min(x)))
}

dataInBayes <-  
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
 # filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(panel_plot = if_else(panel_strat %in% c("Acute", "Pre_New"), "Acute", "Spox")) %>% 
 # filter(panel_plot %in% c("Acute")) %>% 
  select(panel_plot, isotype, panel_detail, sampleID_metadata, analyte, dataIn, childhoodImmuAge) %>% 
  mutate(panel_detail = factor(panel_detail, levels = c("Pre", "MVA", "Mpox"),
                               ordered = TRUE), 
         isotype = factor(isotype, levels = c("IgG", "IgM"), ordered = TRUE),
         childhoodImmuAge = factor(childhoodImmuAge,
                                   levels = c("Yes", "No", "Ambiguous"), ordered = TRUE),
         analyte = str_replace(analyte, "-", "_")) %>% 
  group_by(analyte, isotype) %>% 
  mutate(dataInMinMax = minMax(dataIn),
         dataInMinMax = case_when(dataInMinMax == 0 ~ 0.001,
                                  dataInMinMax == 1 ~ 0.999,
                                  TRUE ~ dataInMinMax))



# Fit the Bayesian logistic regression model
# The formula specifies 'agree_fraq' as the binary outcome variable and 
# 'model', 'serostatus', 'antibody', and 'seropanel' as predictors
model_All <- brm(
  formula = dataInMinMax ~  analyte + panel_detail + childhoodImmuAge + panel_plot + isotype,
  data = dataInBayes ,
 # family = "beta", # gaussian distribution
  prior = c(
    set_prior("normal(0.5, 0.5)", class = "b")
  ),
  chains = 4, # Number of MCMC chains
  iter = 2000, # Number of MCMC iterations
  warmup = 1000, # Number of warmup iterations
  seed = 123 # Setting seed for reproducibility
)


# View summary of model
summary(model_ATI_NAcute)


summary_model_ATI_N <-
  summary(model_ATI_NAcute)

parameters_modelATI_N <- as.data.frame(summary_model_ATI_N$fixed) %>% 
  rownames_to_column("Parameters") %>% 
  mutate(Data = "Training All")

unique(dataInBayesIgG$panel_detail)

## Use tidybayes for extraction and plotting of the results
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

# Define functions
# Function to extract results and calculate median
getMedian <- function(inputModel, inputParam){
  inputModel %>%
    spread_draws((!!sym(inputParam))) %>%
    median_qi((!!sym(inputParam))) %>% 
    mutate(Parameter = str_remove(inputParam, "b_")) %>% 
    rename(Estimator = (!!sym(inputParam)),
           Lower = `.lower`,
           Upper = `.upper`) %>% 
    select(Parameter, Estimator, Lower, Upper)
}

# Function to extract all results for plotting of distributions
getRaw <- function(inputModel, inputParam){
  inputModel %>%
    spread_draws((!!sym(inputParam))) %>%
    # median_qi((!!sym(inputParam))) %>% 
    mutate(Parameter = str_remove(inputParam, "b_")) %>% 
    rename(Estimator = (!!sym(inputParam)),
           Chain = `.chain`,
           Iteration = `.iteration`,
           Draw = `.draw`)
}


# Define function to generate plot based on different datainputs
generate_var_plot <- function(dataInput){
  variables <-
    get_variables(dataInput)[grep("^b_", get_variables(dataInput))]
  
  outlist <- list()
  for(i in c(1:length(variables))){
    outlist[[i]] <- getMedian(dataInput, variables[i])
  }
  
  medianModelChildhoodImmun <- 
    do.call(rbind.data.frame, outlist)
  
  
  outlistRaw <- list()
  for(i in c(1:length(variables))){
    outlistRaw[[i]] <- getRaw(dataInput, variables[i])
  }
  
  rawModelChildhoodImmun <- 
    do.call(rbind.data.frame, outlistRaw)
  
  
  rawModelChildhoodImmunIntercept <-
    rawModelChildhoodImmun %>% 
    filter(Parameter == "Intercept")
  
  rawmodelCombined <-
    rawModelChildhoodImmun %>% 
    left_join(rawModelChildhoodImmunIntercept, 
              by = c("Chain", "Iteration", "Draw"), suffix = c("", "_inter")) %>% 
    mutate(Combined = if_else(Parameter == "Intercept", Estimator,
                              Estimator + Estimator_inter),
           Probability = exp(Combined))
  outPlot <-
    rawmodelCombined %>% 
    # filter(!(Parameter %in% c("seropanelallacute", "seropanelallepi"))) %>% 
    mutate(Factor = case_when(#grepl("model", Parameter) ~ "Model",
      grepl("childhood", Parameter) ~ "Smallpox vacc probable",
      grepl("detail", Parameter) ~ "Serostatus",
      grepl("plot", Parameter)  ~ "Seropanel",
      # grepl("delta", Parameter) ~ "Serostatus",
      grepl("Intercept", Parameter) ~ "Baseline"),
      Factor = factor(Factor, levels = c("Baseline","Serostatus",  "Seropanel", 
                                         "Smallpox vacc probable"), ordered = TRUE),
      Parameter = factor(Parameter, levels = c("Intercept",
                                               "panel_detail.L",
                                               "panel_detail.Q",
                                                "panel_plotSpox",
                                               "childhoodImmuAge.L",
                                               "childhoodImmuAge.Q"),
                         labels = c("Intercept (Pre/Acute/Vacc)",
                                    "Panel Linear",
                                    "Panel Quadratic",
                                    "Spox",
                                    "Vacc Linear",
                                    "Vacc Quadratic"), 
                         ordered = TRUE)) %>% 
    ggplot(aes(x = Parameter, y = Probability, fill = Factor)) +
    geom_violinhalf(scale = "width") +
    stat_halfeye() +
    scale_x_discrete() +
    scale_y_continuous(name = "Probability Correct Prediction (%)") +
    scale_fill_manual(values = colorblind_pal()(8)[2:8]) +
    # coord_flip() +
    geom_hline(yintercept = as.numeric(medianModelChildhoodImmun[1,2]), linetype = "dashed") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(outPlot)
}
#plot_ATI_N <-
generate_var_plot(model_ATI_NAcute)
sort(levels(dataInBayesIgG$panel_detail))






# Plot 1: IgG Results
plotMLpanelAcuteIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat %in% c("Acute", "Pre_New")) %>% 
  # filter(panel_st != "SPox") %>% 
  ggplot(mapping = aes(x =panel_detail , y = dataIn, fill = panel_detail)) +
  
  #  geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  #  geom_violin() +
  geom_boxplot() +
  # geom_jitter(width = 0.1, alpha = 0.2) +
  #  geom_boxplot(outlier.colour = NA, fill = "white", width = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 6) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 5) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,7)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("analyte") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Acute")


plotMLpanelEpiIgG <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat == "SPox") %>% 
  # filter(panel_st != "SPox") %>% 
  ggplot(mapping = aes(x =panel_detail , y = dataIn, fill = panel_detail)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 6) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 5) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,7)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("analyte") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Epi")



# Plot 2: IgM Results
plotMLpanelAcuteIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat %in% c("Acute", "Pre_New")) %>% 
  # filter(panel_st != "SPox") %>% 
  ggplot(mapping = aes(x =panel_detail , y = dataIn, fill = panel_detail)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 6) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 5) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,7)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("analyte") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Acute")

plotMLpanelEpiIgM <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat == "SPox") %>% 
  # filter(panel_st != "SPox") %>% 
  ggplot(mapping = aes(x =panel_detail , y = dataIn, fill = panel_detail)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 6) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 5) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0,7)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("analyte") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Epi")


plotMLpanelsCombined <- 
  ggarrange(plotMLpanelAcuteIgG, plotMLpanelEpiIgG,
            #  plotMLpanelAcuteIgM, plotMLpanelEpiIgM,
            #   ncol = 2,
            ncol = 2, align = "hv")
plotMLpanelsCombinedIgM <- 
  ggarrange(plotMLpanelAcuteIgM, plotMLpanelEpiIgM,
            #  plotMLpanelAcuteIgM, plotMLpanelEpiIgM,
            #   ncol = 2,
            ncol = 2, align = "hv")




##
# Plot ratios for homologue antigens
plotRatiosHomologue <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte %in% c("A27L", "A29L",
                        "L1R", "M1R", 
                        "D8L", "E8L",
                        "A33R", "A35R", 
                        "B5R", "B6R")) %>% 
  # filter(panel_strat == "Acute") %>% 
  # filter(panel_detail != "Pre") %>% 
  mutate(homologue = if_else(analyte %in% c("A27L",
                                            "L1R", 
                                            "D8L", 
                                            "A33R",
                                            "B5R"), "VACV", "MPXV"),
         antigenPair = case_when(analyte %in% c("A27L", "A29L") ~ "A27L/A29L",
                                 analyte %in% c("L1R", "M1R") ~ "L1R/M1R",
                                 analyte %in% c("D8L", "E8L") ~ "D8L/E8L",
                                 analyte %in% c("A33R", "A35R") ~ "A33R/A35R",
                                 analyte %in% c("B5R", "B6R") ~ "B5R/B6")) %>% 
  select(-analyte, -experiment, -plate_assay, -batch, -data, -serostatus, -serostatus_cat,
         serostatus.delta, serostatus_cat.delta, -serostatus_clear) %>% 
  pivot_wider(names_from = homologue, values_from = dataIn) %>% 
  ggplot(mapping = aes(x = VACV, y = MPXV, color = panel_detail)) +
  geom_point(alpha = 0.5) +
  # geom_mark_ellipse(aes( fill =panel_detail), 
  #                    expand = unit(1, "mm"), alpha = 0.2) +
  facet_grid(antigenPair ~ isotype) +
  scale_x_continuous(name = "Multiplex results (VACV proteins)") +
  scale_y_continuous(name = "Multiplex results (MPXV proteins)") +
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  theme_pubr()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



##
# Novel plot ratios -> Instead od scatter, plot ratios as box-plots
# with ratios of 1 as geom_vline
# -> mutate: calclate ratios
# -> plot using the same color sceme and arrangement as in comparison 
# Acute and Epi

dataRatios <- 
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype %in% c("IgG", "IgM")) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte %in% c("A27L", "A29L",
                        "L1R", "M1R", 
                        "D8L", "E8L",
                        "A33R", "A35R", 
                        "B5R", "B6R")) %>% 
  # filter(panel_strat == "Acute") %>% 
  # filter(panel_detail != "Pre") %>% 
  mutate(homologue = if_else(analyte %in% c("A27L",
                                            "L1R", 
                                            "D8L", 
                                            "A33R",
                                            "B5R"), "VACV", "MPXV"),
         antigenPair = case_when(analyte %in% c("A27L", "A29L") ~ "A27L/A29L",
                                 analyte %in% c("L1R", "M1R") ~ "L1R/M1R",
                                 analyte %in% c("D8L", "E8L") ~ "D8L/E8L",
                                 analyte %in% c("A33R", "A35R") ~ "A33R/A35R",
                                 analyte %in% c("B5R", "B6R") ~ "B5R/B6")) %>% 
  select(-analyte, -experiment, -plate_assay, -batch, -data, -serostatus, -serostatus_cat,
         serostatus.delta, serostatus_cat.delta, -serostatus_clear) %>% 
  pivot_wider(names_from = homologue, values_from = dataIn) %>% 
  mutate(ratio = MPXV/VACV)

plotRatiosEpiBoxplotIgG <-
  dataRatios %>% 
  filter(isotype == "IgG") %>% 
  filter(panel_detail != "SPox") %>% 
  filter(panel_strat %in% c("SPox")) %>% 
  ggplot(mapping = aes(x =panel_detail , y = ratio, fill = panel_detail)) +
  # geom_beeswarm(aes(color = panel_detail), alpha = 0.2) +
  #geom_violin() +
  geom_boxplot() +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 2.5) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 2) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0, 3)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("antigenPair") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Ratios Epi")



plotRatiosAcuteBoxplotIgG <-
  dataRatios %>% 
  filter(isotype == "IgG") %>% 
  filter(panel_detail != "SPox") %>% 
  filter(panel_strat %in% c("Acute", "Pre_New")) %>% 
  ggplot(mapping = aes(x =panel_detail , y = ratio, fill = panel_detail)) +
  # geom_beeswarm(aes(color = panel_detail), alpha = 0.2) +
  #geom_violin() +
  geom_boxplot() +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 2.5) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 2) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)", limits = c(0, 3)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  scale_color_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("antigenPair") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Ratios Acute")




plotRatiosBoxplotIgM <-
  dataRatios %>% 
  filter(isotype == "IgM") %>% 
  ggplot(mapping = aes(x =panel_detail , y = ratio, fill = panel_detail)) +
  geom_boxplot() +
  # geom_point(position = position_dodge(width = 0.75
  #  ), alpha = 0.1) +
  stat_compare_means(comparisons = list(c("Mpox", "MVA"),
                                        c("MVA", "Pre")),
                     label = "p.signif", label.y = 2.5) +
  stat_compare_means(comparisons = list(c("Mpox", "Pre")),
                     label = "p.signif", label.y = 2) +
  theme_pubr()+
  scale_y_continuous(name = "Multiplex (log10 VIG µg/mL)") +
  coord_cartesian(ylim = c(0,3)) +
  scale_fill_manual(name = "Panel", values = colorblind_pal()(8)[2:7]) + 
  facet_wrap("antigenPair") +
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  xlab("") +
  labs(title = "Ratios IgM")

plotRatiosIgG <-
  ggarrange(plotRatiosAcuteBoxplotIgG, plotRatiosEpiBoxplotIgG, nrow = 2, align = "hv")





plotFig4 <-
  ggarrange(plotMLpanelsCombined, plotRatiosIgG,
            ncol = 2, widths = c(2,1),
            labels = c("b", "c"))





ggsave("output/Fig_S6_IgM.pdf", plotMLpanelsCombinedIgM,
       width = 10, height = 8, dpi = 600)

ggsave("output/Fig_S6_IgM.jpg", plotMLpanelsCombinedIgM,
       width =10, height = 8, dpi = 600)


##
# Plot timeline in acute panel
plotMPoxTime <-
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype != "IgA") %>% 
  filter(!is.na(serostatus.delta)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  mutate(time_class = case_when(time_passed == 0 ~ 0,
                                time_passed <= 7 ~ 1,
                                time_passed <= 14 ~ 2,
                                time_passed > 14 ~ 3,
                                panel_strat == "SPox" & MPX_diagnose_final == "Yes" ~ 4)) %>% 
  filter(panel_detail == "MPXV") %>% 
  #  mutate(time_passed = if_else(time_passed==0, 0, log10(time_passed))) %>% 
  ggplot(mapping = aes(x =time_class , y = dataIn)) +
  geom_smooth(fullrange = TRUE) +
  geom_jitter(width = 0.2, alpha = 0.2, aes(color = as.factor(serostatus.delta))) +
  scale_fill_manual(name = "Isotype", values = colorblind_pal()(8)[c(4,5,6,7)]) + 
  scale_color_manual(name = "Serostatus", values = colorblind_pal()(8)[c(2,3)], labels = c("Neg", "Pos")) + 
  scale_x_continuous(name = "Time since first diagnosis", limits = c(0,3), breaks = c(0,1,2,3),labels = c("0", "7", "14", "33")) +
  scale_y_continuous(name = "Analyte (log10 µg/mL VIG)", limits = c(0,6)) +
  facet_grid(isotype~analyte) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size = rel(0.5), angle = 45, vjust = 1, hjust=1))



##
# New plot revised fig: seroconversion for each antigen based on bar plot for IgG
# and IgM respectively
# Assuming your data is in a DataFrame named 'data'
# Pre-processing steps (similar to what was done in Python)
# Filter out 'VACV' analyte
# Convert 'serostatus_cat' to binary
# Create extended bins for 'time_passed'


# Combined -> plot bar plots above each other
# Repeat for IgM
data_processed<- dataInput %>%
  filter(panel == "MPXV") %>% 
  filter(analyte != "VACV") %>%
  filter(isotype != "IgA") %>% 
  mutate(
    serostatus_binary = ifelse(serostatus_cat == "positive", 1, 0),
    time_passed_bins = cut(time_passed, breaks = c(0, 1, 7, 11, 16, 21, 35), include.lowest = TRUE, right = TRUE)
  )


# Calculate seropositivity rate and uncertainty IgM
seropositivity <- data_processed %>%
  group_by(analyte, time_passed_bins, isotype) %>%
  summarize(
    mean = mean(serostatus_binary),
    count = n(),
    std_error = sqrt(mean * (1 - mean) / count)
  )

# Plotting using ggplot
plotMPOXacute <-
  ggplot(seropositivity, aes(x = time_passed_bins, y = mean, fill = isotype)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean - std_error, ymax = mean + std_error), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~analyte, ncol = 4) +
  labs(x = "Time Passed Bins", y = "Seropositivity Rate") +
  scale_x_discrete(name = "Time since first sample (days)", labels = c("0", "< 7", "< 11",
                                                                       "< 16", "< 21", "< 35")) +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colorblind_pal()(8)[c(6,5)]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# plot combined seroconversion


# Plot for MVA
data_processed_MVA  <- dataInput %>%
  filter(!is.na(panel_strat)) %>% 
  filter(isotype != "IgA") %>% 
  filter(!is.na(serostatus.delta)) %>% 
  filter(panel != "MPXV") %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat == "Acute") %>% 
  filter(!is.na(serostatus.meta)) %>% 
  filter(serostatus_factor %in% c("Pre", "Prime", "First Boost")) %>% 
  mutate(serostatus_factor = factor(serostatus_factor, levels = c("Pre", "Prime", "First Boost"),
                                    labels = c("Pre", "Prime", "Boost"),
                                    ordered = TRUE)) %>% 
  mutate(
    serostatus_binary = ifelse(serostatus_cat == "positive", 1, 0)
  )




# Calculate seropositivity rate and uncertainty IgM
seropositivity_MVA <- data_processed_MVA %>%
  group_by(analyte, serostatus_factor, isotype) %>%
  summarize(
    mean = mean(serostatus_binary),
    count = n(),
    std_error = sqrt(mean * (1 - mean) / count)
  ) 

# Plotting using ggplot
plotMVAacute <-
  ggplot(seropositivity_MVA, aes(x = serostatus_factor, y = mean, fill = isotype)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean - std_error, ymax = mean + std_error), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~analyte, ncol = 4) +
  labs(y = "Seropositivity Rate") +
  scale_x_discrete(name = "Timepoint immunisation") +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colorblind_pal()(8)[c(6,5)]) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



##
# Combine plots actue phase
plotAcuteCombined <- 
  ggarrange(plotMPOXacute, plotMVAacute, ncol = 2, align = "hv", common.legend = TRUE)






plotMVAtime <- 
  dataInput %>% 
  filter(!is.na(panel_strat)) %>% 
  filter(isotype != "IgA") %>% 
  filter(!is.na(serostatus.delta)) %>% 
  filter(panel_detail != "SPox") %>% 
  filter(analyte != "VACV") %>% 
  filter(panel_strat == "Acute") %>% 
  # filter(!(serostatus_factor %in% c("Pre Second Boost", "Second Boost"))) %>% 
  #  mutate(analyte = factor(analyte, levels = c("A27L", "A29",
  #                                              "D8L", "E8",
  #                                              "H3L", "L1R", "M1",
  #                                              "A33R", "A35R",
  #                                              "B5R", "B6",
  #                                              "A5L", "ATI-C", "ATI-N", "Delta"))) %>%
  mutate(isotype = factor(isotype, levels = c("IgG", "IgM"))) %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(1,0), labels = c("< 1975", ">= 1975"),
                                ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = serostatus.meta, y = (dataIn))) +
  geom_smooth() +
  geom_jitter(width = 0.2, alpha = 0.2, aes(color = as.factor(serostatus.delta))) +
  scale_fill_manual(name = "Isotype", values = colorblind_pal()(8)[c(4,5,6,7)]) + 
  scale_color_manual(name = "Serostatus", values = colorblind_pal()(8)[c(2,3)], labels = c("Neg", "Pos")) + 
  scale_x_continuous(name = "Time since first diagnosis", breaks = c(1:3),labels = c("Pre", "Prime", "Boost")) +
  theme_pubr()+
  scale_y_continuous(name = "Analyte (log10 µg/mL VIG)", limits = c(0,6)) +
  # coord_cartesian(ylim = c(0.5,5.5)) +
  facet_grid(isotype ~ analyte) +
  theme(plot.subtitle=element_text(hjust=0.5),
        axis.text.x = element_text(size = rel(0.5), angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = "none")







plotFig4Combi <-
  ggarrange(plotAcuteCombined, plotFig4, labels = c("a", ""),
            nrow = 2)

ggsave("output/Fig_4_combined.pdf", plotFig4Combi,
       width = 14, height = 18, dpi = 600)

ggsave("output/Fig_4_combined.jpg", plotFig4Combi,
       width = 14, height = 18, dpi = 600)












seroconversionMVA <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "MVA") %>% 
  filter(analyte != "VACV") %>%  
  filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, childhoodImmu,
           serostatus.meta, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1))



plotHeatMapMVA <-
  seroconversionMVA %>%
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("Non immunized", "Childhood immunized")),
         serostatus.meta = factor(serostatus.meta, levels = c(0,1,2,3),
                                  labels = c("Pre", "Prime", "Boost 1", "Boost 2")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(serostatus.meta), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(childhoodImmu~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")

seroconversionMPXV <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "MPXV") %>% 
  filter(analyte != "VACV") %>%  
  filter(!is.na(childhoodImmu)) %>% 
  mutate(time_factor = case_when(time_passed == 0 ~ "0",
                                 time_passed <= 7~ "<= 7",
                                 time_passed <= 14~ "<= 14",
                                 TRUE ~ "> 14")) %>% 
  group_by(isotype, analyte, childhoodImmu,
           time_factor, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1)) 

plotHeatMapMPXV <-
  seroconversionMPXV %>%
  # filter(isotype == "IgG") %>% 
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("Non immunized", "Childhood immunized")),
         time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor( time_factor), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  # scale_fill_continuous_tableau(name = "Seroconversion rate", na.value = 'white') +
  # scale_fill_gradient(
  #  name = "Fraction seropositive",
  #  high = colorblind_pal()(8)[6],
  #  low = "white",
  #  space = "Lab",
  #  na.value = "grey50",
  #  guide = "colourbar",
  #  aesthetics = "fill") +
  theme_bw() +
  # coord_flip() +
  facet_grid(childhoodImmu~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")

unique(dataInput$panel)

seroconversionNK <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "Pre_New") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, childhoodImmu,serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "NC") 


plotHeatMapNK <-
  seroconversionNK %>%
  # filter(isotype == "IgG") %>% 
  mutate(childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("Non immunized", "Childhood immunized")),
         #  time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  # scale_fill_continuous_tableau(name = "Seroconversion rate", na.value = 'white') +
  # scale_fill_gradient(
  #  name = "Fraction seropositive",
  #  high = colorblind_pal()(8)[6],
  #  low = "white",
  #  space = "Lab",
  #  na.value = "grey50",
  #  guide = "colourbar",
  #  aesthetics = "fill") +
  theme_bw() +
  # coord_flip() +
  facet_grid(childhoodImmu~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")




seroconversionNKold <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "Pre") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "NC") 


#plotHeatMapNKOld <-
seroconversionNKold %>%
  # filter(isotype == "IgG") %>% 
  mutate(#childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("Non immunized", "Childhood immunized")),
    #  time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
    sum = stat_0 + stat_1,
    text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  # scale_fill_continuous_tableau(name = "Seroconversion rate", na.value = 'white') +
  # scale_fill_gradient(
  #  name = "Fraction seropositive",
  #  high = colorblind_pal()(8)[6],
  #  low = "white",
  #  space = "Lab",
  #  na.value = "grey50",
  #  guide = "colourbar",
  #  aesthetics = "fill") +
  theme_bw() +
  # coord_flip() +
  facet_grid(.~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")


seroconversionCPXV <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "CPXV") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, childhoodImmu,serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         childhoodImmu = 0,
         panel = "CPXV") 

seroconversionNKCPXV <-
  rbind(seroconversionNK, seroconversionCPXV)

plotHeatMapNKCPXV <-
  seroconversionNKCPXV %>%
  # filter(isotype == "IgG") %>% 
  mutate(panel = factor(panel, levels = c("NC","CPXV")),
         childhoodImmu = factor(childhoodImmu, levels = c(0,1), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor( panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Panel", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(childhoodImmu~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")



# Load SPox panel with complete description
dataInputComplete <- import("../15 Generate Dataframe ML/output/dataInputComplete.csv") %>% 
  mutate(serostatus_clear = if_else(serostatus_cat == "positive", 1, 0))

seroconversionSpox <-
  dataInputComplete %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "SPox") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, POX_VACC, MPX_vac_final_comb, MPX_diagnose_final, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "SPox") 


plotHeatMapSpoxNoMVA <-
  seroconversionSpox%>%
  filter(MPX_vac_final_comb == "No") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(POX_VACC~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")


plotHeatMapSpoxMVA <-
  seroconversionSpox%>%
  filter(MPX_vac_final_comb == "Yes") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(POX_VACC~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")

plotHeatMapSpoxMVALeg <-
  seroconversionSpox%>%
  filter(MPX_vac_final_comb == "Yes") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(POX_VACC~ isotype) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) 

leg <- get_legend(plotHeatMapSpoxMVALeg)

plotSeroconversionPanels <-
  ggarrange(plotHeatMapNK, plotHeatMapMVA, plotHeatMapMPXV,
            plotHeatMapSpoxNoMVA, plotHeatMapSpoxMVA, leg,
            ncol = 3, nrow = 2,
            labels = c("a", "b", "c", "d", "e", ""))


ggsave("output/plotHeatMapMVA.png", plotHeatMapMVA, width = 7, height = 7, dpi = 300)
ggsave("output/plotHeatMapMPXV.png", plotHeatMapMPXV, width = 7, height = 7, dpi = 300)
ggsave("output/plotHeatMapNKCPXV.png", plotHeatMapNKCPXV, width = 5, height = 5, dpi = 300)
ggsave("output/plotSeroconversionPanels.png", plotSeroconversionPanels, width = 16, height = 10, dpi = 300)



##
# Plot non-stratified
seroconversionMVAIso <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "MVA") %>% 
  filter(analyte != "VACV") %>%  
  filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, 
           serostatus.meta, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1))


plotHeatMapMVAIso <-
  seroconversionMVAIso %>%
  mutate(serostatus.meta = factor(serostatus.meta, levels = c(0,1,2,3),
                                  labels = c("Pre", "Prime", "Boost 1", "Boost 2")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(serostatus.meta), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") +
  ggtitle("Imvanex panel") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

seroconversionMPXVIso <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "MPXV") %>% 
  filter(analyte != "VACV") %>%  
  filter(!is.na(childhoodImmu)) %>% 
  mutate(time_factor = case_when(time_passed == 0 ~ "0",
                                 time_passed <= 7~ "<= 7",
                                 time_passed <= 14~ "<= 14",
                                 TRUE ~ "> 14")) %>% 
  group_by(isotype, analyte,
           time_factor, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1)) 

plotHeatMapMPXVIso <-
  seroconversionMPXVIso %>%
  mutate(time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor( time_factor), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") +
  ggtitle("Acute panel") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


seroconversionNKIso <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "Pre_New") %>% 
  filter(analyte != "VACV") %>%  
  group_by(isotype, analyte,serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "NC") 


plotHeatMapNKIso <-
  seroconversionNKIso %>%
  mutate(sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Time point", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank())   + 
  ggtitle("Negative panel") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

seroconversionSpoxIso <-
  dataInputComplete %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "SPox") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, MPX_vac_final_comb, MPX_diagnose_final, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "SPox") 


plotHeatMapSpoxNoMVAIso <-
  seroconversionSpoxIso %>%
  filter(MPX_vac_final_comb == "No") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         #  POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  ggtitle("Epidemiological panel", subtitle = "Non-immunized") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotHeatMapSpoxMVAIso <-
  seroconversionSpoxIso %>%
  filter(MPX_vac_final_comb == "Yes") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         # POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  ggtitle("Epidemiological panel", subtitle = "Imvanex immunized") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
##
# Plot non-stratified data based on more complete dataset
seroconversionSpoxIsoAll <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% "SPox") %>% 
  filter(analyte != "VACV") %>%  
  # filter(!is.na(childhoodImmu)) %>% 
  group_by(isotype, analyte, MPX_vac_final, MPX_diagnose_final, serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "SPox") 


plotHeatMapSpoxNoMVAIsoAll <-
  seroconversionSpoxIsoAll %>%
  filter(MPX_vac_final == "No") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         #  POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  ggtitle("Epidemiological panel", subtitle = "Non-immunized") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plotHeatMapSpoxMVAIsoAll <-
  seroconversionSpoxIsoAll %>%
  filter(MPX_vac_final == "Yes") %>% 
  mutate(MPX_diagnose_final = factor(MPX_diagnose_final, levels = c("No","Yes")),
         # POX_VACC = factor(POX_VACC, levels = c("No", "Yes"), labels = c("Non immunized", "Childhood immunized")),
         #time_factor = factor(time_factor, levels = c("0","<= 7","<= 14", "> 14")),
         sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(MPX_diagnose_final), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Mpox diagnosed", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  # coord_flip() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  ggtitle("Epidemiological panel", subtitle = "Imvanex immunized") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

##
# Plot seroconversion based on both panels (non-stratified)
seroconversionNKIsoAll <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% c("Pre_New", "Pre")) %>% 
  filter(analyte != "VACV") %>%  
  group_by(isotype, analyte,serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "NC") 


plotHeatMapNKIsoAll <-
  seroconversionNKIsoAll %>%
  mutate(sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Panel", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank())   + 
  ggtitle("Negative panel") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



seroconversionCPXVIsoAll <-
  dataInput %>% 
  filter(isotype != "IgA") %>%  
  filter(panel %in% c("CPXV")) %>% 
  filter(analyte != "VACV") %>%  
  group_by(isotype, analyte,serostatus_clear) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = serostatus_clear, names_prefix = "stat_", values_from = n) %>% 
  mutate(stat_0 = if_else(is.na(stat_0), as.integer(0), stat_0),
         stat_1 = if_else(is.na(stat_1), as.integer(0), stat_1),
         fraq_pos = stat_1/(stat_0 + stat_1),
         panel = "CPXV") 


plotHeatMapCPXVIsoAll <-
  seroconversionCPXVIsoAll %>%
  mutate(sum = stat_0 + stat_1,
         text = paste(stat_1, sum, sep="/")) %>% 
  ggplot(mapping = aes(x = as.factor(panel), y = analyte, fill = fraq_pos)) +
  scale_x_discrete(name = "Panel", expand = c(0,0)) +
  scale_y_discrete(name = "Analyte", expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(name = "Fraction seropositive", limits = c(0, 1))+
  theme_bw() +
  facet_grid(isotype ~ .) +
  geom_tile() +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank())   + 
  ggtitle("CPXV infected") +
  geom_text(aes(label=text), color = "grey50", size = 3) + 
  theme(strip.background = element_blank()) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



plotNonStratified <-
  ggarrange(plotHeatMapNKIsoAll, plotHeatMapMVAIso, plotHeatMapCPXVIsoAll, plotHeatMapMPXVIso, plotHeatMapSpoxMVAIsoAll, plotHeatMapSpoxNoMVAIsoAll,
            ncol = 6, align = "hv", common.legend = TRUE)

ggsave("output/plotNonStratified.png", plotNonStratified, width = 16, height = 7, dpi = 300)


